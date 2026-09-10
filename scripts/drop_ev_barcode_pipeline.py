#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DropEV-seq Barcode Extraction and Clustering Pipeline

Core functionality:
1. Quality control and splitting of paired-end FASTQ files using fastp
2. Barcode extraction from R2 reads within a specified window
3. Barcode clustering using starcode
4. Splitting reads into barcode-specific FASTQ files
5. Generation of QC reports and statistics

Requirements:
    - Python >= 3.8
    - fastp (in PATH or specify --fastp-path)
    - starcode (in PATH or specify --starcode-path)
    - Python packages: psutil, openpyxl

Author: DropEV-seq Platform
License: MIT
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import multiprocessing as mp
import os
import re
import shutil
import subprocess
import sys
import traceback
from collections import Counter
from dataclasses import dataclass
from datetime import datetime
from itertools import zip_longest
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

try:
    import psutil
except ImportError as exc:
    raise ImportError("ERROR: psutil package is required. Install via: pip install psutil") from exc

try:
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font
except ImportError as exc:
    raise ImportError("ERROR: openpyxl package is required. Install via: pip install openpyxl") from exc


# === Default Configuration ===
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_ROOT = PROJECT_ROOT / "output"
LOG_ROOT = PROJECT_ROOT / "logs"
README_LOG = PROJECT_ROOT / "README_Log.md"
EXPECTED_CONDA_ENV = "evdna"
TASK_NAME = "drop_ev_barcode_pipeline"
SUMMARY_THRESHOLDS = [500, 1000, 2000, 3000, 5000]
DEFAULT_FASTP_THREADS = 12
DEFAULT_FASTQC_THREADS = 4
DEFAULT_EXTRACTION_WORKERS = 60
DEFAULT_STARCODE_THREADS = 60


def configure_project_paths(custom_project_root: Path | None) -> None:
    """
    Configure project output and log directories based on user-specified root.

    This allows running the pipeline with outputs directed to a different location
    without copying the script itself.
    """
    global PROJECT_ROOT
    global OUTPUT_ROOT
    global LOG_ROOT
    global README_LOG

    if custom_project_root is None:
        return

    PROJECT_ROOT = custom_project_root.resolve()
    OUTPUT_ROOT = PROJECT_ROOT / "output"
    LOG_ROOT = PROJECT_ROOT / "logs"
    README_LOG = PROJECT_ROOT / "README_Log.md"


@dataclass
class ResourcePlan:
    """Resource allocation plan for parallel execution."""
    cpu_count: int
    cpu_usage_percent: float
    total_memory_gb: float
    available_memory_gb: float
    available_memory_ratio: float
    fastp_threads: int
    extraction_workers: int
    starcode_threads: int
    cluster_batch_size: int


@dataclass
class OutputPaths:
    """Container for all output paths to avoid repeated string concatenation."""
    output_dir: Path
    clean_dir: Path
    qc_dir: Path
    fastp_dir: Path
    fastqc_dir: Path
    intermediate_dir: Path
    split_extract_dir: Path
    extracted_dir: Path
    reports_dir: Path
    logs_dir: Path
    root_log_path: Path
    pipeline_log_path: Path
    config_json_path: Path
    summary_tsv_path: Path
    summary_json_path: Path
    sequencing_stats_xlsx_path: Path
    extraction_split_tsv_path: Path
    barcode_length_tsv_path: Path
    cluster_threshold_tsv_path: Path
    cluster_manifest_tsv_path: Path


@dataclass
class SplitFastqPair:
    """Pair of split FASTQ files from fastp."""
    prefix: str
    read1_path: Path
    read2_path: Path


@dataclass
class ExtractionSplitResult:
    """Results from barcode extraction on a single split file."""
    prefix: str
    total_pairs: int
    matched_pairs: int
    barcode_length_counter: Dict[int, int]
    barcode_fastq_path: Path
    read1_fastq_path: Path
    read2_fastq_path: Path


@dataclass
class ClusterBatchResult:
    """Results from processing a batch of clusters."""
    batch_label: str
    cluster_count: int
    written_reads: int
    manifest_rows: List[Dict[str, object]]


def print_and_log(message: str, log_handle) -> None:
    """Print message to both terminal and log file with timestamp."""
    timestamp_text = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp_text}] {message}"
    print(line, flush=True)
    log_handle.write(line + "\n")
    log_handle.flush()


def ensure_conda_environment() -> None:
    """Ensure script is running in the expected conda environment."""
    current_conda_env = os.environ.get("CONDA_DEFAULT_ENV", "").strip()
    if current_conda_env != EXPECTED_CONDA_ENV:
        raise EnvironmentError(
            f"ERROR: Current conda environment is '{current_conda_env or 'not detected'}'. "
            f"Please activate the '{EXPECTED_CONDA_ENV}' environment first using: "
            f"`conda activate {EXPECTED_CONDA_ENV}` or "
            f"`conda run -n {EXPECTED_CONDA_ENV} python {Path(__file__).name}`"
        )


def normalize_task_name(task_name: str) -> str:
    """Normalize task name to snake_case for version detection."""
    normalized_name = re.sub(r"[^a-zA-Z0-9]+", "_", task_name.strip().lower()).strip("_")
    if not normalized_name:
        raise ValueError("ERROR: Task name cannot be empty")
    return normalized_name


def detect_next_version(task_name: str) -> int:
    """
    Detect next version number by scanning both output directory and README log.
    This ensures version numbers never regress even if old directories are deleted.
    """
    normalized_task_name = normalize_task_name(task_name)
    observed_versions: List[int] = []

    try:
        if OUTPUT_ROOT.exists():
            for directory_path in OUTPUT_ROOT.glob(f"*_{normalized_task_name}_v*"):
                if not directory_path.is_dir():
                    continue
                version_match = re.search(
                    rf"_{re.escape(normalized_task_name)}_v(\d+)$",
                    directory_path.name,
                )
                if version_match:
                    observed_versions.append(int(version_match.group(1)))
    except Exception as exc:
        print(f"ERROR: Failed to scan output directory for version numbers: {normalized_task_name}")
        raise exc

    try:
        if README_LOG.exists():
            readme_text = README_LOG.read_text(encoding="utf-8")
            observed_versions.extend(
                int(version_text)
                for version_text in re.findall(
                    rf"{re.escape(normalized_task_name)}_v(\d+)",
                    readme_text,
                )
            )
    except Exception as exc:
        print(f"ERROR: Failed to read README log: {README_LOG}")
        raise exc

    return max(observed_versions, default=0) + 1


def prepare_output_paths(task_name: str, cluster_distance: int) -> Tuple[OutputPaths, int]:
    """Create all output directories and paths for this pipeline run."""
    version_number = detect_next_version(task_name)
    timestamp_text = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_name = f"{timestamp_text}_{normalize_task_name(task_name)}_v{version_number}"

    output_dir = OUTPUT_ROOT / output_name
    clean_dir = output_dir / "clean_reads"
    qc_dir = output_dir / "qc"
    fastp_dir = qc_dir / "fastp"
    fastqc_dir = qc_dir / "fastqc"
    intermediate_dir = output_dir / "intermediate"
    split_extract_dir = intermediate_dir / "barcode_split_outputs"
    extracted_dir = output_dir / "extracted_reads"
    reports_dir = output_dir / "reports"
    logs_dir = output_dir / "logs"
    root_log_path = LOG_ROOT / f"{output_name}.log"
    pipeline_log_path = logs_dir / f"{normalize_task_name(task_name)}_pipeline.log"

    output_paths = OutputPaths(
        output_dir=output_dir,
        clean_dir=clean_dir,
        qc_dir=qc_dir,
        fastp_dir=fastp_dir,
        fastqc_dir=fastqc_dir,
        intermediate_dir=intermediate_dir,
        split_extract_dir=split_extract_dir,
        extracted_dir=extracted_dir,
        reports_dir=reports_dir,
        logs_dir=logs_dir,
        root_log_path=root_log_path,
        pipeline_log_path=pipeline_log_path,
        config_json_path=reports_dir / f"{normalize_task_name(task_name)}_config_v{version_number}.json",
        summary_tsv_path=reports_dir / f"{normalize_task_name(task_name)}_summary_v{version_number}.tsv",
        summary_json_path=reports_dir / f"{normalize_task_name(task_name)}_summary_v{version_number}.json",
        sequencing_stats_xlsx_path=reports_dir / f"{normalize_task_name(task_name)}_sequencing_stats_v{version_number}.xlsx",
        extraction_split_tsv_path=reports_dir / f"{normalize_task_name(task_name)}_extraction_split_summary_v{version_number}.tsv",
        barcode_length_tsv_path=reports_dir / f"{normalize_task_name(task_name)}_barcode_length_distribution_v{version_number}.tsv",
        cluster_threshold_tsv_path=reports_dir / f"{normalize_task_name(task_name)}_cluster_threshold_summary_v{version_number}.tsv",
        cluster_manifest_tsv_path=reports_dir / f"{normalize_task_name(task_name)}_cluster_manifest_d{cluster_distance}_v{version_number}.tsv",
    )

    for directory_path in [
        OUTPUT_ROOT,
        LOG_ROOT,
        output_paths.output_dir,
        output_paths.clean_dir,
        output_paths.qc_dir,
        output_paths.fastp_dir,
        output_paths.fastqc_dir,
        output_paths.intermediate_dir,
        output_paths.split_extract_dir,
        output_paths.extracted_dir,
        output_paths.reports_dir,
        output_paths.logs_dir,
    ]:
        directory_path.mkdir(parents=True, exist_ok=True)

    return output_paths, version_number


def detect_resources(
    split_count: int,
    requested_fastp_threads: int | None = None,
    requested_extraction_workers: int | None = None,
    requested_starcode_threads: int | None = None,
) -> ResourcePlan:
    """
    Auto-detect system resources and determine optimal parallelization strategy.

    Strategy:
    1. When CPU usage <= 75% and available memory >= 25%, use more aggressive parallelization
    2. When system is busy, reduce thread counts and batch sizes to avoid overload
    3. User-specified values override auto-detection
    """
    cpu_count = psutil.cpu_count(logical=True) or 1
    cpu_usage_percent = psutil.cpu_percent(interval=1)
    memory_info = psutil.virtual_memory()
    total_memory_gb = memory_info.total / (1024**3)
    available_memory_gb = memory_info.available / (1024**3)
    available_memory_ratio = memory_info.available / memory_info.total

    if cpu_usage_percent <= 75 and available_memory_ratio >= 0.25:
        fastp_threads = min(max(cpu_count // 4, 4), DEFAULT_FASTP_THREADS)
        extraction_workers = min(split_count, max(cpu_count // 2, 4), DEFAULT_EXTRACTION_WORKERS)
        starcode_threads = min(max(cpu_count // 2, 4), DEFAULT_STARCODE_THREADS)
        cluster_batch_size = 400
    else:
        fastp_threads = 4
        extraction_workers = min(split_count, 4)
        starcode_threads = 4
        cluster_batch_size = 200

    if requested_fastp_threads is not None:
        fastp_threads = max(1, requested_fastp_threads)
    if requested_extraction_workers is not None:
        extraction_workers = max(1, requested_extraction_workers)
    if requested_starcode_threads is not None:
        starcode_threads = max(1, requested_starcode_threads)

    return ResourcePlan(
        cpu_count=cpu_count,
        cpu_usage_percent=cpu_usage_percent,
        total_memory_gb=total_memory_gb,
        available_memory_gb=available_memory_gb,
        available_memory_ratio=available_memory_ratio,
        fastp_threads=fastp_threads,
        extraction_workers=extraction_workers,
        starcode_threads=starcode_threads,
        cluster_batch_size=cluster_batch_size,
    )


def resolve_executable(
    executable_name: str,
    fallback_paths: Sequence[Path] | None = None,
    prefer_fallback: bool = False,
) -> str | None:
    """Resolve path to external executable, with optional fallback paths."""
    candidate_paths = list(fallback_paths or [])

    if prefer_fallback:
        for fallback_path in candidate_paths:
            if fallback_path.exists() and os.access(fallback_path, os.X_OK):
                return str(fallback_path)

    executable_path = shutil.which(executable_name)
    if executable_path:
        return executable_path

    for fallback_path in candidate_paths:
        if fallback_path.exists() and os.access(fallback_path, os.X_OK):
            return str(fallback_path)
    return None


def ensure_required_tools(run_fastqc: bool, fastp_path: str | None, starcode_path: str | None, fastqc_path: str | None) -> Dict[str, str | None]:
    """Check availability of required external tools."""
    known_fastp_paths = [
        Path(sys.executable).resolve().parent / "fastp",
    ]
    known_starcode_paths: List[Path] = []
    known_fastqc_paths = [
        Path(sys.executable).resolve().parent / "fastqc",
    ]

    if fastp_path:
        known_fastp_paths.insert(0, Path(fastp_path))
    if starcode_path:
        known_starcode_paths.insert(0, Path(starcode_path))
    if fastqc_path:
        known_fastqc_paths.insert(0, Path(fastqc_path))

    tool_paths: Dict[str, str | None] = {
        "fastp": resolve_executable("fastp", known_fastp_paths, prefer_fallback=bool(fastp_path)),
        "starcode": resolve_executable("starcode", known_starcode_paths, prefer_fallback=bool(starcode_path)),
        "fastqc": resolve_executable("fastqc", known_fastqc_paths, prefer_fallback=bool(fastqc_path)) if run_fastqc else None,
    }
    required_missing = [tool_name for tool_name in ("fastp", "starcode") if tool_paths[tool_name] is None]
    if required_missing:
        raise FileNotFoundError(f"ERROR: Required tools not found: {', '.join(required_missing)}")
    return tool_paths


def infer_sample_id(read1_path: Path) -> str:
    """Infer sample ID from R1 filename when not explicitly provided."""
    filename = read1_path.name
    for suffix in ["_R1.fastq.gz", "_R1.fq.gz", "_R1_001.fastq.gz", "_R1_001.fq.gz"]:
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    raise ValueError(
        "ERROR: Cannot infer sample_id from R1 filename. Please provide --sample-id explicitly."
    )


def build_barcode_output_sample_id(sample_id: str) -> str:
    """
    Build sample ID prefix for barcode-level FASTQ filenames.

    Removes '_merged' suffix if present to avoid including technical processing state
    in barcode identifiers.
    """
    if sample_id.endswith("_merged"):
        return sample_id[: -len("_merged")]
    return sample_id


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="DropEV-seq barcode extraction and clustering pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python drop_ev_barcode_pipeline.py \\
      --read1 sample_R1.fastq.gz \\
      --read2 sample_R2.fastq.gz \\
      --barcode-file barcodes.txt \\
      --sample-id sample001 \\
      --output-dir results/

Output structure:
  results/
    ├── clean_reads/           # Quality-filtered and split reads
    ├── qc/                    # QC reports (fastp, fastqc)
    ├── intermediate/          # Intermediate processing files
    ├── extracted_reads/       # Barcode-split FASTQ files
    ├── reports/               # Summary statistics and manifests
    └── logs/                  # Execution logs
        """
    )
    parser.add_argument("--read1", required=True, help="Path to R1 FASTQ.gz file")
    parser.add_argument("--read2", required=True, help="Path to R2 FASTQ.gz file")
    parser.add_argument(
        "--barcode-file",
        required=True,
        help="Path to barcode whitelist file (one barcode per line)",
    )
    parser.add_argument(
        "--sample-id",
        default=None,
        help="Sample identifier; if not provided, will attempt to infer from R1 filename",
    )
    parser.add_argument(
        "--project-root",
        default=None,
        help="Optional project root directory for output/logs/README",
    )
    parser.add_argument(
        "--existing-clean-dir",
        default=None,
        help="If provided, skip fastp and reuse existing clean split directory",
    )
    parser.add_argument(
        "--barcode-window-start",
        type=int,
        default=13,
        help="R2 barcode search window start position (0-based, inclusive). Default: 13",
    )
    parser.add_argument(
        "--barcode-window-end",
        type=int,
        default=24,
        help="R2 barcode search window end position (0-based, exclusive). Default: 24",
    )
    parser.add_argument(
        "--split-count",
        type=int,
        default=60,
        help="Number of split files for fastp output. Default: 60",
    )
    parser.add_argument(
        "--min-read-length",
        type=int,
        default=120,
        help="Minimum read length to retain after fastp filtering. Default: 120",
    )
    parser.add_argument(
        "--cluster-distance",
        type=int,
        default=0,
        help="Starcode clustering distance parameter (-d). Default: 0",
    )
    parser.add_argument(
        "--min-cluster-reads",
        type=int,
        default=1000,
        help="Minimum reads per cluster (strictly greater than). Default: 1000",
    )
    parser.add_argument(
        "--disable-fastqc",
        action="store_true",
        help="Skip fastqc step",
    )
    parser.add_argument("--fastp-threads", type=int, default=None, help="Override fastp thread count")
    parser.add_argument(
        "--extraction-workers",
        type=int,
        default=None,
        help="Override barcode extraction worker count",
    )
    parser.add_argument(
        "--starcode-threads",
        type=int,
        default=None,
        help="Override starcode thread count",
    )
    parser.add_argument(
        "--fastp-path",
        type=str,
        default=None,
        help="Path to fastp executable (if not in PATH)",
    )
    parser.add_argument(
        "--starcode-path",
        type=str,
        default=None,
        help="Path to starcode executable (if not in PATH)",
    )
    parser.add_argument(
        "--fastqc-path",
        type=str,
        default=None,
        help="Path to fastqc executable (if not in PATH)",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> Tuple[Path, Path, Path, Path | None, Path | None]:
    """Validate user input paths and key parameters."""
    read1_path = Path(args.read1).resolve()
    read2_path = Path(args.read2).resolve()
    barcode_file = Path(args.barcode_file).resolve()
    existing_clean_dir = Path(args.existing_clean_dir).resolve() if args.existing_clean_dir else None
    project_root = Path(args.project_root).resolve() if args.project_root else None

    if not read1_path.exists():
        raise FileNotFoundError(f"ERROR: R1 file not found: {read1_path}")
    if not read2_path.exists():
        raise FileNotFoundError(f"ERROR: R2 file not found: {read2_path}")
    if not barcode_file.exists():
        raise FileNotFoundError(f"ERROR: Barcode file not found: {barcode_file}")
    if existing_clean_dir is not None and not existing_clean_dir.exists():
        raise FileNotFoundError(f"ERROR: Existing clean directory not found: {existing_clean_dir}")
    if project_root is not None and project_root.is_file():
        raise NotADirectoryError(f"ERROR: --project-root cannot point to a file: {project_root}")
    if args.barcode_window_start < 0 or args.barcode_window_end <= args.barcode_window_start:
        raise ValueError("ERROR: Invalid barcode window parameters. Check --barcode-window-start / --barcode-window-end")
    if args.split_count < 1:
        raise ValueError("ERROR: --split-count must be >= 1")
    if args.min_cluster_reads < 0:
        raise ValueError("ERROR: --min-cluster-reads must be >= 0")

    return read1_path, read2_path, barcode_file, existing_clean_dir, project_root


def run_command(command: Sequence[str], log_file: Path) -> None:
    """Execute external command and write stdout/stderr to log file."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with log_file.open("a", encoding="utf-8") as log_handle:
            log_handle.write("[COMMAND] " + " ".join(command) + "\n")
            subprocess.run(
                list(command),
                check=True,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
            )
    except FileNotFoundError as exc:
        print(f"ERROR: Command not found: {' '.join(command)}")
        raise exc
    except subprocess.CalledProcessError as exc:
        print(f"ERROR: Command failed. Check log file: {log_file}")
        raise exc


def read_barcodes(barcode_file: Path) -> List[str]:
    """Read barcode list from file, preserving original order."""
    try:
        barcodes = [line.strip().upper() for line in barcode_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    except Exception as exc:
        print(f"ERROR: Failed to read barcode file: {barcode_file}")
        raise exc

    if not barcodes:
        raise ValueError(f"ERROR: Barcode file is empty: {barcode_file}")
    return barcodes


def write_tsv(path: Path, rows: Iterable[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    """Write TSV table for downstream analysis or manual inspection."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: Path, payload: Dict[str, object]) -> None:
    """Write JSON file for storing configuration and summary data."""
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def write_xlsx(path: Path, rows: Iterable[Dict[str, object]], fieldnames: Sequence[str], sheet_name: str = "sheet1") -> None:
    """
    Write Excel workbook with basic formatting.

    Provides consistent interface with write_tsv:
    1. fieldnames determines column order
    2. Header row is formatted for readability
    """
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = sheet_name[:31] if sheet_name else "sheet1"
    worksheet.freeze_panes = "A2"

    worksheet.append(list(fieldnames))
    header_font = Font(bold=True)
    header_alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    body_alignment = Alignment(vertical="top", wrap_text=True)

    for column_index, fieldname in enumerate(fieldnames, start=1):
        header_cell = worksheet.cell(row=1, column=column_index)
        header_cell.font = header_font
        header_cell.alignment = header_alignment
        worksheet.column_dimensions[header_cell.column_letter].width = min(max(len(str(fieldname)) + 2, 14), 48)

    for row in rows:
        worksheet.append([row.get(fieldname, "") for fieldname in fieldnames])

    for row in worksheet.iter_rows(min_row=2):
        for cell in row:
            cell.alignment = body_alignment

    workbook.save(path)


def run_fastp(
    fastp_executable: str,
    read1_path: Path,
    read2_path: Path,
    sample_id: str,
    output_paths: OutputPaths,
    resource_plan: ResourcePlan,
    min_read_length: int,
    split_count: int,
) -> Tuple[Path, Path, Path]:
    """Run fastp for quality control and read splitting."""
    clean_r1_path = output_paths.clean_dir / f"{sample_id}_R1.clean.fq.gz"
    clean_r2_path = output_paths.clean_dir / f"{sample_id}_R2.clean.fq.gz"
    fastp_json_path = output_paths.fastp_dir / f"{sample_id}.json"
    fastp_html_path = output_paths.fastp_dir / f"{sample_id}.html"
    fastp_log_path = output_paths.logs_dir / f"{sample_id}_fastp.log"

    fastp_command = [
        fastp_executable,
        "-i",
        str(read1_path),
        "-I",
        str(read2_path),
        "-o",
        str(clean_r1_path),
        "-O",
        str(clean_r2_path),
        "--compression",
        "4",
        "--json",
        str(fastp_json_path),
        "--html",
        str(fastp_html_path),
        "-w",
        str(resource_plan.fastp_threads),
        "--length_required",
        str(min_read_length),
        "-s",
        str(split_count),
    ]
    run_command(fastp_command, fastp_log_path)
    return clean_r1_path, clean_r2_path, fastp_json_path


def run_fastqc_if_needed(
    fastqc_executable: str | None,
    read1_path: Path,
    read2_path: Path,
    sample_id: str,
    output_paths: OutputPaths,
    skip_fastqc: bool,
) -> str:
    """Run fastqc if available and not disabled."""
    if skip_fastqc:
        return "disabled_by_user"
    if not fastqc_executable:
        return "fastqc_not_found_skip"

    fastqc_log_path = output_paths.logs_dir / f"{sample_id}_fastqc.log"
    fastqc_command = [
        fastqc_executable,
        "--threads",
        str(DEFAULT_FASTQC_THREADS),
        "--outdir",
        str(output_paths.fastqc_dir),
        str(read1_path),
        str(read2_path),
    ]
    run_command(fastqc_command, fastqc_log_path)
    return "completed"


def collect_split_fastq_pairs(clean_dir: Path, sample_id: str) -> List[SplitFastqPair]:
    """Collect all split FASTQ file pairs from fastp output."""
    r1_pattern = re.compile(rf"^(\d+)\.{re.escape(sample_id)}_R1\.clean\.fq\.gz$")
    r2_pattern = re.compile(rf"^(\d+)\.{re.escape(sample_id)}_R2\.clean\.fq\.gz$")

    r1_map: Dict[str, Path] = {}
    r2_map: Dict[str, Path] = {}

    try:
        for path_item in clean_dir.iterdir():
            if not path_item.is_file():
                continue
            r1_match = r1_pattern.match(path_item.name)
            if r1_match:
                r1_map[r1_match.group(1)] = path_item
                continue
            r2_match = r2_pattern.match(path_item.name)
            if r2_match:
                r2_map[r2_match.group(1)] = path_item
    except Exception as exc:
        print(f"ERROR: Failed to scan clean split files: {clean_dir}")
        raise exc

    shared_prefixes = sorted(set(r1_map).intersection(r2_map), key=lambda item: int(item))
    if not shared_prefixes:
        raise FileNotFoundError(f"ERROR: No paired clean FASTQ files found in: {clean_dir}")

    return [
        SplitFastqPair(prefix=prefix, read1_path=r1_map[prefix], read2_path=r2_map[prefix])
        for prefix in shared_prefixes
    ]


def iterate_fastq_records(fastq_path: Path) -> Iterator[Tuple[str, str, str]]:
    """Stream FASTQ records to avoid loading entire file into memory."""
    opener = gzip.open if fastq_path.suffix == ".gz" else open
    try:
        with opener(fastq_path, "rt", encoding="utf-8", errors="replace") as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                sequence = handle.readline()
                plus_line = handle.readline()
                quality = handle.readline()
                if not sequence or not plus_line or not quality:
                    raise ValueError(f"Incomplete FASTQ record in: {fastq_path}")
                header = header.rstrip("\n\r")
                sequence = sequence.rstrip("\n\r")
                plus_line = plus_line.rstrip("\n\r")
                quality = quality.rstrip("\n\r")
                if not header.startswith("@"):
                    raise ValueError(f"Invalid FASTQ header in: {fastq_path} -> {header}")
                if not plus_line.startswith("+"):
                    raise ValueError(f"Invalid FASTQ plus line in: {fastq_path} -> {plus_line}")
                # Extract read name (first whitespace-delimited field)
                yield header[1:].split()[0], sequence, quality
    except Exception as exc:
        print(f"ERROR: Failed to read FASTQ file: {fastq_path}")
        raise exc


def find_barcode_end(window_sequence: str, ordered_barcodes: Sequence[str]) -> int | None:
    """Find first matching barcode in window and return its end position."""
    for barcode_sequence in ordered_barcodes:
        match_index = window_sequence.find(barcode_sequence)
        if match_index != -1:
            return match_index + len(barcode_sequence)
    return None


def extract_barcodes_from_split(task_payload: Tuple[SplitFastqPair, Sequence[str], int, int, str]) -> ExtractionSplitResult:
    """
    Extract barcodes from a single fastp split file.

    Process:
    1. Search for barcode matches within R2[window_start:window_end]
    2. Output barcode as R2[0:barcode_end] (prefix from start to barcode end)
    3. Keep R1 intact
    4. Output R2_filter as R2[barcode_end:] (remainder after barcode)
    """
    split_pair, ordered_barcodes, window_start, window_end, split_output_root = task_payload
    split_output_dir = Path(split_output_root) / split_pair.prefix
    split_output_dir.mkdir(parents=True, exist_ok=True)

    barcode_fastq_path = split_output_dir / "barcode.fq"
    read1_fastq_path = split_output_dir / "R1_filter.fq"
    read2_fastq_path = split_output_dir / "R2_filter.fq"

    total_pairs = 0
    matched_pairs = 0
    barcode_length_counter: Counter[int] = Counter()

    try:
        with barcode_fastq_path.open("w", encoding="utf-8") as barcode_handle, \
            read1_fastq_path.open("w", encoding="utf-8") as read1_handle, \
            read2_fastq_path.open("w", encoding="utf-8") as read2_handle:
            for read1_record, read2_record in zip_longest(
                iterate_fastq_records(split_pair.read1_path),
                iterate_fastq_records(split_pair.read2_path),
            ):
                if read1_record is None or read2_record is None:
                    raise ValueError(
                        f"Mismatched R1/R2 read counts: {split_pair.read1_path} vs {split_pair.read2_path}"
                    )

                total_pairs += 1
                read1_name, read1_sequence, read1_quality = read1_record
                read2_name, read2_sequence, read2_quality = read2_record

                read2_sequence_upper = read2_sequence.upper()
                window_sequence = read2_sequence_upper[window_start:window_end]
                barcode_relative_end = find_barcode_end(window_sequence, ordered_barcodes)
                if barcode_relative_end is None:
                    continue

                barcode_absolute_end = window_start + barcode_relative_end
                matched_pairs += 1
                barcode_length_counter[barcode_absolute_end] += 1

                barcode_handle.write(
                    f"@{read2_name}\n{read2_sequence[:barcode_absolute_end]}\n+\n{read2_quality[:barcode_absolute_end]}\n"
                )
                read1_handle.write(f"@{read1_name}\n{read1_sequence}\n+\n{read1_quality}\n")
                read2_handle.write(
                    f"@{read2_name}\n{read2_sequence[barcode_absolute_end:]}\n+\n{read2_quality[barcode_absolute_end:]}\n"
                )
    except Exception as exc:
        print(f"ERROR: Barcode extraction failed for split={split_pair.prefix}")
        raise exc

    return ExtractionSplitResult(
        prefix=split_pair.prefix,
        total_pairs=total_pairs,
        matched_pairs=matched_pairs,
        barcode_length_counter=dict(barcode_length_counter),
        barcode_fastq_path=barcode_fastq_path,
        read1_fastq_path=read1_fastq_path,
        read2_fastq_path=read2_fastq_path,
    )


def run_barcode_extraction(
    split_pairs: Sequence[SplitFastqPair],
    ordered_barcodes: Sequence[str],
    output_paths: OutputPaths,
    resource_plan: ResourcePlan,
    window_start: int,
    window_end: int,
    log_handle,
) -> List[ExtractionSplitResult]:
    """Run barcode extraction on all split files in parallel."""
    task_payloads = [
        (split_pair, ordered_barcodes, window_start, window_end, str(output_paths.split_extract_dir))
        for split_pair in split_pairs
    ]

    extraction_results: List[ExtractionSplitResult] = []
    print_and_log(
        f"[EXTRACT] Processing {len(task_payloads)} split files with {resource_plan.extraction_workers} workers",
        log_handle,
    )

    if resource_plan.extraction_workers <= 1:
        for task_payload in task_payloads:
            extraction_results.append(extract_barcodes_from_split(task_payload))
    else:
        with mp.Pool(processes=resource_plan.extraction_workers) as pool:
            for extraction_result in pool.imap_unordered(extract_barcodes_from_split, task_payloads):
                extraction_results.append(extraction_result)

    extraction_results.sort(key=lambda item: int(item.prefix))
    print_and_log(
        f"[EXTRACT] Completed barcode extraction for {len(extraction_results)} splits",
        log_handle,
    )
    return extraction_results


def merge_split_outputs(
    extraction_results: Sequence[ExtractionSplitResult],
    sample_id: str,
    output_paths: OutputPaths,
    log_handle,
) -> Tuple[Path, Path, Path, Dict[int, int], List[Dict[str, object]]]:
    """Merge intermediate files from all splits in order to maintain read synchronization."""
    merged_barcode_fastq_path = output_paths.extracted_dir / f"{sample_id}_S1_L001_bc_001.fastq"
    merged_read1_fastq_path = output_paths.extracted_dir / f"{sample_id}_S1_L001_R1_001.fastq"
    merged_read2_fastq_path = output_paths.extracted_dir / f"{sample_id}_S1_L001_R2_001.fastq"

    aggregate_length_counter: Counter[int] = Counter()
    extraction_summary_rows: List[Dict[str, object]] = []

    try:
        with merged_barcode_fastq_path.open("w", encoding="utf-8") as merged_barcode_handle, \
            merged_read1_fastq_path.open("w", encoding="utf-8") as merged_read1_handle, \
            merged_read2_fastq_path.open("w", encoding="utf-8") as merged_read2_handle:
            for extraction_result in extraction_results:
                extraction_summary_rows.append(
                    {
                        "split_prefix": extraction_result.prefix,
                        "total_pairs": extraction_result.total_pairs,
                        "matched_pairs": extraction_result.matched_pairs,
                        "matched_ratio_percent": round(
                            extraction_result.matched_pairs / extraction_result.total_pairs * 100, 4
                        )
                        if extraction_result.total_pairs
                        else 0.0,
                    }
                )
                aggregate_length_counter.update(extraction_result.barcode_length_counter)

                for source_path, target_handle in [
                    (extraction_result.barcode_fastq_path, merged_barcode_handle),
                    (extraction_result.read1_fastq_path, merged_read1_handle),
                    (extraction_result.read2_fastq_path, merged_read2_handle),
                ]:
                    with source_path.open("r", encoding="utf-8") as source_handle:
                        shutil.copyfileobj(source_handle, target_handle)
    except Exception as exc:
        print("ERROR: Failed to merge split intermediate files")
        raise exc

    print_and_log(
        f"[MERGE] Generated merged barcode/R1/R2 FASTQ files: {merged_barcode_fastq_path.name}",
        log_handle,
    )
    return (
        merged_barcode_fastq_path,
        merged_read1_fastq_path,
        merged_read2_fastq_path,
        dict(aggregate_length_counter),
        extraction_summary_rows,
    )


def load_fastp_summary(fastp_json_path: Path) -> Dict[str, int]:
    """
    Extract read pair counts from fastp JSON output.

    For paired-end data, read pair count = read1_total_reads = read2_total_reads
    This avoids double-counting when using summary.total_reads which includes both R1+R2.
    """
    try:
        fastp_payload = json.loads(fastp_json_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"ERROR: Failed to read fastp JSON: {fastp_json_path}")
        raise exc

    try:
        read1_before_reads = int(fastp_payload["read1_before_filtering"]["total_reads"])
        read1_after_reads = int(fastp_payload["read1_after_filtering"]["total_reads"])
        read2_before_reads = int(fastp_payload["read2_before_filtering"]["total_reads"])
        read2_after_reads = int(fastp_payload["read2_after_filtering"]["total_reads"])
    except Exception as exc:
        try:
            summary_before_reads = int(fastp_payload["summary"]["before_filtering"]["total_reads"])
            summary_after_reads = int(fastp_payload["summary"]["after_filtering"]["total_reads"])
            sequencing_mode = str(fastp_payload.get("summary", {}).get("sequencing", "")).lower()
        except Exception:
            print(f"ERROR: fastp JSON missing required fields: {fastp_json_path}")
            raise exc

        # Fallback for JSONs missing read1/read2 details
        if "paired end" in sequencing_mode:
            if summary_before_reads % 2 != 0 or summary_after_reads % 2 != 0:
                raise ValueError(
                    f"ERROR: fastp summary total_reads is not even, cannot safely divide by 2: {fastp_json_path}"
                )
            return {
                "raw_read_pairs": summary_before_reads // 2,
                "clean_read_pairs": summary_after_reads // 2,
            }
        return {"raw_read_pairs": summary_before_reads, "clean_read_pairs": summary_after_reads}

    if read1_before_reads != read2_before_reads or read1_after_reads != read2_after_reads:
        raise ValueError(
            f"ERROR: Inconsistent read1/read2 total_reads in fastp JSON: {fastp_json_path}"
        )

    return {"raw_read_pairs": read1_before_reads, "clean_read_pairs": read1_after_reads}


def infer_existing_fastp_json(existing_clean_dir: Path, sample_id: str) -> Path | None:
    """
    Attempt to locate fastp JSON from existing clean directory.

    Compatible with directory structure:
    <date>_<sample>/
      ├── 1.Clean.reads/<sample_id>/
      └── 0.fastqc/fastp/<sample_id>.json
    """
    candidate_path = existing_clean_dir.parent.parent / "0.fastqc" / "fastp" / f"{sample_id}.json"
    if candidate_path.exists():
        return candidate_path
    return None


def run_starcode(
    barcode_fastq_path: Path,
    extracted_dir: Path,
    cluster_distance: int,
    resource_plan: ResourcePlan,
    starcode_executable: str,
) -> Path:
    """Run starcode clustering on barcode FASTQ with sequence ID preservation."""
    clr_path = extracted_dir / f"d{cluster_distance}.clr"
    starcode_log_path = extracted_dir.parent / "logs" / f"starcode_d{cluster_distance}.log"

    starcode_command = [
        starcode_executable,
        "-d",
        str(cluster_distance),
        "-t",
        str(resource_plan.starcode_threads),
        "-i",
        str(barcode_fastq_path),
        "-o",
        str(clr_path),
        "-s",
        "--seq-id",
    ]
    run_command(starcode_command, starcode_log_path)
    return clr_path


def summarize_cluster_thresholds(clr_path: Path) -> List[Dict[str, int]]:
    """Compute cluster statistics at different read count thresholds."""
    threshold_summary = {
        threshold: {"cluster_count": 0, "reads_in_clusters": 0}
        for threshold in SUMMARY_THRESHOLDS
    }

    try:
        with clr_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                cluster_size = int(parts[1])
                for threshold in SUMMARY_THRESHOLDS:
                    if cluster_size > threshold:
                        threshold_summary[threshold]["cluster_count"] += 1
                        threshold_summary[threshold]["reads_in_clusters"] += cluster_size
    except Exception as exc:
        print(f"ERROR: Failed to read starcode clr file: {clr_path}")
        raise exc

    return [
        {
            "reads_threshold_gt": threshold,
            "cluster_count": threshold_summary[threshold]["cluster_count"],
            "reads_in_clusters": threshold_summary[threshold]["reads_in_clusters"],
        }
        for threshold in SUMMARY_THRESHOLDS
    ]


def filter_clr_by_cluster_size(clr_path: Path, min_cluster_reads: int) -> Tuple[Path, int, int]:
    """Filter clusters to retain only those with read count strictly greater than threshold."""
    filtered_path = clr_path.with_name(f"{clr_path.name}_filtered")
    total_cluster_count = 0
    filtered_cluster_count = 0

    try:
        with clr_path.open("r", encoding="utf-8") as source_handle, filtered_path.open("w", encoding="utf-8") as target_handle:
            for line in source_handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                total_cluster_count += 1
                cluster_size = int(parts[1])
                if cluster_size > min_cluster_reads:
                    target_handle.write(line)
                    filtered_cluster_count += 1
    except Exception as exc:
        print(f"ERROR: Failed to filter clr file: {clr_path}")
        raise exc

    return filtered_path, total_cluster_count, filtered_cluster_count


def parse_filtered_clr(filtered_clr_path: Path) -> Tuple[Dict[int, str], Dict[str, int]]:
    """
    Parse filtered cluster file to build index-to-cluster mapping.

    Starcode --seq-id outputs 1-based indices; we convert to 0-based internally.
    Cluster sizes are validated against the declared size in column 2.
    """
    index_to_cluster: Dict[int, str] = {}
    declared_cluster_size_map: Dict[str, int] = {}
    observed_cluster_size_counter: Counter[str] = Counter()

    try:
        with filtered_clr_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                cluster_key = parts[0]
                declared_cluster_size = int(parts[1])
                declared_cluster_size_map[cluster_key] = declared_cluster_size
                for index_text in parts[2].split(","):
                    if not index_text:
                        continue
                    index_to_cluster[int(index_text) - 1] = cluster_key
                    observed_cluster_size_counter[cluster_key] += 1
    except Exception as exc:
        print(f"ERROR: Failed to parse filtered clr file: {filtered_clr_path}")
        raise exc

    cluster_size_map: Dict[str, int] = {}
    for cluster_key, observed_cluster_size in observed_cluster_size_counter.items():
        declared_cluster_size = declared_cluster_size_map.get(cluster_key)
        if declared_cluster_size is None:
            raise ValueError(f"ERROR: Cluster missing declared size: {cluster_key}")
        if declared_cluster_size != observed_cluster_size:
            raise ValueError(
                f"ERROR: Cluster size mismatch for {cluster_key}: "
                f"declared={declared_cluster_size}, observed={observed_cluster_size}"
            )
        cluster_size_map[cluster_key] = observed_cluster_size

    return index_to_cluster, cluster_size_map


def iterate_cluster_sync_records(
    merged_barcode_fastq_path: Path,
    merged_read1_fastq_path: Path,
    merged_read2_fastq_path: Path,
) -> Iterator[Tuple[int, Tuple[str, str, str], Tuple[str, str, str], Tuple[str, str, str]]]:
    """
    Synchronously iterate over barcode/R1/R2 merged FASTQ files.

    This ensures starcode indices are correctly bound to barcode FASTQ record order,
    and R1/R2 splitting advances in lockstep with barcode records.
    """
    for record_index, (barcode_record, read1_record, read2_record) in enumerate(
        zip_longest(
            iterate_fastq_records(merged_barcode_fastq_path),
            iterate_fastq_records(merged_read1_fastq_path),
            iterate_fastq_records(merged_read2_fastq_path),
        )
    ):
        if barcode_record is None or read1_record is None or read2_record is None:
            raise ValueError(
                f"ERROR: Merged barcode/R1/R2 FASTQ files have mismatched record counts: "
                f"{merged_barcode_fastq_path} | {merged_read1_fastq_path} | {merged_read2_fastq_path}"
            )
        yield record_index, barcode_record, read1_record, read2_record


def process_cluster_batch(
    batch_label: str,
    barcode_output_sample_id: str,
    cluster_items: Sequence[Tuple[str, int]],
    index_to_cluster: Dict[int, str],
    merged_barcode_fastq_path: Path,
    merged_read1_fastq_path: Path,
    merged_read2_fastq_path: Path,
    cluster_output_dir: Path,
) -> ClusterBatchResult:
    """
    Process a batch of clusters by splitting reads into barcode-specific FASTQ files.

    Barcode-level filename format: <sample_id>-<cluster_size>-<cluster_key>_S1_L001_R{1,2}_001.fastq.gz
    This format preserves sample origin and cluster information for downstream analysis.
    """
    manifest_rows: List[Dict[str, object]] = []
    written_reads = 0

    if not cluster_items:
        return ClusterBatchResult(
            batch_label=batch_label,
            cluster_count=0,
            written_reads=0,
            manifest_rows=[],
        )

    cluster_output_dir.mkdir(parents=True, exist_ok=True)
    current_batch_keys = {cluster_key for cluster_key, _cluster_size in cluster_items}
    file_handles: Dict[str, Tuple[gzip.GzipFile, gzip.GzipFile]] = {}

    try:
        for cluster_key, cluster_size in cluster_items:
            barcode_file_prefix = f"{barcode_output_sample_id}-{cluster_size}-{cluster_key}"
            read1_output_path = cluster_output_dir / f"{barcode_file_prefix}_S1_L001_R1_001.fastq.gz"
            read2_output_path = cluster_output_dir / f"{barcode_file_prefix}_S1_L001_R2_001.fastq.gz"
            manifest_rows.append(
                {
                    "sample_id": barcode_output_sample_id,
                    "cluster_batch_label": batch_label,
                    "cluster_size": cluster_size,
                    "cluster_key": cluster_key,
                    "barcode_file_prefix": barcode_file_prefix,
                    "read1_fastq_gz": str(read1_output_path),
                    "read2_fastq_gz": str(read2_output_path),
                }
            )
            file_handles[cluster_key] = (
                gzip.open(read1_output_path, "wt", compresslevel=6, encoding="utf-8"),
                gzip.open(read2_output_path, "wt", compresslevel=6, encoding="utf-8"),
            )

        for record_index, _barcode_record, read1_record, read2_record in iterate_cluster_sync_records(
            merged_barcode_fastq_path=merged_barcode_fastq_path,
            merged_read1_fastq_path=merged_read1_fastq_path,
            merged_read2_fastq_path=merged_read2_fastq_path,
        ):
            cluster_key = index_to_cluster.get(record_index)
            if cluster_key is None or cluster_key not in current_batch_keys:
                continue

            read1_handle, read2_handle = file_handles[cluster_key]
            read1_name, read1_sequence, read1_quality = read1_record
            read2_name, read2_sequence, read2_quality = read2_record
            read1_handle.write(f"@{read1_name}\n{read1_sequence}\n+\n{read1_quality}\n")
            read2_handle.write(f"@{read2_name}\n{read2_sequence}\n+\n{read2_quality}\n")
            written_reads += 1
    except Exception as exc:
        print(f"ERROR: Failed to split FASTQ by cluster batch: {batch_label}")
        raise exc
    finally:
        for read1_handle, read2_handle in file_handles.values():
            read1_handle.close()
            read2_handle.close()

    return ClusterBatchResult(
        batch_label=batch_label,
        cluster_count=len(cluster_items),
        written_reads=written_reads,
        manifest_rows=manifest_rows,
    )


def run_cluster_split(
    filtered_clr_path: Path,
    barcode_output_sample_id: str,
    merged_barcode_fastq_path: Path,
    merged_read1_fastq_path: Path,
    merged_read2_fastq_path: Path,
    cluster_distance: int,
    output_paths: OutputPaths,
    resource_plan: ResourcePlan,
    log_handle,
) -> Tuple[Path, List[ClusterBatchResult]]:
    """
    Split reads by cluster assignment into barcode-level FASTQ files.

    Strategy:
    1. Load entire filtered cluster file
    2. Sort clusters by size (descending)
    3. Process in batches to balance memory and I/O
    4. Each batch synchronously scans barcode/R1/R2 FASTQ files
    """
    cluster_output_dir = output_paths.extracted_dir / f"d{cluster_distance}"
    cluster_output_dir.mkdir(parents=True, exist_ok=True)

    index_to_cluster, cluster_size_map = parse_filtered_clr(filtered_clr_path)
    cluster_items = sorted(cluster_size_map.items(), key=lambda item: item[1], reverse=True)
    total_batches = (len(cluster_items) + resource_plan.cluster_batch_size - 1) // resource_plan.cluster_batch_size

    print_and_log(
        f"[SPLIT] Processing clusters in batches (batch_size={resource_plan.cluster_batch_size}, "
        f"total_batches={total_batches})",
        log_handle,
    )

    cluster_batch_results: List[ClusterBatchResult] = []
    batch_label_width = max(2, len(str(max(total_batches, 1))))

    for batch_start in range(0, len(cluster_items), resource_plan.cluster_batch_size):
        batch_index = batch_start // resource_plan.cluster_batch_size + 1
        current_batch_items = cluster_items[batch_start : batch_start + resource_plan.cluster_batch_size]
        batch_label = f"batch_{batch_index:0{batch_label_width}d}"
        print_and_log(
            f"[SPLIT] Processing cluster batch {batch_index}/{total_batches}: {batch_label}",
            log_handle,
        )
        cluster_batch_result = process_cluster_batch(
            batch_label=batch_label,
            barcode_output_sample_id=barcode_output_sample_id,
            cluster_items=current_batch_items,
            index_to_cluster=index_to_cluster,
            merged_barcode_fastq_path=merged_barcode_fastq_path,
            merged_read1_fastq_path=merged_read1_fastq_path,
            merged_read2_fastq_path=merged_read2_fastq_path,
            cluster_output_dir=cluster_output_dir,
        )
        cluster_batch_results.append(cluster_batch_result)
        print_and_log(
            f"[SPLIT] Completed {batch_label}: {cluster_batch_result.cluster_count} clusters, "
            f"{cluster_batch_result.written_reads} reads written",
            log_handle,
        )

    return cluster_output_dir, cluster_batch_results


def build_pipeline_summary(
    sample_id: str,
    barcode_file: Path,
    barcodes: Sequence[str],
    fastp_counts: Dict[str, int],
    extraction_results: Sequence[ExtractionSplitResult],
    barcode_length_counter: Dict[int, int],
    threshold_rows: Sequence[Dict[str, int]],
    total_cluster_count: int,
    filtered_cluster_count: int,
    cluster_output_dir: Path,
    args: argparse.Namespace,
    resource_plan: ResourcePlan,
    fastqc_status: str,
) -> Dict[str, object]:
    """Build comprehensive pipeline summary with all key statistics."""
    matched_pairs = sum(item.matched_pairs for item in extraction_results)
    total_pairs = sum(item.total_pairs for item in extraction_results)
    generated_cluster_fastq_count = len(list(cluster_output_dir.glob("*_R1_001.fastq.gz")))
    raw_read_pairs = int(fastp_counts.get("raw_read_pairs", 0) or 0)
    clean_read_pairs = int(fastp_counts.get("clean_read_pairs", 0) or 0)

    summary_payload: Dict[str, object] = {
        "sample_id": sample_id,
        "read1_path": str(args.read1),
        "read2_path": str(args.read2),
        "barcode_file": str(barcode_file),
        "barcode_count": len(barcodes),
        "barcode_window_start_0based": args.barcode_window_start,
        "barcode_window_end_0based": args.barcode_window_end,
        "barcode_window_human_readable": f"R2 positions {args.barcode_window_start + 1}-{args.barcode_window_end} (1-based)",
        "raw_read_pairs": raw_read_pairs,
        "clean_read_pairs": clean_read_pairs,
        "clean_read_retention_percent": round(
            clean_read_pairs / raw_read_pairs * 100, 4
        )
        if raw_read_pairs
        else 0.0,
        "barcode_matched_pairs": matched_pairs,
        "barcode_match_percent_of_clean": round(
            matched_pairs / clean_read_pairs * 100, 4
        )
        if clean_read_pairs
        else 0.0,
        "processed_clean_read_pairs": total_pairs,
        "barcode_match_percent_of_processed": round(matched_pairs / total_pairs * 100, 4) if total_pairs else 0.0,
        "observed_barcode_prefix_lengths": ",".join(str(item) for item in sorted(barcode_length_counter)),
        "starcode_distance": args.cluster_distance,
        "min_cluster_reads_strict_gt": args.min_cluster_reads,
        "total_clusters_from_starcode": total_cluster_count,
        "filtered_clusters_after_threshold": filtered_cluster_count,
        "generated_barcode_fastq_pairs": generated_cluster_fastq_count,
        "split_count": args.split_count,
        "fastqc_status": fastqc_status,
        "cpu_count": resource_plan.cpu_count,
        "cpu_usage_percent": round(resource_plan.cpu_usage_percent, 2),
        "available_memory_gb": round(resource_plan.available_memory_gb, 2),
        "available_memory_ratio_percent": round(resource_plan.available_memory_ratio * 100, 2),
        "fastp_threads": resource_plan.fastp_threads,
        "extraction_workers": resource_plan.extraction_workers,
        "starcode_threads": resource_plan.starcode_threads,
        "cluster_batch_size": resource_plan.cluster_batch_size,
        "fastp_split_strategy": "fastp_native_s_split",
        "cluster_split_strategy": "whole_filtered_clr_batched_barcode_r1_r2_sync_scan",
        "cluster_output_dir": str(cluster_output_dir),
    }

    for threshold_row in threshold_rows:
        threshold_value = threshold_row["reads_threshold_gt"]
        summary_payload[f"clusters_gt_{threshold_value}"] = threshold_row["cluster_count"]
        summary_payload[f"reads_in_clusters_gt_{threshold_value}"] = threshold_row["reads_in_clusters"]

    return summary_payload


def build_sequencing_stats_row(summary_payload: Dict[str, object]) -> Dict[str, object]:
    """
    Build sequencing statistics table row compatible with standard summary format.

    Column names include descriptions for clarity when viewing directly.
    """
    threshold_cluster_count_map = {
        int(summary_payload_key.split("_")[-1]): summary_payload[summary_payload_key]
        for summary_payload_key in summary_payload
        if summary_payload_key.startswith("clusters_gt_")
    }
    threshold_cluster_reads_map = {
        int(summary_payload_key.split("_")[-1]): summary_payload[summary_payload_key]
        for summary_payload_key in summary_payload
        if summary_payload_key.startswith("reads_in_clusters_gt_")
    }

    stats_row: Dict[str, object] = {
        "sample_id": summary_payload["sample_id"],
        "raw_read_num (original read pairs before fastp)": summary_payload["raw_read_pairs"],
        "clean_read_num (read pairs retained after fastp)": summary_payload["clean_read_pairs"],
        "clean_reads_percentage (retention rate)": summary_payload["clean_read_retention_percent"],
        "barcode_num (read pairs with extracted barcode)": summary_payload["barcode_matched_pairs"],
        "barcode_percentage (barcode extraction rate from clean)": summary_payload["barcode_match_percent_of_clean"],
    }

    for threshold in SUMMARY_THRESHOLDS:
        stats_row[
            f"barcode{threshold}_num (clusters with >{threshold} reads)"
        ] = threshold_cluster_count_map.get(threshold, 0)
    for threshold in SUMMARY_THRESHOLDS:
        stats_row[
            f"barcode{threshold}_num_reads (total reads in clusters >{threshold})"
        ] = threshold_cluster_reads_map.get(threshold, 0)
    return stats_row


def main() -> int:
    """Main pipeline entry point."""
    args = parse_args()

    try:
        ensure_conda_environment()
        read1_path, read2_path, barcode_file, existing_clean_dir, project_root = validate_args(args)
        configure_project_paths(project_root)
        sample_id = args.sample_id or infer_sample_id(read1_path)
        barcode_output_sample_id = build_barcode_output_sample_id(sample_id)
        output_paths, version_number = prepare_output_paths(TASK_NAME, args.cluster_distance)

        with output_paths.root_log_path.open("w", encoding="utf-8") as log_handle:
            print_and_log(f"[START] Task: {TASK_NAME}, Version: v{version_number}", log_handle)
            print_and_log(f"[START] Sample: {sample_id}", log_handle)
            print_and_log(f"[START] Barcode output sample prefix: {barcode_output_sample_id}", log_handle)
            print_and_log(f"[START] Output directory: {output_paths.output_dir}", log_handle)

            resource_plan = detect_resources(
                split_count=args.split_count,
                requested_fastp_threads=args.fastp_threads,
                requested_extraction_workers=args.extraction_workers,
                requested_starcode_threads=args.starcode_threads,
            )
            print_and_log(
                f"[RESOURCE] CPU={resource_plan.cpu_count}, usage={resource_plan.cpu_usage_percent:.2f}%, "
                f"total_mem={resource_plan.total_memory_gb:.2f}GB, avail_mem={resource_plan.available_memory_gb:.2f}GB, "
                f"fastp_threads={resource_plan.fastp_threads}, extraction_workers={resource_plan.extraction_workers}, "
                f"starcode_threads={resource_plan.starcode_threads}, cluster_batch_size={resource_plan.cluster_batch_size}",
                log_handle,
            )

            tool_paths = ensure_required_tools(
                run_fastqc=not args.disable_fastqc,
                fastp_path=args.fastp_path,
                starcode_path=args.starcode_path,
                fastqc_path=args.fastqc_path,
            )
            print_and_log(
                f"[TOOLS] fastp={tool_paths['fastp']} | starcode={tool_paths['starcode']} | fastqc={tool_paths['fastqc'] or 'disabled/not found'}",
                log_handle,
            )

            barcodes = read_barcodes(barcode_file)
            print_and_log(
                f"[BARCODE] Loaded {len(barcodes)} barcodes from whitelist",
                log_handle,
            )

            fastqc_status = run_fastqc_if_needed(
                fastqc_executable=tool_paths["fastqc"],
                read1_path=read1_path,
                read2_path=read2_path,
                sample_id=sample_id,
                output_paths=output_paths,
                skip_fastqc=args.disable_fastqc,
            )
            print_and_log(f"[FASTQC] Status: {fastqc_status}", log_handle)

            if existing_clean_dir is None:
                _clean_r1_path, _clean_r2_path, fastp_json_path = run_fastp(
                    fastp_executable=tool_paths["fastp"],
                    read1_path=read1_path,
                    read2_path=read2_path,
                    sample_id=sample_id,
                    output_paths=output_paths,
                    resource_plan=resource_plan,
                    min_read_length=args.min_read_length,
                    split_count=args.split_count,
                )
                fastp_counts = load_fastp_summary(fastp_json_path)
                split_source_dir = output_paths.clean_dir
                print_and_log("[FASTP] Quality control and splitting completed", log_handle)
            else:
                split_source_dir = existing_clean_dir
                existing_fastp_json = infer_existing_fastp_json(existing_clean_dir, sample_id)
                if existing_fastp_json is not None:
                    fastp_counts = load_fastp_summary(existing_fastp_json)
                    print_and_log(
                        f"[FASTP] Skipped, reusing existing clean split directory and fastp JSON: {existing_fastp_json}",
                        log_handle,
                    )
                else:
                    fastp_counts = {"raw_read_pairs": 0, "clean_read_pairs": 0}
                    print_and_log(
                        f"[FASTP] Skipped, reusing existing clean split directory (no fastp JSON found): {existing_clean_dir}",
                        log_handle,
                    )

            split_pairs = collect_split_fastq_pairs(split_source_dir, sample_id)
            print_and_log(f"[FASTP] Detected {len(split_pairs)} split file pairs", log_handle)

            extraction_results = run_barcode_extraction(
                split_pairs=split_pairs,
                ordered_barcodes=barcodes,
                output_paths=output_paths,
                resource_plan=resource_plan,
                window_start=args.barcode_window_start,
                window_end=args.barcode_window_end,
                log_handle=log_handle,
            )

            (
                merged_barcode_fastq_path,
                merged_read1_fastq_path,
                merged_read2_fastq_path,
                barcode_length_counter,
                extraction_summary_rows,
            ) = merge_split_outputs(
                extraction_results=extraction_results,
                sample_id=sample_id,
                output_paths=output_paths,
                log_handle=log_handle,
            )

            if int(fastp_counts.get("clean_read_pairs", 0) or 0) == 0:
                fastp_counts["clean_read_pairs"] = sum(item.total_pairs for item in extraction_results)

            clr_path = run_starcode(
                barcode_fastq_path=merged_barcode_fastq_path,
                extracted_dir=output_paths.extracted_dir,
                cluster_distance=args.cluster_distance,
                resource_plan=resource_plan,
                starcode_executable=str(tool_paths["starcode"]),
            )
            print_and_log(f"[STARCODE] Clustering completed: {clr_path.name}", log_handle)

            threshold_rows = summarize_cluster_thresholds(clr_path)
            filtered_clr_path, total_cluster_count, filtered_cluster_count = filter_clr_by_cluster_size(
                clr_path=clr_path,
                min_cluster_reads=args.min_cluster_reads,
            )
            print_and_log(
                f"[CLUSTER] Total clusters: {total_cluster_count}, filtered clusters: {filtered_cluster_count} (threshold >{args.min_cluster_reads})",
                log_handle,
            )

            cluster_output_dir, cluster_batch_results = run_cluster_split(
                filtered_clr_path=filtered_clr_path,
                barcode_output_sample_id=barcode_output_sample_id,
                merged_barcode_fastq_path=merged_barcode_fastq_path,
                merged_read1_fastq_path=merged_read1_fastq_path,
                merged_read2_fastq_path=merged_read2_fastq_path,
                cluster_distance=args.cluster_distance,
                output_paths=output_paths,
                resource_plan=resource_plan,
                log_handle=log_handle,
            )

            barcode_length_rows = [
                {
                    "barcode_prefix_length": barcode_length,
                    "read_pairs": barcode_length_counter[barcode_length],
                }
                for barcode_length in sorted(barcode_length_counter)
            ]
            cluster_manifest_rows = [
                row
                for cluster_batch_result in cluster_batch_results
                for row in cluster_batch_result.manifest_rows
            ]

            write_tsv(
                output_paths.extraction_split_tsv_path,
                extraction_summary_rows,
                ["split_prefix", "total_pairs", "matched_pairs", "matched_ratio_percent"],
            )
            write_tsv(
                output_paths.barcode_length_tsv_path,
                barcode_length_rows,
                ["barcode_prefix_length", "read_pairs"],
            )
            write_tsv(
                output_paths.cluster_threshold_tsv_path,
                threshold_rows,
                ["reads_threshold_gt", "cluster_count", "reads_in_clusters"],
            )
            write_tsv(
                output_paths.cluster_manifest_tsv_path,
                cluster_manifest_rows,
                [
                    "sample_id",
                    "cluster_batch_label",
                    "cluster_size",
                    "cluster_key",
                    "barcode_file_prefix",
                    "read1_fastq_gz",
                    "read2_fastq_gz",
                ],
            )

            summary_payload = build_pipeline_summary(
                sample_id=sample_id,
                barcode_file=barcode_file,
                barcodes=barcodes,
                fastp_counts=fastp_counts,
                extraction_results=extraction_results,
                barcode_length_counter=barcode_length_counter,
                threshold_rows=threshold_rows,
                total_cluster_count=total_cluster_count,
                filtered_cluster_count=filtered_cluster_count,
                cluster_output_dir=cluster_output_dir,
                args=args,
                resource_plan=resource_plan,
                fastqc_status=fastqc_status,
            )
            sequencing_stats_row = build_sequencing_stats_row(summary_payload)
            write_tsv(output_paths.summary_tsv_path, [summary_payload], list(summary_payload.keys()))
            write_json(output_paths.summary_json_path, summary_payload)
            write_xlsx(
                output_paths.sequencing_stats_xlsx_path,
                [sequencing_stats_row],
                list(sequencing_stats_row.keys()),
                sheet_name="sequencing_stats",
            )
            write_json(
                output_paths.config_json_path,
                {
                    "task_name": TASK_NAME,
                    "version_number": version_number,
                    "sample_id": sample_id,
                    "parameters": {
                        "read1": str(read1_path),
                        "read2": str(read2_path),
                        "project_root": str(project_root) if project_root is not None else None,
                        "barcode_file": str(barcode_file),
                        "existing_clean_dir": str(existing_clean_dir) if existing_clean_dir is not None else None,
                        "barcode_window_start": args.barcode_window_start,
                        "barcode_window_end": args.barcode_window_end,
                        "split_count": args.split_count,
                        "min_read_length": args.min_read_length,
                        "cluster_distance": args.cluster_distance,
                        "min_cluster_reads": args.min_cluster_reads,
                        "disable_fastqc": args.disable_fastqc,
                    },
                    "resource_plan": {
                        "cpu_count": resource_plan.cpu_count,
                        "cpu_usage_percent": resource_plan.cpu_usage_percent,
                        "total_memory_gb": resource_plan.total_memory_gb,
                        "available_memory_gb": resource_plan.available_memory_gb,
                        "available_memory_ratio": resource_plan.available_memory_ratio,
                        "fastp_threads": resource_plan.fastp_threads,
                        "extraction_workers": resource_plan.extraction_workers,
                        "starcode_threads": resource_plan.starcode_threads,
                        "cluster_batch_size": resource_plan.cluster_batch_size,
                        "fastp_split_strategy": "fastp_native_s_split",
                        "cluster_split_strategy": "whole_filtered_clr_batched_barcode_r1_r2_sync_scan",
                    },
                },
            )

            print_and_log(
                f"[COMPLETE] raw={fastp_counts['raw_read_pairs']}, clean={fastp_counts['clean_read_pairs']}, "
                f"barcode_matched={summary_payload['barcode_matched_pairs']}, total_clusters={total_cluster_count}, "
                f"filtered_clusters={filtered_cluster_count}, final_fastq_pairs={summary_payload['generated_barcode_fastq_pairs']}",
                log_handle,
            )
            print_and_log(
                f"[COMPLETE] Sequencing stats: {output_paths.sequencing_stats_xlsx_path}",
                log_handle,
            )
            print_and_log(f"[COMPLETE] Results directory: {output_paths.output_dir}", log_handle)

        return 0
    except Exception as exc:
        print(f"ERROR: Pipeline failed: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
