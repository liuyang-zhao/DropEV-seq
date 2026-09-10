# DropEV-seq analysis pipelines

This repository contains the code used for two DropEV-seq barcode-level
analysis workflows for paired-end sequencing data:

1. `drop_ev_barcode_pipeline.py` performs read QC, barcode extraction,
   barcode clustering, barcode-specific FASTQ splitting, and summary report
   generation.
2. `dropev_alignment_purity.py` aligns barcode-split reads to a combined
   reference, calculates KP/EC alignment summaries, and performs a read-count
   threshold sweep.

Only source code and documentation are included. No FASTQ, BAM, reference
genome, sample table, or generated result is distributed in this repository.

## Repository layout

```text
.
├── scripts/
│   ├── drop_ev_barcode_pipeline.py
│   └── dropev_alignment_purity.py
├── .gitignore
├── environment.yml
├── LICENSE
└── README.md
```

The scripts are kept under `scripts/` because their default project-root
resolution assumes that layout. Use `--project-root` when a different output
location is required.

## Software requirements

- Linux or another POSIX-like environment
- Python 3.10 or newer
- Conda environment named `evdna` (see `environment.yml`)
- Python packages: `pandas`, `openpyxl`, `psutil`, and `pysam`
- Command-line tools: `fastp`, `fastqc`, `starcode`, `bwa`, and `samtools`

The two workflows are independent: the first produces barcode-split FASTQ
files, while the second consumes those files and a combined reference FASTA.
Paths to sequencing data and references are supplied at runtime and are not
hard-coded in this repository.

## Installation

Use the project environment rather than the system Python installation:

```bash
conda env create -f environment.yml
conda activate evdna
```

Install the external bioinformatics tools through Conda or the software
management system approved by your laboratory. Record their versions in the
run log before processing data.

## Workflow 1: barcode extraction and clustering

Run from the repository root:

```bash
conda run -n evdna python scripts/drop_ev_barcode_pipeline.py \
  --read1 /path/to/sample_R1.fastq.gz \
  --read2 /path/to/sample_R2.fastq.gz \
  --barcode-file /path/to/barcode_whitelist.txt \
  --sample-id sample001 \
  --project-root /path/to/analysis_project
```

Important defaults include an R2 barcode search slice of zero-based indices
13 through 23 (`[13:24]`), 60 split files, a minimum retained read length of
120, and a strict cluster read cutoff of more than 1,000 reads. Use `--help`
to inspect all available parameters. The pipeline writes versioned output and
log directories below the selected project root.

## Workflow 2: alignment and purity threshold sweep

```bash
conda run -n evdna python scripts/dropev_alignment_purity.py \
  --barcode-fastq-dir /path/to/barcode_split_fastq \
  --reference-fasta /path/to/combined_reference.fasta \
  --sample-id d0 \
  --mapq-threshold 60 \
  --min-total-reads 10 \
  --purity-threshold-percent 85 \
  --reads-thresholds 200
```

The alignment workflow uses strict defaults of MAPQ >= 60, alignment identity >= 0.95, aligned fraction >= 0.80, and a minimum of 10 total mapped reads per barcode. Sorted per-barcode BAM files are enabled by default and can be disabled with `--no-save-barcode-bam` to reduce disk usage. Use `--help` for the complete parameter list. The combined reference must have a `.fai` index and BWA index files (`.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`). Reference contig names are expected to carry `kp_` or `ec_` prefixes so that species assignment and length summaries can be calculated correctly.

## Validation

The source files pass a Python syntax compilation check:

```bash
conda run -n evdna python -m py_compile scripts/*.py
```

Syntax validation does not replace a complete end-to-end test with the
required external tools and representative input data. Validate those tools,
reference indexes, and resource limits in the target compute environment
before processing a production dataset.

## Data and privacy policy

Do not commit raw reads, assembled references, alignment files, run logs,
machine-specific configuration, or generated reports. The `.gitignore` file
contains conservative rules for common sequencing and analysis artifacts, but
review `git status` before every commit. Runtime logs may contain absolute
input paths; inspect and redact them before sharing.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.
