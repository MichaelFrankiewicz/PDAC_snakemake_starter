# PDAC RNA-seq Snakemake Starter

This repo is a minimal, production-ready pipeline that runs **QC → alignment (STAR) → BAM index → bigWig** and aggregates reports with **MultiQC**. It is tailored for Slurm clusters and wet+dry PDAC projects.

## Quick start (local)
```bash
# 1) Create base env with snakemake
mamba create -n pdac -c conda-forge -c bioconda snakemake
conda activate pdac

# 2) Edit config/config.yaml and config/samples.csv
#    - Provide paths to hg38 FASTA + GTF or download them, and set star_index dir
#    - Put your FASTQ paths or SRA accessions in samples.csv

# 3) Dry-run
snakemake -n

# 4) Run (local cores)
snakemake --cores 8
```

## Quick start (Slurm)
```bash
# Ensure the logs folder exists
mkdir -p logs/slurm

# Use the included Slurm profile
snakemake --profile profiles/slurm
```

## Inputs
- `config/config.yaml`: global settings, reference files, STAR index location
- `config/samples.csv`: sample sheet (supports either FASTQ paths or SRA accessions)

## Outputs
- `results/fastqc/*_fastqc.html(.zip)`
- `results/bam/{sample}.bam` and `.bai`
- `results/bigwig/{sample}.bw`
- `results/multiqc/multiqc_report.html`

## Notes
- If `reference.star_index` is empty, the pipeline will build it from `reference.fasta` + `reference.gtf`.
- If `samples.csv` provides an `sra` accession (e.g., SRR123...), the pipeline will fetch FASTQs with `fasterq-dump`.
- If you already have FASTQs, leave `sra` empty and set `fastq_1` (and `fastq_2` if paired-end).

## Recommended references (hg38)
- FASTA: Ensembl GRCh38 primary assembly
- GTF: Ensembl gene annotation (versioned, e.g., 109)

> Tip: Commit this repo, then open a PR for each change (parameters, new rules, etc.) to showcase reproducible research skills.
