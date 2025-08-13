#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs/slurm
snakemake --profile profiles/slurm
