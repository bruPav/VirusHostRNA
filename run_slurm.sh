#!/bin/bash
# =============================================================================
# Submit the VirusHostRNA pipeline to a SLURM cluster.
#
# Prerequisites:
#   1. Edit config.yaml to match your data and compute environment.
#   2. Adjust profiles/slurm/config.yaml for your cluster (partition, defaults).
#
# Usage:
#   chmod +x run_slurm.sh
#   ./run_slurm.sh
# =============================================================================
set -euo pipefail

echo "Starting VirusHostRNA pipeline on SLURM cluster..."
echo "Profile:  profiles/slurm"
echo "Max jobs: 50"

snakemake \
    --profile profiles/slurm \
    --use-conda \
    --jobs 50 \
    --latency-wait 60 \
    --restart-times 2

echo "Pipeline submission complete. Check status with: squeue -u $USER"
echo "Or monitor with: snakemake --profile profiles/slurm --report results/report.html"
