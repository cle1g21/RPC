#!/bin/bash
#SBATCH --job-name=conf_tiers
#SBATCH --partition=amd_serial
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/home/cle1g21/RPC/confidence_pipeline/logs/conf_tiers_%j.out
#SBATCH --error=/home/cle1g21/RPC/confidence_pipeline/logs/conf_tiers_%j.err

set -eo pipefail
mkdir -p /home/cle1g21/RPC/confidence_pipeline/logs

echo "=== Update confidence tiers + assemble masters ==="
echo "Job ID:  ${SLURM_JOB_ID:-local}"
echo "Node:    $(hostname)"
echo "Start:   $(date)"

python3 /home/cle1g21/RPC/confidence_pipeline/update_confidence_tiers.py

echo "End:     $(date)"
echo "=== Done ==="
