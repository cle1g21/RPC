#!/bin/bash
#SBATCH --job-name=proteome_filt
#SBATCH --partition=amd_serial
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/home/cle1g21/RPC/confidence_pipeline/logs/proteome_filt_%j.out
#SBATCH --error=/home/cle1g21/RPC/confidence_pipeline/logs/proteome_filt_%j.err

set -eo pipefail
mkdir -p /home/cle1g21/RPC/confidence_pipeline/logs

echo "=== Proteome filter root InstaNovo CSVs ==="
echo "Job ID:  ${SLURM_JOB_ID:-local}"
echo "Node:    $(hostname)"
echo "Start:   $(date)"

python3 /home/cle1g21/RPC/confidence_pipeline/filter_root_instanovo_vs_proteome.py

echo "End:     $(date)"
echo "=== Done ==="
