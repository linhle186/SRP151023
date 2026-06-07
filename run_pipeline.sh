#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --mail-user=gyu@genzentrum.lmu.de
#SBATCH --mail-type=fail
#SBATCH -o %j.log
#SBATCH -e %j.err
#SBATCH -J RNAvelo.pp.pipeline
#SBATCH -t 7-00:00:00
# 
# This script automates the complete single-cell RNA-seq workflow:
# 1. Downloads, processes, and builds STAR indexes for both human and mouse.
# 2. Downloads raw SRA data, performs alignment, velocyto run, and QC for SRP151023.

# Exit immediately if any command exits with a non-zero status
set -e

# Path to the preconfigured cluster profile
PROFILE_PATH="/home/share/ag_klughammer/cluster_setup/snakemake_cluster_configs/slurm_nice" # SFB
PROFILE_PATH="/raid/ag_klughammer/shared/cluster_setup/snakemake_cluster_configs/slurm_nice" # DGX

echo "=========================================================="
echo "         STARTING AUTOMATED SINGLE-CELL PIPELINE          "
echo "=========================================================="
echo "Execution Date: $(date)"
echo "Working Directory: $(pwd)"
echo "SLURM Profile: $PROFILE_PATH"
echo "=========================================================="
echo ""

# --- STAGE -1: ENVIRONMENT INITIALIZATION ---
echo "Checking for required Conda environment..."
# Check if scanpy_env is already installed on the server
if ! conda info --envs | grep -q "scanpy_env"; then
    echo "Environment 'scanpy_env' not found. Creating it from environment.yml..."
    conda env create -f environment.yml -n scanpy_env
    echo "Environment successfully created."
else
    echo "Environment 'scanpy_env' already exists. Skipping creation."
fi
echo ""

# --- STAGE 0: UNLOCK WORKING DIRECTORIES ---
echo "1. Unlocking Snakemake working directories..."
snakemake -s Snakefile_prep --unlock 2>/dev/null || true
snakemake -s Snakefile --unlock 2>/dev/null || true
echo "Unlock completed."
echo ""

# --- STAGE 1: GENOME REFERENCE PREPARATION ---
echo "----------------------------------------------------------"
echo "STAGE 1: Preparing Genome References (Human & Mouse)..."
echo "----------------------------------------------------------"
echo "Downloading reference sequences, GTFs, whitelists, and"
echo "generating STAR indexes. This may take several hours."
echo ""

snakemake -s Snakefile_prep \
          --profile "$PROFILE_PATH" \
          --jobs 4 \
          --resources mem_mb=120000 \
          --latency-wait 120

echo ""
echo "Reference preparation successfully completed."
echo ""

# --- STAGE 2: MAIN SINGLE-CELL RNA-SEQ PROCESSING ---
echo "----------------------------------------------------------"
echo "STAGE 2: Executing Main Single-Cell RNA-seq Processing..."
echo "----------------------------------------------------------"
echo "Processing dataset SRP151023 (Throttled download, alignment,"
echo "velocyto execution, and downstream matrix aggregation)..."
echo ""

snakemake -s Snakefile \
          --profile "$PROFILE_PATH" \
          --jobs 6 \
          --resources download_slots=4 star_slots=2 velocyto_slots=2 \
          --latency-wait 120

echo ""
echo "=========================================================="
echo "         PIPELINE RUN SUCCESSFULLY COMPLETED!             "
echo "=========================================================="
echo "Completion Date: $(date)"