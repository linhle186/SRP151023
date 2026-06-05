#!/bin/bash
# run_pipeline.sh
# 
# This script automates the complete single-cell RNA-seq workflow:
# 1. Downloads, processes, and builds STAR indexes for both human and mouse.
# 2. Downloads raw SRA data, performs alignment, velocyto run, and QC for SRP151023.

# Exit immediately if any command exits with a non-zero status
set -e

# Path to the preconfigured cluster profile
PROFILE_PATH="/home/share/ag_klughammer/cluster_setup/snakemake_cluster_configs/slurm_nice" # <-- UPDATE THIS PATH TO THE ACTUAL SNAKEMAKE PROFILE

echo "=========================================================="
echo "         STARTING AUTOMATED SINGLE-CELL PIPELINE          "
echo "=========================================================="
echo "Execution Date: $(date)"
echo "Working Directory: $(pwd)"
echo "SLURM Profile: $PROFILE_PATH"
echo "=========================================================="
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
          --resources gsm_download_slots=4 star_slots=2 velocyto_slots=2 \
          --latency-wait 120

echo ""
echo "=========================================================="
echo "         PIPELINE RUN SUCCESSFULLY COMPLETED!             "
echo "=========================================================="
echo "Completion Date: $(date)"