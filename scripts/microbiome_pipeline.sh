#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Check input arguments
if [ $# -lt 1 ]; then
  echo "Usage: $0 <PROJECT_DIR>"
  exit 1
fi

# Set paths
PROJECT_DIR=$1
LOG_FILE_MICRO="$PROJECT_DIR/logs/microbiome_pipeline.log"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate AMPscope

# ===================== Step 1: Prepare metadata =====================
rm -f "$LOG_FILE_MICRO"
rm -f "$PROJECT_DIR/data/interim/nfcore_metadata.tsv"
python "$PROJECT_DIR/src/preprocess/HMP2/prepare_metadata.py" \
  --input_csv "$PROJECT_DIR/data/raw/HMP2/hmp2_metadata_2018-08-20.csv" \
  --output_tsv "$PROJECT_DIR/data/interim/nfcore_metadata.tsv"

# ===================== Step 2: Microbiome composition analysis =====================
rm -rf "$PROJECT_DIR/.nextflow.log"
rm -rf "$PROJECT_DIR/experiments/HMP2"
rm -rf "$PROJECT_DIR/work"
rm -rf "$PROJECT_DIR/.nextflow"
python "$PROJECT_DIR/src/main.py" \
  --stage microbiome_composition \
  --log_path "$LOG_FILE_MICRO" \
  --microbiome_input_csv "$PROJECT_DIR/data/raw/HMP2/biopsy_16S/samplesheet.csv" \
  --microbiome_metadata "$PROJECT_DIR/data/interim/nfcore_metadata.tsv" \
  --microbiome_silva_train "$PROJECT_DIR/data/external/silva/nr99_v138.1_wSpecies_train_set.fa.gz" \
  --microbiome_silva_species "$PROJECT_DIR/data/external/silva/species_assignment_v138.1.fa.gz" \
  --microbiome_output_dir "$PROJECT_DIR/experiments/HMP2" \
  --microbiome_dada2_threads 8 \
  --microbiome_fastp_threads 8 \
  --microbiome_cutadapt_threads 8 \
  --microbiome_max_cpus 24 \
  --microbiome_max_memory 256.GB \
  --microbiome_max_time 48.h

# Deactivate the Conda environment
conda deactivate