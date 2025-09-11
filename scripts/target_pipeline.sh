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
LOG_FILE_BACDIVE="$PROJECT_DIR/logs/resolve_target_by_bacdive.log"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate AMPscope

# ===================== Step 1: Run BLAST-based TaxID Resolution =====================
rm -f "$LOG_FILE_BACDIVE"
rm -rf "$PROJECT_DIR/data/interim/resolve_targets"
python "$PROJECT_DIR/src/main.py" \
  --stage clean_targets \
  --log_path "$LOG_FILE_BACDIVE" \
  --targets_input_path "$PROJECT_DIR/data/raw/dbAMP/dbAMP3_pepinfo.xlsx" \
  --targets_output_csv "$PROJECT_DIR/data/interim/resolve_targets/targets_clean.csv"

# ===================== Step 2: Run Target Mapping via BacDive =====================
rm -rf "$PROJECT_DIR/data/processed/dbAMP/resolved_manual_targets.csv"
rm -rf "$PROJECT_DIR/data/processed/dbAMP/checked_targets.csv"
python "$PROJECT_DIR/src/main.py" \
  --stage targets_mapping \
  --log_path "$LOG_FILE_BACDIVE" \
  --targetsmap_input_yaml "$PROJECT_DIR/data/manual/targets_mapping/targets_mapping.yml" \
  --targetsmap_bacdive_config "$PROJECT_DIR/configs/bacdive.json" \
  --targetsmap_output_csv "$PROJECT_DIR/data/processed/dbAMP/resolved_manual_targets.csv" \
  --targetsmap_amp_input_csv "$PROJECT_DIR/data/interim/resolve_targets/targets_clean.csv" \
  --targetsmap_amp_output_csv "$PROJECT_DIR/data/processed/dbAMP/checked_targets.csv"


# Deactivate the Conda environment
conda deactivate
