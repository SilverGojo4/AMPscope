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
LOG_FILE="$PROJECT_DIR/logs/taxonomy_pipeline.log"

# Get the base directory of the Conda installation
CONDA_BASE=$(conda info --base)

# Source the Conda initialization script to enable 'conda' commands
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Activate the Conda environment
conda activate AMPscope

# ===================== Step 1: Run BLAST-based TaxID Resolution =====================
rm -f "$LOG_FILE"
rm -rf "$PROJECT_DIR/data/interim/resolve_taxid"
python "$PROJECT_DIR/src/main.py" \
  --stage resolve_taxid_by_blast \
  --log_path "$LOG_FILE" \
  --resolve_input_dbamp "$PROJECT_DIR/data/raw/dbAMP/dbAMP3_pepinfo.xlsx" \
  --resolve_output_dir "$PROJECT_DIR/data/interim/resolve_taxid" \
  --resource_dir "$PROJECT_DIR/data/external/uniprot" \
  --resolve_force_download \
  --resolve_force_build \
  --resolve_blast_evalue 1e-8 \
  --resolve_blast_max_hits 50 \
  --resolve_blast_max_hsps 1 \
  --resolve_blast_pident_min 50 \
  --resolve_blast_qcovs_min 70 \
  --resolve_blast_threads 24 \
  --resolve_blast_pe_max 2

# ===================== Step 2: Run Manual Mapping-based Taxonomy Supplement =====================
rm -rf "$PROJECT_DIR/data/processed/dbAMP/resolved_manual_taxonomy.csv"
rm -rf "$PROJECT_DIR/data/interim/resolve_taxid/02_tax_check/notna_nomatch_filled.csv"
python "$PROJECT_DIR/src/main.py" \
  --stage taxonomy_nomatch_mapping \
  --log_path "$LOG_FILE" \
  --nomatch_input_yaml "$PROJECT_DIR/data/manual/taxid_mapping/notna_nomatch_mapping.yml" \
  --nomatch_output_csv "$PROJECT_DIR/data/processed/dbAMP/resolved_manual_taxonomy.csv" \
  --nomatch_amp_input_csv "$PROJECT_DIR/data/interim/resolve_taxid/02_tax_check/notna_nomatch.csv" \
  --nomatch_amp_output_csv "$PROJECT_DIR/data/interim/resolve_taxid/02_tax_check/notna_nomatch_filled.csv"

# ===================== Step 3: Merge Checked Taxonomies =====================
rm -rf "$PROJECT_DIR/data/processed/dbAMP/checked_taxonomies.csv"
rm -rf "$PROJECT_DIR/data/processed/dbAMP/checked_taxonomies_with_blast.csv"
python "$PROJECT_DIR/src/main.py" \
  --stage merge_checked_taxonomies \
  --log_path "$LOG_FILE" \
  --merge_multi_path "$PROJECT_DIR/data/interim/resolve_taxid/02_tax_check/notna_multi.csv" \
  --merge_single_path "$PROJECT_DIR/data/interim/resolve_taxid/02_tax_check/notna_single.csv" \
  --merge_filled_path "$PROJECT_DIR/data/interim/resolve_taxid/02_tax_check/notna_nomatch_filled.csv" \
  --merge_na_with_blast_path "$PROJECT_DIR/data/interim/resolve_taxid/06_merge/na_with_blast.csv" \
  --merge_output_path "$PROJECT_DIR/data/processed/dbAMP/checked_taxonomies.csv"

# Deactivate the Conda environment
conda deactivate
