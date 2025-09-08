# pylint: disable=line-too-long, import-error, wrong-import-position, too-many-locals
"""
Merge checked AMP taxonomy tables (single, multi, filled) into one file,
and optionally append sequence-alignment derived taxonomies (from na_with_blast.csv).
"""
# ============================== Standard Library Imports ==============================
import argparse
import os
import sys

# ============================== Third-Party Library Imports ==============================
import pandas as pd

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
SRC_PATH = os.path.join(BASE_PATH, "src")

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)
if SRC_PATH not in sys.path:
    sys.path.append(SRC_PATH)

# ============================== Project-Specific Imports ==============================
from utils.io_utils import directory_exists, file_exists, load_dataframe_by_columns


# ============================== Main Merge Function ==============================
def merge_checked_tables(
    multi_path: str,
    single_path: str,
    filled_path: str,
    na_with_blast_path: str,
    output_path: str,
) -> None:
    """
    Merge AMP taxonomy check tables from different sources into one,
    keeping only required columns in a fixed order, renaming Tax_Checked to Tax,
    and filtering out rows where Tax is missing. Also appends sequence-alignment
    derived Tax rows from na_with_blast.csv.

    Parameters
    ----------
    multi_path : str
        Path to notna_multi.csv.
    single_path : str
        Path to notna_single.csv.
    filled_path : str
        Path to notna_nomatch_filled.csv.
    na_with_blast_path : str
        Path to na_with_blast.csv (Tax resolved by sequence alignment).
    output_path : str
        Path to save the merged CSV file (without BLAST entries).

    Returns
    -------
    None
    """
    # Validate all input files
    for path in [multi_path, single_path, filled_path, na_with_blast_path]:
        if not file_exists(path):
            raise FileNotFoundError(f"File not found: '{path}'")

    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if not directory_exists(output_dir):
        os.makedirs(output_dir)

    # Required columns per source
    req_checked_cols = ["dbAMP_ID", "Tax_Checked", "Targets", "Seq"]
    req_blast_cols = ["dbAMP_ID", "Tax", "Targets", "Seq"]

    # Load each table using io_utils loader (enforces required columns & handles csv/tsv/xlsx)
    df_multi = load_dataframe_by_columns(multi_path, required_columns=req_checked_cols)
    df_single = load_dataframe_by_columns(
        single_path, required_columns=req_checked_cols
    )
    df_filled = load_dataframe_by_columns(
        filled_path, required_columns=req_checked_cols
    )

    # Merge and clean (checked-only)
    df_merged = pd.concat([df_multi, df_single, df_filled], ignore_index=True)

    # Keep only required columns (already ensured by loader) and rename Tax_Checked -> Tax
    df_merged = df_merged[req_checked_cols].rename(columns={"Tax_Checked": "Tax"})

    # Drop rows with missing Tax and enforce final column order
    final_cols = ["dbAMP_ID", "Tax", "Targets", "Seq"]
    df_merged = df_merged.dropna(subset=["Tax"])[final_cols]

    # Save merged checked result
    df_merged.to_csv(output_path, index=False)

    # Handle na_with_blast.csv (BLAST-derived)
    df_blast = load_dataframe_by_columns(
        na_with_blast_path, required_columns=req_blast_cols
    )
    df_blast = df_blast.dropna(subset=["Tax"])[final_cols]

    # Final merged output (checked + blast)
    df_final = pd.concat([df_merged, df_blast], ignore_index=True)
    df_final = df_final.dropna(subset=["Tax"])[final_cols]  # extra safety

    # Auto-generate new output path with "_with_blast" suffix
    base, ext = os.path.splitext(output_path)
    final_output_path = f"{base}_with_blast{ext}"

    # Save final result with BLAST-based rows included
    df_final.to_csv(final_output_path, index=False)


# ============================== CLI Entry Point ==============================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge AMP taxonomy tables (multi, single, filled) and optionally add BLAST-derived taxonomies"
    )
    parser.add_argument(
        "--multi_path", type=str, required=True, help="Path to notna_multi.csv"
    )
    parser.add_argument(
        "--single_path", type=str, required=True, help="Path to notna_single.csv"
    )
    parser.add_argument(
        "--filled_path",
        type=str,
        required=True,
        help="Path to notna_nomatch_filled.csv",
    )
    parser.add_argument(
        "--na_with_blast_path",
        type=str,
        required=True,
        help="Path to na_with_blast.csv",
    )
    parser.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Path to save the merged CSV file (without BLAST rows)",
    )
    args = parser.parse_args()

    merge_checked_tables(
        multi_path=args.multi_path,
        single_path=args.single_path,
        filled_path=args.filled_path,
        na_with_blast_path=args.na_with_blast_path,
        output_path=args.output_path,
    )
