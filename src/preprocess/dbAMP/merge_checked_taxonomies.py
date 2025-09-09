# pylint: disable=import-error, wrong-import-position, too-many-arguments, too-many-positional-arguments, too-many-locals
"""
Merge checked AMP taxonomy tables (single, multi, filled) into one file,
and optionally append sequence-alignment derived taxonomies (from na_with_blast.csv).
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time

# ============================== Third-Party Library Imports ==============================
import pandas as pd

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
LOGGING_PATH = os.path.join(BASE_PATH, "src/utils/logging_toolkit/src/python")

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)
if LOGGING_PATH not in sys.path:
    sys.path.append(LOGGING_PATH)

# ============================== Project-Specific Imports ==============================
# Submodule imports (external tools integrated into project)
from setup_logging import CustomLogger

# Local AMPscope utility modules
from utils.io_utils import directory_exists, file_exists, load_dataframe_by_columns
from utils.log_utils import get_pipeline_completion_message, get_task_completion_message


# ============================== Main Merge Function ==============================
def merge_checked_tables(
    multi_path: str,
    single_path: str,
    filled_path: str,
    na_with_blast_path: str,
    output_path: str,
    logger: CustomLogger,
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
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Merge taxonomy tables (checked + optional BLAST)")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Required columns per source
        req_checked_cols = ["dbAMP_ID", "Tax_Checked", "Targets", "Seq"]
        req_blast_cols = ["dbAMP_ID", "Tax", "Targets", "Seq"]

        # Load each table using io_utils loader (enforces required columns & handles csv/tsv/xlsx)
        df_multi = load_dataframe_by_columns(
            multi_path, required_columns=req_checked_cols
        )
        df_single = load_dataframe_by_columns(
            single_path, required_columns=req_checked_cols
        )
        df_filled = load_dataframe_by_columns(
            filled_path, required_columns=req_checked_cols
        )

        # Merge and clean (checked-only)
        df_merged_raw = pd.concat([df_multi, df_single, df_filled], ignore_index=True)
        before_drop = df_merged_raw.shape[0]

        # Keep only required columns (already ensured by loader) and rename Tax_Checked -> Tax
        df_merged = df_merged_raw[req_checked_cols].rename(
            columns={"Tax_Checked": "Tax"}
        )

        # Drop rows with missing Tax and enforce final column order
        final_cols = ["dbAMP_ID", "Tax", "Targets", "Seq"]
        df_merged = df_merged.dropna(subset=["Tax"])[final_cols]
        after_drop = df_merged.shape[0]

        # Save merged checked result
        df_merged.to_csv(output_path, index=False)

        # Handle na_with_blast.csv (BLAST-derived)
        df_blast_raw = load_dataframe_by_columns(
            na_with_blast_path, required_columns=req_blast_cols
        )
        blast_before = df_blast_raw.shape[0]
        df_blast = df_blast_raw.dropna(subset=["Tax"])[final_cols]
        blast_after = df_blast.shape[0]
        blast_dropped = blast_before - blast_after

        # Final merged output (checked + blast)
        df_final = pd.concat([df_merged, df_blast], ignore_index=True)
        df_final = df_final.dropna(subset=["Tax"])[final_cols]  # extra safety
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Merge summary ]\n"
                f"▸ Checked-only merge (Single + Multi + Filled):\n"
                f"  - Rows before filter : {before_drop}\n"
                f"  - Rows dropped (Tax NA): {before_drop - after_drop}\n"
                f"  - Rows after filter  : {after_drop}\n"
                f"▸ BLAST-derived table:\n"
                f"  - Rows loaded        : {blast_before}\n"
                f"  - Rows dropped (Tax NA): {blast_dropped}\n"
                f"  - Rows retained      : {blast_after}\n"
                f"▸ Final combined output:\n"
                f"  - Total rows (checked + BLAST): {df_final.shape[0]}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Auto-generate new output path with "_with_blast" suffix
        base, ext = os.path.splitext(output_path)
        final_output_path = f"{base}_with_blast{ext}"

        # Save final result with BLAST-based rows included
        df_final.to_csv(final_output_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"CSV saved:\n'{output_path}'\n'{final_output_path}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Summary
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'merge_checked_tables()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_merge_checked_taxonomies(
    base_path: str, logger: CustomLogger, **kwargs
) -> None:
    """
    Merge AMP taxonomy tables (single, multi, filled) with optional
    BLAST-derived taxonomies. Saves both the merged checked table
    and the final combined table including BLAST rows.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance to track progress.
    **kwargs : dict
        Optional overrides for input/output file paths.

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        # -------------------- Retrieve input parameters (CLI or default) --------------------
        multi_path = kwargs.get(
            "merge_multi_path",
            os.path.join(
                base_path, "data/interim/resolve_taxid/02_tax_check/notna_multi.csv"
            ),
        )
        single_path = kwargs.get(
            "merge_single_path",
            os.path.join(
                base_path, "data/interim/resolve_taxid/02_tax_check/notna_single.csv"
            ),
        )
        filled_path = kwargs.get(
            "merge_filled_path",
            os.path.join(
                base_path,
                "data/interim/resolve_taxid/02_tax_check/notna_nomatch_filled.csv",
            ),
        )
        na_with_blast_path = kwargs.get(
            "merge_na_with_blast_path",
            os.path.join(base_path, "data/interim/blast/06_merge/na_with_blast.csv"),
        )
        output_path = kwargs.get(
            "merge_output_path",
            os.path.join(base_path, "data/processed/dbAMP/checked_taxonomies.csv"),
        )

        # Validate all input files
        for path in [multi_path, single_path, filled_path, na_with_blast_path]:
            if not file_exists(path):
                raise FileNotFoundError(f"File not found: '{path}'")

        # Ensure output directory exists
        output_dir = os.path.dirname(output_path)
        if not directory_exists(output_dir):
            os.makedirs(output_dir)

        # -------------------- Pipeline Execution --------------------
        merge_checked_tables(
            multi_path=multi_path,
            single_path=single_path,
            filled_path=filled_path,
            na_with_blast_path=na_with_blast_path,
            output_path=output_path,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_merge_checked_taxonomies()'")
        raise

    # Final summary block
    logger.info("[ 'Pipeline Execution Summary' ]")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="=")
    logger.log_with_borders(
        level=logging.INFO,
        message="\n".join(get_pipeline_completion_message(start_time)),
        border="║",
        length=120,
    )
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="=")
