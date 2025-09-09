# pylint: disable=import-error, wrong-import-position, too-many-locals, too-many-statements
"""
Resolve unmatched AMP taxonomies via manual mapping + NCBI lookup
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time
from typing import Optional, Tuple

# ============================== Third-Party Library Imports ==============================
import numpy as np
import pandas as pd
import yaml
from ete3 import NCBITaxa

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
from preprocess.dbAMP.resolve_taxid_by_blast import get_taxonomy_lineage_from_name
from utils.io_utils import directory_exists, file_exists, load_dataframe_by_columns
from utils.log_utils import get_pipeline_completion_message, get_task_completion_message


# ============================== Custom Functions ==============================
def resolve_manual_taxonomy(
    input_yaml: str,
    output_csv: str,
    amp_input_csv: str,
    amp_output_csv: str,
    logger: CustomLogger,
) -> None:
    """
    Resolve unmatched AMP taxonomy using manually mapped names and NCBI taxonomy.
    Then apply this mapping back to a given AMP table and export updated results.

    Parameters
    ----------
    input_yaml : str
        Path to the YAML file with manual mappings.
    output_csv : str
        Path to save resolved taxonomy results as CSV (mapping table).
    amp_input_csv : str
        Path to AMP data file (e.g., notna_nomatch.csv).
    amp_output_csv : str
        Path to save AMP table with filled-in taxonomy information.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Resolve unmatched taxonomy via manual mapping + NCBI")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Load manual mappings from YAML into DataFrame
        with open(input_yaml, "r", encoding="utf-8") as f:
            manual_mapping = yaml.safe_load(f)
        df_map = pd.DataFrame(
            manual_mapping.items(), columns=["original_name", "mapped_name"]
        )

        # Prepare query stats (before NCBI)
        skip_keywords = {"unknown", "unclassified", "na", "not available", "n/a", ""}
        total_rows = df_map.shape[0]
        to_query_mask = (
            ~df_map["mapped_name"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(skip_keywords)
        )
        num_to_query = int(to_query_mask.sum())
        num_skipped = int(total_rows - num_to_query)

        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Manual mapping query plan ]\n"
                f"▸ total rows        : {total_rows}\n"
                f"  - to query (NCBI)   : {num_to_query}\n"
                f"  - skipped (placeholder): {num_skipped}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Initialize NCBI
        ncbi = NCBITaxa()
        taxonomy_cache: dict[str, Tuple[Optional[dict], Optional[str], int]] = {}

        def cached_tax_lookup(name: str) -> pd.Series:
            """
            Wrapper to cache the results of NCBI taxonomy matching for a given taxon name,
            while skipping known invalid values like 'Unknown' or 'unclassified'.

            Parameters
            ----------
            name : str
                A taxon name (mapped manually) to query from NCBI.

            Returns
            -------
            pd.Series
                A Series of (Tax_Checked, Matched_TaxIDs, TaxID_Count)
            """
            # Define lowercase skip keywords
            skip_keywords = {
                "unknown",
                "unclassified",
                "na",
                "not available",
                "n/a",
                "",
            }
            name_clean = str(name).strip()
            name_lower = name_clean.lower()

            # Skip known invalid placeholders to save time and avoid lookup errors
            if name_lower in skip_keywords:
                return pd.Series([None, None, 0])

            # Return cached result if available
            if name_clean in taxonomy_cache:
                lineage_dict, taxids, count = taxonomy_cache[name_clean]
            else:
                lineage_dict, taxids, count = get_taxonomy_lineage_from_name(
                    name=name_clean, ncbi=ncbi
                )
                taxonomy_cache[name_clean] = (lineage_dict, taxids, count)

            # Convert lineage dict to comma-separated taxonomy string
            if lineage_dict is not None:
                tax_str = ",".join(lineage_dict.values())
            else:
                tax_str = None

            return pd.Series([tax_str, taxids, count])

        # Perform lookup for each mapped name
        df_map[["Tax_Checked", "Matched_TaxIDs", "TaxID_Count"]] = df_map[
            "mapped_name"
        ].apply(cached_tax_lookup)
        df_map["TaxID_Count"] = pd.to_numeric(
            df_map["TaxID_Count"], errors="coerce"
        ).astype("Int64")

        # Post-lookup result stats
        cnt0 = int((df_map["TaxID_Count"] == 0).sum())
        cnt1 = int((df_map["TaxID_Count"] == 1).sum())
        cntm = int((df_map["TaxID_Count"] > 1).sum())
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ NCBI lookup results ]\n"
                f"▸ Total records processed: {df_map.shape[0]}\n"
                f"  - TaxID_Count == 0 (no match): {cnt0}\n"
                f"  - TaxID_Count == 1 (single match): {cnt1}\n"
                f"  - TaxID_Count > 1 (multiple matches): {cntm}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Save result to CSV
        df_map.to_csv(output_csv, index=False)

        # Load AMP CSV via io_utils loader (keep all columns; validates existence)
        df_amp = load_dataframe_by_columns(amp_input_csv, required_columns=None)

        # Rename df_map columns to avoid suffix ambiguity
        df_map_renamed = df_map.rename(
            columns={
                "original_name": "Tax_original",
                "mapped_name": "Tax_mapped",
                "Tax_Checked": "Tax_Checked_mapped",
                "Matched_TaxIDs": "Matched_TaxIDs_mapped",
                "TaxID_Count": "TaxID_Count_mapped",
            }
        )

        # Merge AMP table with mapping table
        df_amp = df_amp.merge(
            df_map_renamed,
            how="left",
            left_on="Tax",
            right_on="Tax_original",
        )
        df_amp["Tax_Checked"] = df_amp["Tax_Checked"].astype("object")
        df_amp["Matched_TaxIDs"] = df_amp["Matched_TaxIDs"].astype("object")

        # Apply replacements only where mapping exists
        has_mapping = df_amp["Tax_mapped"].notna()
        df_amp.loc[has_mapping, "Tax"] = df_amp.loc[has_mapping, "Tax_mapped"]
        df_amp.loc[has_mapping, "Tax_Checked"] = df_amp.loc[
            has_mapping, "Tax_Checked_mapped"
        ]
        df_amp.loc[has_mapping, "Matched_TaxIDs"] = df_amp.loc[
            has_mapping, "Matched_TaxIDs_mapped"
        ]
        df_amp.loc[has_mapping, "TaxID_Count"] = df_amp.loc[
            has_mapping, "TaxID_Count_mapped"
        ]

        # Drop temporary columns
        df_amp.drop(
            columns=[
                "Tax_original",
                "Tax_mapped",
                "Tax_Checked_mapped",
                "Matched_TaxIDs_mapped",
                "TaxID_Count_mapped",
            ],
            inplace=True,
            errors="ignore",
        )

        # Replace 'Unknown' values in Tax with NaN
        df_amp["Tax"] = df_amp["Tax"].replace("Unknown", np.nan)

        # Save result
        df_amp.to_csv(amp_output_csv, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"CSV saved:\n'{output_csv}'\n'{amp_output_csv}'",
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
        logger.exception("Unexpected error in 'resolve_manual_taxonomy()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_prepare_nomatch_mapping(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Resolve unmatched AMP taxonomies using manual mappings (YAML)
    combined with NCBI taxonomy lookup. Saves both the resolved
    mapping table and the updated AMP table.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance for tracking progress.
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
        input_yaml = kwargs.get(
            "nomatch_input_yaml",
            os.path.join(
                base_path, "data/manual/taxid_mapping/notna_nomatch_mapping.yml"
            ),
        )
        output_csv = kwargs.get(
            "nomatch_output_csv",
            os.path.join(
                base_path, "data/processed/dbAMP/resolved_manual_taxonomy.csv"
            ),
        )
        amp_input_csv = kwargs.get(
            "nomatch_amp_input_csv",
            os.path.join(
                base_path, "data/interim/resolve_taxid/02_tax_check/notna_nomatch.csv"
            ),
        )
        amp_output_csv = kwargs.get(
            "nomatch_amp_output_csv",
            os.path.join(
                base_path,
                "data/interim/resolve_taxid/02_tax_check/notna_nomatch_filled.csv",
            ),
        )

        # Check input files exist
        for path in [input_yaml, amp_input_csv]:
            if not file_exists(file_path=path):
                raise FileNotFoundError(f"File not found: '{path}'")

        # Ensure output directories exist
        for path in [output_csv, amp_output_csv]:
            output_dir = os.path.dirname(path)
            if not directory_exists(dir_path=output_dir):
                os.makedirs(name=output_dir)

        # -------------------- Pipeline Execution --------------------
        resolve_manual_taxonomy(
            input_yaml=input_yaml,
            output_csv=output_csv,
            amp_input_csv=amp_input_csv,
            amp_output_csv=amp_output_csv,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_prepare_nomatch_mapping()'")
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
