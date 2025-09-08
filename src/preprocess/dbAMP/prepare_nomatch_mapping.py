# pylint: disable=import-error, wrong-import-position
"""
Resolve unmatched AMP taxonomies via manual mapping + NCBI lookup
"""
# ============================== Standard Library Imports ==============================
import argparse
import os
import sys
from typing import Optional, Tuple

# ============================== Third-Party Library Imports ==============================
import numpy as np
import pandas as pd
import yaml
from ete3 import NCBITaxa

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
SRC_PATH = os.path.join(BASE_PATH, "src")

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)
if SRC_PATH not in sys.path:
    sys.path.append(SRC_PATH)


# ============================== Project-Specific Imports ==============================
from src.data.dbAMP.resolve_taxid_by_blast import get_taxonomy_lineage_from_name

# Local AMPscope utility modules
from utils.io_utils import directory_exists, file_exists


# ============================== Main Processing Function ==============================
def resolve_manual_taxonomy(
    input_yaml: str,
    output_csv: str,
    amp_input_csv: str,
    amp_output_csv: str,
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

    Returns
    -------
    None
    """
    # Check input files exist
    for path in [input_yaml, amp_input_csv]:
        if not file_exists(file_path=path):
            raise FileNotFoundError(f"File not found: '{path}'")

    # Ensure output directories exist
    for path in [output_csv, amp_output_csv]:
        output_dir = os.path.dirname(path)
        if not directory_exists(dir_path=output_dir):
            os.makedirs(name=output_dir)

    # Load manual mappings from YAML into DataFrame
    with open(input_yaml, "r", encoding="utf-8") as f:
        manual_mapping = yaml.safe_load(f)
    df_map = pd.DataFrame(
        manual_mapping.items(), columns=["original_name", "mapped_name"]
    )

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
        skip_keywords = {"unknown", "unclassified", "na", "not available", "n/a", ""}
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

    # Save result to TSV
    df_map.to_csv(output_csv, index=False)

    # Load AMP CSV
    df_amp = pd.read_csv(amp_input_csv)

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
    )

    # Replace 'Unknown' values in Tax with 'NA'
    df_amp["Tax"] = df_amp["Tax"].replace("Unknown", np.nan)

    # Save result
    df_amp.to_csv(amp_output_csv, index=False)


# ============================== CLI Entry Point ==============================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Resolve unmatched AMP taxonomy using manual mapping + NCBI"
    )
    parser.add_argument(
        "--input_yaml",
        type=str,
        required=True,
        help="Path to the manually mapped taxonomy YAML file",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        required=True,
        help="Path to save the resolved taxonomy CSV file (mapping table)",
    )
    parser.add_argument(
        "--amp_input_csv",
        type=str,
        required=True,
        help="Path to the AMP table (e.g., notna_nomatch.csv) to apply mapping to",
    )
    parser.add_argument(
        "--amp_output_csv",
        type=str,
        required=True,
        help="Path to save the updated AMP table with mapped taxonomy",
    )
    args = parser.parse_args()

    resolve_manual_taxonomy(
        input_yaml=args.input_yaml,
        output_csv=args.output_csv,
        amp_input_csv=args.amp_input_csv,
        amp_output_csv=args.amp_output_csv,
    )
