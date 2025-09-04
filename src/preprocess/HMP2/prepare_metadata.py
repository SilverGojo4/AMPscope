# pylint: disable=import-error, wrong-import-position
"""
Prepare HMP2 metadata for nf-core/ampliseq
"""
# ============================== Standard Library Imports ==============================
import argparse
import os
import sys

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
SRC_PATH = os.path.join(BASE_PATH, "src")

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)
if SRC_PATH not in sys.path:
    sys.path.append(SRC_PATH)

# ============================== Project-Specific Imports ==============================
# Local AMPscope utility modules
from utils.io_utils import directory_exists, file_exists, load_dataframe_by_columns


# ============================== Main Processing Function ==============================
def prepare_metadata(input_csv: str, output_tsv: str) -> None:
    """
    Filter HMP2 metadata for biopsy_16S samples and export standardized sample info.

    Parameters
    ----------
    input_csv : str
        Path to the input raw HMP2 metadata CSV.
    output_tsv : str
        Path to save the filtered metadata TSV.
    """
    # Check input file exists
    if not file_exists(file_path=input_csv):
        raise FileNotFoundError(f"File not found: '{input_csv}'")

    # Ensure output directory exists
    output_dir = os.path.dirname(output_tsv)
    if not directory_exists(dir_path=output_dir):
        os.makedirs(name=output_dir)

    # Load relevant columns
    required_columns = [
        "Project",
        "External ID",
        "Participant ID",
        "data_type",
        "biopsy_location",
    ]
    df = load_dataframe_by_columns(
        file_path=input_csv, required_columns=required_columns
    )

    # Filter for biopsy_16S
    df_filtered = df[df["data_type"] == "biopsy_16S"]

    # Keep only required columns
    df_final = df_filtered[["External ID", "biopsy_location"]].copy()

    # Add prefix 's' to External ID values
    df_final["External ID"] = df_final["External ID"].apply(lambda x: f"s{x}")

    # Rename columns
    df_final.rename(
        columns={
            "External ID": "sample_name",
            "biopsy_location": "sample_index",
        },
        inplace=True,
    )

    # Save to TSV
    df_final.to_csv(output_tsv, sep="\t", index=False)
    print(f"Filtered metadata saved to: {output_tsv}")


# ============================== CLI Entry Point ==============================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare HMP2 metadata for nf-core/ampliseq"
    )
    parser.add_argument(
        "--input_csv",
        type=str,
        default=os.path.join(BASE_PATH, "data/raw/HMP2/hmp2_metadata_2018-08-20.csv"),
        help="Path to the raw HMP2 metadata CSV file",
    )
    parser.add_argument(
        "--output_tsv",
        type=str,
        default=os.path.join(BASE_PATH, "data/interim/nfcore_metadata.tsv"),
        help="Path to save the filtered metadata TSV",
    )
    args = parser.parse_args()

    prepare_metadata(input_csv=args.input_csv, output_tsv=args.output_tsv)
