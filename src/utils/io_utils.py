# pylint: disable=
"""
I/O Utility Module
"""
# ============================== Standard Library Imports ==============================
import os
from typing import List, Optional

# ============================== Third-Party Library Imports ==============================
import pandas as pd


# ============================== Custom Functions ==============================
def file_exists(file_path: str) -> bool:
    """
    Check whether the specified file exists.

    Parameters
    ----------
    file_path : str
        Path to the file to check.

    Returns
    -------
    bool
        True if the file exists, False otherwise.
    """
    return os.path.isfile(file_path)


def directory_exists(dir_path: str) -> bool:
    """
    Check whether the specified directory exists.

    Parameters
    ----------
    dir_path : str
        Path to the directory to check.

    Returns
    -------
    bool
        True if the directory exists, False otherwise.
    """
    return os.path.isdir(dir_path)


def get_missing_columns(df: pd.DataFrame, required_columns: List[str]) -> List[str]:
    """
    Return a list of missing columns from the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to validate.
    required_columns : list of str
        A list of required column names.

    Returns
    -------
    List[str]
        List of missing column names. Empty if all required columns are present.
    """
    return [col for col in required_columns if col not in df.columns]


def load_dataframe_by_columns(
    file_path: str,
    required_columns: Optional[List[str]] = None,
    has_header: bool = True,
) -> pd.DataFrame:
    """
    Load a DataFrame from a file (CSV, TSV, or Excel) and return only the specified columns.
    Optionally handles files without headers.

    Parameters
    ----------
    file_path : str
        Path to the input file (supports .csv, .tsv, .xlsx, .xls).
    required_columns : list of str, optional
        List of required columns to extract. If None, return all columns.
    has_header : bool, optional
        Whether the file contains a header row (default: True).

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the specified or all columns.
    """
    # Check file exists
    if not file_exists(file_path=file_path):
        raise FileNotFoundError(f"File not found: '{file_path}'")

    # Detect file type
    ext = os.path.splitext(file_path)[-1].lower()

    # Header mode
    header_opt = 0 if has_header else None

    # Load DataFrame
    if ext == ".csv":
        df = pd.read_csv(file_path, header=header_opt)
    elif ext == ".tsv":
        df = pd.read_csv(file_path, sep="\t", header=header_opt)
    elif ext in [".xlsx", ".xls"]:
        df = pd.read_excel(file_path, header=header_opt)
    else:
        raise ValueError(f"Unsupported file extension: '{ext}'")

    # If no header, assign temporary column names
    if not has_header:
        df.columns = [f"column_{i}" for i in range(df.shape[1])]

    # Use all columns if none specified
    if required_columns is None:
        return df.copy()

    # Check required columns
    missing = get_missing_columns(df=df, required_columns=required_columns)
    if missing:
        raise ValueError(f"Missing required columns: '{str(missing)}'")

    return df[required_columns].copy()


def write_readme(
    dir_path: str, title: str, description: str, extra_info: Optional[dict] = None
) -> None:
    """
    Create a README.md file in the specified directory with step information.

    Parameters
    ----------
    dir_path : str
        Path to the step directory.
    title : str
        Title of the step.
    description : str
        Short description of what this step does.
    extra_info : dict, optional
        Extra key-value info to append (e.g., source URL, output location, date).
    """
    # Check folder exists
    if not directory_exists(dir_path=dir_path):
        os.makedirs(dir_path)

    readme_path = os.path.join(dir_path, "README.md")
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write(f"# {title}\n\n")
        f.write(f"{description}\n\n")
        if extra_info:
            for key, value in extra_info.items():
                f.write(f"**{key}:** {value}\n")
