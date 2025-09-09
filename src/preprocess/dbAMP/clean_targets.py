# pylint: disable=line-too-long, import-error, wrong-import-position, too-many-locals
"""
Clean AMP 'Targets' text into a normalized long table (text cleaning only).
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import re
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

# ============================== Precompiled Regex Patterns ==============================
_BRACKETS_RE = re.compile(r"\s*\(.*?\)")  # remove bracketed notes e.g. (MIC=25.4µM)
_PREFIX_RE = re.compile(
    r"^[A-Za-z\s\-]+:\s*"
)  # remove leading category like 'Bacteria: '
_MULTI_SPACE_RE = re.compile(r"\s+")  # collapse spaces
_SPACE_BEFORE_COMMA_RE = re.compile(r"\s+,")  # normalize " ," -> ","
_DOUBLE_COMMA_RE = re.compile(r",,+")  # collapse ",,"
_SPACE_AROUND_PAREN_RE = re.compile(
    r"\s*\(\s*|\s*\)\s*"
)  # normalize spaces around parentheses


# ============================== Custom Functions ==============================
def clean_target_string(target_str: str) -> str:
    """
    Normalize a single target string.

    Steps
    -----
    - Trim spaces
    - Replace 'μ' with 'µ'
    - Remove bracketed notes like '(MIC=25.4µM)'
    - Remove category prefixes like 'Bacteria: '
    - Collapse excessive spaces
    - Normalize spacing around punctuation and strip trailing punctuation

    Parameters
    ----------
    target_str : str
        Raw target string.

    Returns
    -------
    str
        Cleaned target string (may be empty).
    """
    s = str(target_str).strip()
    if not s:
        return ""

    # Normalize micro sign
    if "μ" in s:
        s = s.replace("μ", "µ")

    # Normalize punctuation variants
    if "；" in s:
        s = s.replace("；", ";")

    # Remove bracketed notes (e.g., MIC/IC50 descriptions)
    s = _BRACKETS_RE.sub("", s)

    # Remove leading category prefixes like 'Bacteria: '
    s = _PREFIX_RE.sub("", s)

    # Collapse excessive spaces
    s = _MULTI_SPACE_RE.sub(" ", s).strip()

    # Normalize spaces around punctuation/parens
    s = _SPACE_BEFORE_COMMA_RE.sub(",", s)  # " ," -> ","
    s = _DOUBLE_COMMA_RE.sub(",", s)  # ",," -> ","
    s = _SPACE_AROUND_PAREN_RE.sub(lambda m: "(" if "(" in m.group(0) else ")", s)

    # Final tidy
    s = s.strip().rstrip(",.;")

    return s


def split_and_clean_targets_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Split and clean the 'Targets' column into a normalized long table.

    Splitting order
    ---------------
    1) '&&'
    2) ','
    3) ';'
    4) ' and ' (lowercase; conservative split)

    After each split, strings are cleaned via clean_target_string,
    and empty rows are removed.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing columns ['dbAMP_ID', 'Targets'].

    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns ['dbAMP_ID', 'Targets'] (deduplicated).
    """
    if "Targets" not in df.columns or "dbAMP_ID" not in df.columns:
        raise KeyError("Input DataFrame must contain columns ['dbAMP_ID', 'Targets'].")

    # 1) split by '&&'
    df1 = (
        df.assign(Targets=df["Targets"].astype(str).str.split("&&"))
        .explode("Targets")
        .reset_index(drop=True)
    )

    # Clean once and 2) split by comma
    df1["Targets"] = df1["Targets"].apply(clean_target_string)
    df2 = (
        df1.assign(Targets=df1["Targets"].str.split(","))
        .explode("Targets")
        .reset_index(drop=True)
    )

    # 3) split by semicolon ';'
    df2["Targets"] = df2["Targets"].apply(clean_target_string)
    df3a = (
        df2.assign(Targets=df2["Targets"].str.split(";"))
        .explode("Targets")
        .reset_index(drop=True)
    )

    # 4) split by ' and '
    def _split_and(x: str):
        x = str(x)
        return x.split(" and ") if " and " in x else [x]

    df3 = (
        df3a.assign(Targets=df3a["Targets"].apply(_split_and))
        .explode("Targets")
        .reset_index(drop=True)
    )

    # Final cleaning and pruning
    df3["Targets"] = df3["Targets"].apply(clean_target_string)
    df3 = df3[df3["Targets"].notna()]
    df3 = df3[df3["Targets"].astype(str).str.strip() != ""]
    df3 = df3[["dbAMP_ID", "Targets"]].drop_duplicates(ignore_index=True)

    return df3


def clean_targets_table(
    input_path: str,
    output_path: str,
    logger: CustomLogger,
) -> None:
    """
    Load AMP data, keep rows with non-empty 'Targets', text-clean and split 'Targets'
    into a normalized long table, then save as CSV. Also auto-exports a unique
    list of Targets to '<base>_unique<ext>'.

    Parameters
    ----------
    input_path : str
        Path to the AMP table (Excel/CSV/TSV supported by loader) with columns ['dbAMP_ID', 'Targets'].
        Typically: data/raw/dbAMP/dbAMP3_pepinfo.xlsx
    output_path : str
        Path to save the cleaned long-format CSV.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Clean Targets into normalized long table")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Ensure output directory
        outdir = os.path.dirname(output_path)
        if outdir and not directory_exists(outdir):
            os.makedirs(outdir)

        # Load input
        required_columns = ["dbAMP_ID", "Targets"]
        df = load_dataframe_by_columns(
            file_path=input_path, required_columns=required_columns, has_header=True
        )

        total = df.shape[0]

        # Keep only non-empty Targets
        df["Targets"] = df["Targets"].apply(
            lambda x: str(x).strip() if pd.notna(x) else x
        )
        df = df[df["Targets"].notna()]
        df = df[df["Targets"].astype(str).str.len() > 0]
        kept = df.shape[0]
        dropped = total - kept

        # Split & clean
        df_long = split_and_clean_targets_df(df)
        final_rows = df_long.shape[0]

        # Save main long-format
        df_long.to_csv(output_path, index=False)

        # Auto-generate unique list path and save
        base, ext = os.path.splitext(output_path)
        unique_path = f"{base}_unique{ext}"
        unique_df = (
            pd.DataFrame({"Targets": df_long["Targets"].dropna().unique()})
            .sort_values(by="Targets", key=lambda s: s.str.casefold(), kind="mergesort")
            .reset_index(drop=True)
        )
        unique_df.to_csv(unique_path, index=False)

        # Summary
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Targets cleaning summary ]\n"
                f"▸ Input records              : {total}\n"
                f"▸ Dropped (Targets NA/empty) : {dropped}\n"
                f"▸ Retained for cleaning      : {kept}\n"
                f"▸ Output rows (long-format)  : {final_rows}\n"
                f"▸ Unique targets             : {unique_df.shape[0]}\n"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message=f"CSV saved:\n'{output_path}'\n'{unique_path}'",
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
        logger.exception("Unexpected error in 'clean_targets_table()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_clean_targets(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Clean AMP 'Targets' from raw dbAMP table and export a normalized long CSV.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance for structured logging.
    **kwargs : dict
        Optional overrides for input/output file paths.

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        input_path = kwargs.get(
            "targets_input_path",
            os.path.join(base_path, "data/raw/dbAMP/dbAMP3_pepinfo.xlsx"),
        )
        output_path = kwargs.get(
            "targets_output_csv",
            os.path.join(base_path, "data/processed/dbAMP/targets_clean.csv"),
        )

        # Validate input file
        if not file_exists(input_path):
            raise FileNotFoundError(f"File not found: '{input_path}'")

        # Execute
        clean_targets_table(
            input_path=input_path,
            output_path=output_path,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_clean_targets()'")
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
