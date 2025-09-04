# pylint: disable=superfluous-parens, line-too-long, too-many-lines, import-error, wrong-import-position, broad-exception-caught, too-many-locals, too-many-arguments, too-many-positional-arguments, too-many-statements
"""
Resolve AMP TaxID via BLAST Alignment to UniProt
"""
# ============================== Standard Library Imports ==============================
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from typing import Optional, Tuple

# ============================== Third-Party Library Imports ==============================
import pandas as pd
from Bio import SeqIO
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
from utils.io_utils import (
    directory_exists,
    file_exists,
    load_dataframe_by_columns,
    write_readme,
)
from utils.log_utils import get_pipeline_completion_message, get_task_completion_message


# ============================== Custom Functions ==============================
def prepare_fasta_for_blast(
    input_path: str, output_dir: str, logger: CustomLogger
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load dbAMP dataset, split by TaxID presence, and export missing-TaxID sequences to FASTA.

    Parameters
    ----------
    input_path : str
        Path to the input dbAMP Excel file.
    output_dir : str
        Directory to save all output files.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        df_tax_na: DataFrame with missing TaxID but valid sequences.
        df_tax_notna: DataFrame with valid TaxID entries.
    """
    logger.info("/ Task: Load dbAMP dataset and split by TaxID presence.")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Ensure output directory exists
        step1_dir = os.path.join(output_dir, "01_split")
        if not directory_exists(dir_path=step1_dir):
            os.makedirs(name=step1_dir)

        # Define output file paths
        df_tax_na_output_path = os.path.join(step1_dir, "na.csv")
        df_tax_notna_output_path = os.path.join(step1_dir, "notna.csv")
        fasta_output_path = os.path.join(step1_dir, "na.fasta")

        # Load dbAMP dataset
        required_columns = ["dbAMP_ID", "Tax", "Targets", "Seq"]
        df = load_dataframe_by_columns(
            file_path=input_path, required_columns=required_columns, has_header=True
        )

        # Split rows by Tax presence
        df_tax_na = df[df["Tax"].isna()].copy()
        df_tax_notna = df[df["Tax"].notna()].copy()
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ The AMPs in dbAMP ]\n"
                f"▸ Number of records: {df.shape[0]}\n"
                f"  - TaxID is 'NA': {df_tax_na.shape[0]}\n"
                f"  - TaxID is 'not NA': {df_tax_notna.shape[0]}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Save CSV files
        df_tax_na.to_csv(path_or_buf=df_tax_na_output_path, index=False)
        df_tax_notna.to_csv(path_or_buf=df_tax_notna_output_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"CSV saved:\n'{df_tax_na_output_path}'\n'{df_tax_notna_output_path}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Write FASTA file
        with open(fasta_output_path, "w", encoding="utf-8") as f:
            for _, row in df_tax_na.iterrows():
                f.write(f">{row['dbAMP_ID']}\n{row['Seq']}\n")
        logger.log_with_borders(
            level=logging.INFO,
            message=f"FASTA saved:\n'{fasta_output_path}'",
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

        return df_tax_na, df_tax_notna

    except Exception:
        logger.exception("Unexpected error in 'prepare_fasta_for_blast()'")
        raise


def clean_tax_string(tax_str: str) -> str:
    """
    Clean a full lineage string by:
    - Removing leading/trailing spaces from each taxonomic level
    - Removing trailing period or space from the entire lineage string

    Parameters
    ----------
    tax_str : str
        A comma-separated taxonomic lineage string.

    Returns
    -------
    str
        Cleaned lineage string.
    """
    cleaned = ",".join([level.strip() for level in str(tax_str).split(",")])
    cleaned = cleaned.rstrip(". ").strip()
    return cleaned


def split_and_clean_tax_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Split the 'Tax' column by '&&' and clean each lineage string.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least a 'Tax' column.

    Returns
    -------
    pd.DataFrame
        Processed DataFrame with 'Tax' split into multiple rows and cleaned.
    """
    # Split 'Tax' column into multiple rows using '&&' as a delimiter
    df_split = (
        df.assign(Tax=df["Tax"].astype(str).str.split("&&"))
        .explode("Tax")
        .reset_index(drop=True)
    )
    # Clean each lineage string
    df_split["Tax"] = df_split["Tax"].apply(clean_tax_string)
    # Remove rows with empty or null 'Tax' values
    df_split = df_split[df_split["Tax"].notna()]
    df_split = df_split[df_split["Tax"].str.len() > 0]

    return df_split


def get_taxonomy_lineage_from_name(
    name: str, ncbi: NCBITaxa
) -> Tuple[Optional[dict], Optional[str], int]:
    """
    Query NCBI with a taxon name and return lineage dict, taxids string, and count.

    Parameters
    ----------
    name : str
        Taxon name to search in NCBI.
    ncbi : NCBITaxa
        Initialized ETE3 NCBITaxa instance.

    Returns
    -------
    tuple
        lineage_dict (dict or None): Mapping {rank: name}, None if not found.
        matched_taxids (str or None): Comma-separated string of matching TaxIDs.
        match_count (int): Number of matching TaxIDs.
    """
    try:
        taxids = ncbi.get_name_translator([name])[name]
        matched_taxids = ",".join(map(str, taxids))
        match_count = len(taxids)

        taxid = taxids[0]
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        lineage_dict = {ranks[t]: names[t] for t in lineage if ranks[t] != "no rank"}  # type: ignore
        return lineage_dict, matched_taxids, match_count

    except Exception:
        return None, None, 0


def check_and_correct_taxonomy(
    original_tax: str, ncbi: NCBITaxa
) -> Tuple[Optional[str], Optional[str], int]:
    """
    Compare the original lineage with the NCBI lineage and record changes.

    Parameters
    ----------
    original_tax : str
        Original taxonomy string.
    ncbi : NCBITaxa
        Initialized ETE3 NCBITaxa instance.

    Returns
    -------
    tuple
        checked_tax (str or None): Corrected lineage from NCBI or None if not found.
        matched_taxids (str or None): Comma-separated TaxIDs.
        match_count (int): Number of matching TaxIDs.
    """
    last_level = original_tax.split(",")[-1].strip()
    ncbi_lineage, taxid_str, count = get_taxonomy_lineage_from_name(
        name=last_level, ncbi=ncbi
    )

    if ncbi_lineage is None:
        return None, taxid_str, count

    ncbi_tax_str = ",".join(ncbi_lineage.values())
    return ncbi_tax_str, taxid_str, count


def annotate_taxonomy_with_ncbi(output_dir: str, logger: CustomLogger) -> None:
    """
    Annotate a DataFrame taxonomy column using NCBI lineage check.

    Parameters
    ----------
    output_dir : str
        Directory to save all output files.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Annotating taxonomy with NCBI lineage")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    ncbi = NCBITaxa()
    start_time = time.time()

    try:
        # Ensure output directory exists
        step1_dir = os.path.join(output_dir, "01_split")
        step2_dir = os.path.join(output_dir, "02_tax_check")
        if not directory_exists(dir_path=step1_dir):
            os.makedirs(name=step1_dir)
        if not directory_exists(dir_path=step2_dir):
            os.makedirs(name=step2_dir)

        # Define output file paths
        amp_source_path = os.path.join(step1_dir, "notna.csv")
        df_no_match_output_path = os.path.join(step2_dir, "notna_nomatch.csv")
        nomatch_unique_path = os.path.join(step2_dir, "notna_nomatch_unique.csv")
        df_single_match_output_path = os.path.join(step2_dir, "notna_single.csv")
        df_multi_match_output_path = os.path.join(step2_dir, "notna_multi.csv")

        # Load Dataframe
        df_amp = load_dataframe_by_columns(file_path=amp_source_path)

        # Split and clean 'Tax' column
        df_amp = split_and_clean_tax_df(df_amp)

        # Initialize in-memory cache to avoid redundant NCBI lookups
        taxonomy_cache: dict[str, tuple[str, list[int], int]] = {}

        def cached_check_taxonomy(x: str) -> pd.Series:
            """
            Wrapper to cache the results of NCBI taxonomy matching for a given lineage string.

            Parameters
            ----------
            x : str
                The raw lineage string from 'Tax' column.

            Returns
            -------
            pd.Series
                A Series of (Tax_Checked, Matched_TaxIDs, TaxID_Count)
            """
            x_clean = str(x).strip()
            if x_clean in taxonomy_cache:
                return pd.Series(taxonomy_cache[x_clean])
            result = check_and_correct_taxonomy(original_tax=x_clean, ncbi=ncbi)
            taxonomy_cache[x_clean] = result  # type: ignore
            return pd.Series(result)

        # Apply taxonomy checking with cache acceleration
        df_out = df_amp.copy()
        df_out[["Tax_Checked", "Matched_TaxIDs", "TaxID_Count"]] = df_out["Tax"].apply(
            cached_check_taxonomy
        )
        df_out["TaxID_Count"] = pd.to_numeric(
            df_out["TaxID_Count"], errors="coerce"
        ).astype("Int64")

        # Categorize results
        df_no_match = df_out[df_out["TaxID_Count"] == 0]
        df_single_match = df_out[df_out["TaxID_Count"] == 1]
        df_multi_match = df_out[df_out["TaxID_Count"] > 1]

        # Extract unique unresolved lineage strings
        unique_tax_values = (
            pd.DataFrame(df_no_match["Tax"].dropna().unique(), columns=["Tax"])
            .sort_values("Tax", kind="mergesort")
            .reset_index(drop=True)
        )

        # Log statistics for each category
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"[ The AMPs' TaxID is not NA in dbAMP after NCBI annotation ]\n"
                f"▸ Total records processed: {df_out.shape[0]}\n"
                f"  - TaxID_Count == 0 (no match): {df_no_match.shape[0]} (Unique lineage count: {unique_tax_values.shape[0]})\n"
                f"  - TaxID_Count == 1 (single match): {df_single_match.shape[0]}\n"
                f"  - TaxID_Count > 1 (multiple matches): {df_multi_match.shape[0]}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Save each category to separate CSV files
        df_no_match.to_csv(path_or_buf=df_no_match_output_path, index=False)
        df_single_match.to_csv(path_or_buf=df_single_match_output_path, index=False)
        df_multi_match.to_csv(path_or_buf=df_multi_match_output_path, index=False)
        unique_tax_values.to_csv(nomatch_unique_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"CSV saved:\n'{df_no_match_output_path}'\n'{nomatch_unique_path}'\n'{df_single_match_output_path}'\n'{df_multi_match_output_path}'",
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
        logger.exception("Unexpected error in 'annotate_taxonomy_with_ncbi()'")
        raise


def download_uniprot_sprot(
    resource_dir: str, force_download: bool, logger: CustomLogger
) -> None:
    """
    Download and decompress the UniProt Swiss-Prot FASTA file if not already present,
    or force re-download if specified.

    Parameters
    ----------
    resource_dir : str
        Base directory where external resources (e.g., UniProt FASTA) are stored.
    force_download : bool
        If True, force re-download even if file already exists.
    logger : CustomLogger
        Logger instance to log progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Download and decompress UniProt Swiss-Prot FASTA")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Define variables
        url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        fasta_gz_path = os.path.join(resource_dir, "uniprot_sprot.fasta.gz")
        fasta_path = os.path.join(resource_dir, "uniprot_sprot.fasta")

        # If file already exists and no force flag
        if file_exists(file_path=fasta_path) and not force_download:
            logger.log_with_borders(
                level=logging.INFO,
                message="UniProt FASTA already exists. Skipping download.",
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        else:
            if force_download:
                logger.log_with_borders(
                    level=logging.INFO,
                    message="'Force re-download enabled'. Removing existing UniProt directory.",
                    border="|",
                    length=120,
                )
                logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
                if directory_exists(dir_path=resource_dir):
                    shutil.rmtree(resource_dir)

            # Re-create directory
            os.makedirs(resource_dir, exist_ok=True)

            # Download .gz file
            logger.log_with_borders(
                level=logging.INFO,
                message=f"Downloading 'UniProt Swiss-Prot' from:\n{url}",
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            urllib.request.urlretrieve(url, fasta_gz_path)

            # Decompress .gz to .fasta
            with gzip.open(fasta_gz_path, "rb") as f_in:
                with open(fasta_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(fasta_gz_path)
            logger.log_with_borders(
                level=logging.INFO,
                message=f"FASTA decompressed:\n'{fasta_path}'",
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
        logger.exception("Unexpected error in 'download_uniprot_sprot()'")
        raise


def make_blast_database(
    resource_dir: str, force_rebuild: bool, logger: CustomLogger
) -> None:
    """
    Create a BLAST protein database from the UniProt FASTA file using makeblastdb.

    Parameters
    ----------
    resource_dir : str
        Base directory where external resources (e.g., UniProt FASTA) are stored.
    force_rebuild : bool
        Whether to force rebuild the BLAST database even if it already exists.
    logger : CustomLogger
        Logger instance for logging.

    Returns
    -------
    None
    """
    logger.info("/ Task: Create BLAST database from UniProt FASTA")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Define variables
        fasta_path = os.path.join(resource_dir, "uniprot_sprot.fasta")
        db_prefix = os.path.splitext(fasta_path)[0]
        psq_file = db_prefix + ".psq"
        blast_db_extensions = [
            ".phr",
            ".pin",
            ".psq",
            ".pdb",
            ".pog",
            ".pos",
            ".pjs",
            ".pot",
            ".ptf",
            ".pto",
        ]

        # Check file exists
        if not file_exists(file_path=fasta_path):
            raise FileNotFoundError(f"File not found: '{fasta_path}'")

        # If file already exists and no force flag
        if file_exists(file_path=psq_file) and not force_rebuild:
            logger.log_with_borders(
                level=logging.INFO,
                message="BLAST database already exists. Skipping makeblastdb.",
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        else:
            if force_rebuild:
                logger.log_with_borders(
                    level=logging.INFO,
                    message="'Force rebuild enabled'. Removing existing BLAST database files.",
                    border="|",
                    length=120,
                )
                logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
                for ext in blast_db_extensions:
                    file_to_remove = db_prefix + ext
                    if file_exists(file_path=file_to_remove):
                        os.remove(file_to_remove)

            # Run makeblastdb
            logger.log_with_borders(
                level=logging.INFO,
                message="Running makeblastdb to build BLAST protein database.",
                border="|",
                length=120,
            )
            cmd = [
                "makeblastdb",
                "-in",
                fasta_path,
                "-dbtype",
                "prot",
                "-parse_seqids",
                "-out",
                db_prefix,
                "-title",
                "uniprot_sprot",
            ]
            subprocess.run(cmd, check=True)
            logger.log_with_borders(
                level=logging.INFO,
                message="BLAST database created successfully.",
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
        logger.exception("Unexpected error in 'make_blast_database()'")
        raise


def split_query_fasta_for_blast(
    input_fasta: str,
    clean_fasta: str,
    invalid_fasta: str,
    logger: CustomLogger,
) -> tuple[int, int]:
    """
    Split query FASTA sequences into valid and invalid sets.

    Valid sequences:
      - Only contain allowed amino-acid codes
      - Not composed entirely of ambiguous/unknown residues (X, B, Z)

    Invalid sequences:
      - Contain any illegal characters
      - OR are entirely composed of ambiguous/unknown residues

    Parameters
    ----------
    input_fasta : str
        Path to the original query FASTA file.
    clean_fasta : str
        Output path for valid sequences (to be used for BLAST).
    invalid_fasta : str
        Output path for invalid sequences.
    logger : CustomLogger
        Logger instance for recording process details.

    Returns
    -------
    tuple[int, int]
        (num_valid, num_invalid)
    """
    allowed = set("ACDEFGHIKLMNPQRSTVWYXBZU")  # note: 'O' is disallowed
    ambiguous = set("XBZ")  # ambiguity codes

    clean_records, invalid_records = [], []
    invalid_meta: list[tuple[str, str]] = []  # (query_id, reason)

    for rec in SeqIO.parse(input_fasta, "fasta"):
        seq = str(rec.seq).upper()
        qid = rec.id

        if not seq:  # empty
            invalid_records.append(rec)
            invalid_meta.append((qid, "empty"))
        elif any(ch not in allowed for ch in seq):  # illegal chars
            illegal = sorted({ch for ch in seq if ch not in allowed})
            invalid_records.append(rec)
            invalid_meta.append((qid, f"illegal_char:{''.join(illegal)}"))
        elif all(ch in ambiguous for ch in seq):  # all ambiguous
            invalid_records.append(rec)
            invalid_meta.append((qid, "all_ambiguous(XBZ)"))
        else:  # valid
            clean_records.append(rec)

    # Write valid sequences
    with open(clean_fasta, "w", encoding="utf-8") as f:
        SeqIO.write(clean_records, f, "fasta")

    # Write invalid sequences + CSV
    invalid_csv_path = os.path.splitext(invalid_fasta)[0] + ".csv"
    if invalid_records:
        with open(invalid_fasta, "w", encoding="utf-8") as f:
            SeqIO.write(invalid_records, f, "fasta")
        pd.DataFrame(invalid_meta, columns=["query_id", "reason"]).to_csv(
            invalid_csv_path, index=False
        )
    else:
        # Always create CSV for consistency
        pd.DataFrame(columns=["query_id", "reason"]).to_csv(
            invalid_csv_path, index=False
        )

    # Log summary
    reason_counts = {}
    for _, r in invalid_meta:
        reason_counts[r] = reason_counts.get(r, 0) + 1
    breakdown = "\n".join(
        f"  - '{k}': {v}" for k, v in sorted(reason_counts.items(), key=lambda x: -x[1])
    )
    logger.log_with_borders(
        level=logging.INFO,
        message=(
            "[ Query sequence QC in 05_blast ]\n"
            f"▸ Clean sequences: {len(clean_records)}\n"
            f"▸ Invalid sequences: {len(invalid_records)}"
            + (f"\n▸ Invalid reasons:\n{breakdown}" if reason_counts else "")
        ),
        border="|",
        length=120,
    )
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    return len(clean_records), len(invalid_records)


def extract_fields_from_stitle(stitle: str) -> pd.Series:
    """
    Extract structured fields from the BLAST stitle field.

    Parameters
    ----------
    stitle : str
        The stitle field from BLAST output (e.g., "Protein name OS=Organism OX=1234 ...")

    Returns
    -------
    pd.Series
        A pandas Series with extracted fields:
        - protein_name
        - organism
        - taxid
        - gene
        - PE
        - SV
    """
    pattern = r"(.*?) OS=(.*?) OX=(\d+)(?: GN=(\S+))? PE=(\d) SV=(\d)"
    match = re.match(pattern, stitle)
    if match:
        return pd.Series(
            {
                "protein_name": match.group(1).strip(),
                "organism": match.group(2).strip(),
                "taxid": int(match.group(3)),
                "gene": match.group(4) if match.group(4) else None,
                "PE": int(match.group(5)),
                "SV": int(match.group(6)),
            }
        )

    return pd.Series(
        {
            "protein_name": None,
            "organism": None,
            "taxid": None,
            "gene": None,
            "PE": None,
            "SV": None,
        }
    )


def blastp_alignment(
    output_dir: str,
    resource_dir: str,
    logger: CustomLogger,
    evalue: float = 1e-5,
    max_target_seqs: int = 1,
    max_hsps: int = 1,
    pident_min: float = 30.0,
    qcovs_min: float = 60.0,
    num_threads: int = 4,
    pe_max: int = 2,
    outfmt: str = "6 qseqid sseqid pident length evalue bitscore qcovs stitle",
) -> bool:
    """
    Run BLASTP alignment for the given query FASTA against the specified BLAST database.

    Parameters
    ----------
    output_dir : str
        Directory to save pipeline step outputs.
    resource_dir : str
        Directory containing UniProt FASTA and BLAST database files.
    logger : CustomLogger
        Project-specific logger for structured logging.
    evalue : float, optional
        E-value cutoff for BLAST (default: 1e-5).
    max_target_seqs : int, optional
        Max number of target sequences to return per query (default: 1).
    max_hsps : int, optional
        Max number of HSPs (High-Scoring Segment Pairs) per subject (default: 1).
    pident_min : float, optional
        Minimum percent identity threshold (default: 30.0).
    qcovs_min : float, optional
        Minimum query coverage threshold (default: 60.0).
    outfmt : str, optional
        BLAST output format string (default includes fields needed for filtering).
    num_threads : int, optional
        Number of CPU threads to use for BLAST (default: 4).

    Returns
    -------
    bool
        True if final filtered BLAST results contain at least one unique TaxID hit per query, else False.
    """
    logger.info("/ Task: Run 'BLASTP' alignment")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Validate PE max value
        if not isinstance(pe_max, int) or not (1 <= pe_max <= 5):
            logger.log_with_borders(
                level=logging.ERROR,
                message=(
                    f"Invalid PE max value: {pe_max}. Must be an integer between 1 and 5.\n"
                    "Please check the argument 'pe_max' in your pipeline call."
                ),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message="\n".join(get_task_completion_message(start_time)),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            return False

        # Ensure output directory exists
        step1_dir = os.path.join(output_dir, "01_split")
        step5_dir = os.path.join(output_dir, "05_blast")
        if not directory_exists(dir_path=step1_dir):
            os.makedirs(name=step1_dir)
        if not directory_exists(dir_path=step5_dir):
            os.makedirs(name=step5_dir)

        # Define input and output file paths
        query_fasta_path = os.path.join(step1_dir, "na.fasta")
        clean_query_fasta_path = os.path.join(step5_dir, "na_clean.fasta")
        invalid_query_fasta_path = os.path.join(step5_dir, "na_invalid.fasta")
        uniprot_fasta_path = os.path.join(resource_dir, "uniprot_sprot.fasta")
        blast_db_prefix = os.path.splitext(uniprot_fasta_path)[0]
        psq_file = blast_db_prefix + ".psq"
        blast_output_tsv_path = os.path.join(step5_dir, "na.tsv")
        blast_filt_csv_path = os.path.join(step5_dir, "na_filt.csv")
        blast_final_csv_path = os.path.join(step5_dir, "na_final.csv")

        # Validate files
        if not file_exists(file_path=query_fasta_path):
            raise FileNotFoundError(f"File not found: '{query_fasta_path}'")
        if not file_exists(file_path=psq_file):
            raise FileNotFoundError(f"File not found: '{psq_file}'")

        # Split query FASTA into clean/invalid in 05_blast
        num_clean, _ = split_query_fasta_for_blast(
            input_fasta=query_fasta_path,
            clean_fasta=clean_query_fasta_path,
            invalid_fasta=invalid_query_fasta_path,
            logger=logger,
        )
        if num_clean == 0:
            logger.log_with_borders(
                level=logging.WARNING,
                message=(
                    "No valid query sequences remain after QC. "
                    f"Invalid sequences have been written to: '{invalid_query_fasta_path}'. "
                    "Skipping BLASTP."
                ),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message="\n".join(get_task_completion_message(start_time)),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            return False

        # Use the cleaned FASTA for BLAST
        query_fasta_path = clean_query_fasta_path

        # Log parameters for reproducibility
        num_threads = max(1, min(num_threads, os.cpu_count() or 1))
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "Running BLASTP with the following parameters:\n"
                f"- 'Query': '{query_fasta_path}'\n"
                f"- 'Database': '{blast_db_prefix}'\n"
                f"- 'E-value': {evalue}\n"
                f"- 'Max target seqs': {max_target_seqs}\n"
                f"- 'Max HSPs per subject': {max_hsps}\n"
                f"- 'pident_min': {pident_min}\n"
                f"- 'qcovs_min': {qcovs_min}\n"
                f"- 'PE max (UniProt)': {pe_max if pe_max is not None else 'disabled'}\n"
                f"- 'num_threads': {num_threads}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Run BLASTP
        cmd = [
            "blastp",
            "-query",
            query_fasta_path,
            "-db",
            blast_db_prefix,
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            str(max_target_seqs),
            "-max_hsps",
            str(max_hsps),
            "-outfmt",
            outfmt,
            "-out",
            blast_output_tsv_path,
            "-num_threads",
            str(num_threads),
        ]
        subprocess.run(cmd, check=True)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Raw TSV saved:\n'{blast_output_tsv_path}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # -------------------- Parse BLAST TSV --------------------
        df = load_dataframe_by_columns(
            file_path=blast_output_tsv_path, has_header=False
        )

        # Check if BLASTP returned any alignments at all
        if df.empty:
            logger.log_with_borders(
                level=logging.WARNING,
                message="BLASTP completed but returned no alignments. Skipping.",
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message="\n".join(get_task_completion_message(start_time)),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            return False

        # Ensure numeric types for filtering
        df.columns = [
            "dbAMP_ID",
            "subject_id",
            "pident",
            "align_len",
            "evalue",
            "bitscore",
            "qcovs",
            "stitle",
        ]
        for col in ["pident", "align_len", "evalue", "bitscore", "qcovs"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        # Count query statistics before filtering
        hit_queries = df["dbAMP_ID"].nunique()
        total_queries = sum(1 for _ in SeqIO.parse(query_fasta_path, "fasta"))
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Raw BLAST alignment statistics ]\n"
                f"▸ Total query sequences: {total_queries}\n"
                f"  - Queries with hits: {hit_queries}\n"
                f"  - Queries without hits: {total_queries - hit_queries}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # -------------------- Extract stitle fields --------------------
        parsed_df = df.join(df["stitle"].apply(extract_fields_from_stitle))
        parsed_df["taxid"] = pd.to_numeric(parsed_df["taxid"], errors="coerce").astype(
            "Int64"
        )

        # -------------------- Apply filtering --------------------
        filt_df = parsed_df[
            (parsed_df["pident"] >= pident_min) & (parsed_df["qcovs"] >= qcovs_min)
        ].copy()
        filt_df = filt_df[(filt_df["PE"].notna()) & (filt_df["PE"] <= pe_max)]
        filt_df.to_csv(blast_filt_csv_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ After applying thresholds ]\n"
                f"▸ Retained hits: {filt_df.shape[0]} "
                f"(pident >= {pident_min}, qcovs >= {qcovs_min})"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # -------------------- Select unique TaxID per query --------------------
        final_records = []
        for _, group in filt_df.groupby("dbAMP_ID"):
            unique_taxids = group["taxid"].dropna().unique()
            if len(unique_taxids) == 1:
                # Keep only queries that map to exactly one unique TaxID
                best_hit = group.sort_values(
                    by=["bitscore", "pident", "qcovs"], ascending=[False, False, False]
                ).iloc[0]
                final_records.append(best_hit)
        final_df = pd.DataFrame(final_records)

        # Check if any hits passed filtering AND mapped to exactly one TaxID
        if final_df.empty:
            logger.log_with_borders(
                level=logging.WARNING,
                message=(
                    "No queries retained after filtering for unique TaxID. "
                    "BLASTP result cannot be used for downstream steps."
                ),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            logger.log_with_borders(
                level=logging.INFO,
                message="\n".join(get_task_completion_message(start_time)),
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            return False
        final_df.to_csv(blast_final_csv_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Final BLAST alignment selection ]\n"
                f"▸ Queries retained with unique TaxID: {final_df['dbAMP_ID'].nunique()}\n"
                f"▸ Final hits saved to:\n'{blast_final_csv_path}'"
            ),
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

        return True

    except Exception:
        logger.exception("Unexpected error during 'blastp_alignment()'")
        raise


def get_taxonomy_lineage_from_taxid(taxid: int, ncbi: NCBITaxa) -> Optional[str]:
    """
    Retrieve full NCBI taxonomy lineage as comma-separated names, excluding 'no rank'.

    Parameters
    ----------
    taxid : int
        NCBI TaxID to query.
    ncbi : NCBITaxa
        Initialized ETE3 NCBITaxa instance.

    Returns
    -------
    str or None
        Comma-separated taxonomy lineage names, or None if failed.
    """
    try:
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)

        name_list = [names[t] for t in lineage if ranks[t] != "no rank" and t in names]  # type: ignore

        return ",".join(name_list)

    except Exception:
        return None


def integrate_amp_taxonomy_info(output_dir: str, logger: CustomLogger) -> None:
    """
    Merge BLAST output (with TaxID) and AMP metadata, append full taxonomy lineage.

    Parameters
    ----------
    output_dir : str
        Directory to save all output files.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info(
        "/ Task: Merge BLAST TaxID results with AMP metadata and resolve NCBI lineages"
    )
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    ncbi = NCBITaxa()
    start_time = time.time()

    try:
        # Ensure output directory exists
        step1_dir = os.path.join(output_dir, "01_split")
        step5_dir = os.path.join(output_dir, "05_blast")
        step6_dir = os.path.join(output_dir, "06_merge")
        if not directory_exists(dir_path=step1_dir):
            os.makedirs(name=step1_dir)
        if not directory_exists(dir_path=step5_dir):
            os.makedirs(name=step5_dir)
        if not directory_exists(dir_path=step6_dir):
            os.makedirs(name=step6_dir)

        # Define output file paths
        parsed_blast_path = os.path.join(step5_dir, "na_final.csv")
        amp_source_path = os.path.join(step1_dir, "na.csv")
        output_path = os.path.join(step6_dir, "na_with_blast.csv")

        # Load Dataframe
        df_blast = load_dataframe_by_columns(file_path=parsed_blast_path)
        df_amp = load_dataframe_by_columns(file_path=amp_source_path)

        # Merge Targets & Seq by dbAMP_ID
        df_merged = df_blast.merge(df_amp, on="dbAMP_ID", how="left")

        # Resolve taxonomy lineage from taxid
        logger.log_with_borders(
            level=logging.INFO,
            message="Resolving NCBI lineage from TaxID.",
            border="|",
            length=120,
        )
        df_merged["Tax"] = df_merged["taxid"].apply(
            lambda tid: (
                get_taxonomy_lineage_from_taxid(int(tid), ncbi)
                if not pd.isna(tid)
                else None
            )
        )

        # Log save location
        df_final = df_merged[["dbAMP_ID", "Tax", "Targets", "Seq"]]
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Number of merged AMP entries: {df_final.shape[0]}",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        df_final.to_csv(output_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Final CSV saved:\n'{output_path}'",
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
        logger.exception("Unexpected error in 'integrate_amp_taxonomy_info()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_resolve_taxid_by_blast(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Use UniProt BLAST results to resolve missing TaxIDs in dbAMP3_pepinfo.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance for tracking progress.

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        # -------------------- Retrieve input parameters (CLI or default) --------------------
        input_path = kwargs.get(
            "resolve_input_dbamp",
            os.path.join(base_path, "data/raw/dbAMP/dbAMP3_pepinfo.xlsx"),
        )
        force_download = kwargs.get("resolve_force_download", False)
        force_build = kwargs.get("resolve_force_build", False)
        resolve_output_dir = kwargs.get(
            "resolve_output_dir", os.path.join(base_path, "data/interim/blast")
        )
        resource_dir = kwargs.get(
            "resource_dir", os.path.join(base_path, "data/external/uniprot")
        )
        resolve_blast_evalue = kwargs.get("resolve_blast_evalue", 1e-5)
        resolve_blast_max_hits = kwargs.get("resolve_blast_max_hits", 25)
        resolve_blast_max_hsps = kwargs.get("resolve_blast_max_hsps", 1)
        resolve_blast_pident_min = kwargs.get("resolve_blast_pident_min", 30.0)
        resolve_blast_qcovs_min = kwargs.get("resolve_blast_qcovs_min", 60.0)
        resolve_blast_threads = kwargs.get("resolve_blast_threads", 4)
        resolve_blast_pe_max = kwargs.get("resolve_blast_pe_max", 2)

        # Check folders exist
        if not directory_exists(dir_path=resolve_output_dir):
            os.makedirs(resolve_output_dir)
        if not directory_exists(dir_path=resource_dir):
            os.makedirs(resource_dir)

        # -------------------- Pipeline Execution --------------------
        # Step 1: Load dbAMP file, split by whether TaxID is present, and export the missing ones to FASTA.
        prepare_fasta_for_blast(
            input_path=input_path, output_dir=resolve_output_dir, logger=logger
        )
        write_readme(
            os.path.join(resolve_output_dir, "01_split"),
            "Step 1 - Input Split",
            "Splits dbAMP dataset into NA/notNA and exports NA to FASTA.",
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 2: Annotate and check taxonomy for AMPs with existing Tax (not NA)
        annotate_taxonomy_with_ncbi(
            output_dir=resolve_output_dir,
            logger=logger,
        )
        write_readme(
            os.path.join(resolve_output_dir, "02_tax_check"),
            "Step 2 - Taxonomy Check",
            "Checks existing Tax entries against NCBI lineage.",
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 3: Download and extract UniProt Swiss-Prot FASTA file.
        download_uniprot_sprot(
            resource_dir=resource_dir, force_download=force_download, logger=logger
        )
        write_readme(
            os.path.join(resolve_output_dir, "03_uniprot"),
            "Step 3 - UniProt Download",
            "Downloads and decompresses UniProt Swiss-Prot FASTA.",
            {
                "Output Location": resource_dir,
                "URL": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            },
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 4: Build BLAST database using makeblastdb.
        make_blast_database(
            resource_dir=resource_dir, force_rebuild=force_build, logger=logger
        )
        write_readme(
            os.path.join(resolve_output_dir, "04_blastdb"),
            "Step 4 - Make BLAST Database",
            "Builds BLAST protein DB from UniProt FASTA.",
            {"Output Location": resource_dir},
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 5: Run BLASTP using query FASTA and UniProt BLAST DB.
        blast_success = blastp_alignment(
            output_dir=resolve_output_dir,
            resource_dir=resource_dir,
            logger=logger,
            evalue=resolve_blast_evalue,
            max_target_seqs=resolve_blast_max_hits,
            max_hsps=resolve_blast_max_hsps,
            pident_min=resolve_blast_pident_min,
            qcovs_min=resolve_blast_qcovs_min,
            num_threads=resolve_blast_threads,
            pe_max=resolve_blast_pe_max,
        )
        write_readme(
            os.path.join(resolve_output_dir, "05_blast"),
            "Step 5 - BLASTP Alignment",
            "Runs BLASTP for NA entries against UniProt DB.",
        )
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 6: Merge BLAST result with AMP metadata and resolve full NCBI taxonomy lineages.
        if blast_success:
            integrate_amp_taxonomy_info(
                output_dir=resolve_output_dir,
                logger=logger,
            )
            write_readme(
                os.path.join(resolve_output_dir, "06_merge"),
                "Step 6 - Final Merge",
                "Merges BLAST TaxIDs with AMP metadata and resolves full lineage.",
            )
            logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_resolve_taxid_by_blast()'")
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
