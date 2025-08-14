# pylint: disable=import-error, wrong-import-position, line-too-long, broad-exception-caught
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
from ete3 import NCBITaxa

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
LOGGING_PATH = os.path.join(BASE_PATH, "src/utils/logging_toolkit/src/python")
print(BASE_PATH)

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


# ============================== Custom Functions ==============================
def prepare_fasta_for_blast(
    input_path: str, logger: CustomLogger
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load dbAMP dataset, split by TaxID presence, and export missing-TaxID sequences to FASTA.

    Parameters
    ----------
    input_path : str
        Path to the input dbAMP Excel file.
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
        # Define variables
        interim_dir_path = os.path.join(BASE_PATH, "data/interim")
        df_tax_na_output_path = os.path.join(interim_dir_path, "dbAMP_tax_NA.csv")
        df_tax_notna_output_path = os.path.join(interim_dir_path, "dbAMP_tax_notNA.csv")
        fasta_output_path = os.path.join(interim_dir_path, "dbAMP_tax_NA.fasta")
        required_columns = ["dbAMP_ID", "Tax", "Targets", "Seq"]

        # Load dbAMP dataset
        df = load_dataframe_by_columns(
            file_path=input_path, required_columns=required_columns, has_header=True
        )

        # Split rows by Tax presence and reoport statistics
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
        df_tax_na.to_csv(path_or_buf=df_tax_na_output_path, index=False)
        df_tax_notna.to_csv(path_or_buf=df_tax_notna_output_path, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"CSV saved:\n'{df_tax_na_output_path}'\n'{df_tax_notna_output_path}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Write a FASTA file
        if not directory_exists(dir_path=interim_dir_path):
            os.makedirs(name=interim_dir_path)
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

        # Final summary block
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


def download_uniprot_sprot(force_download: bool, logger: CustomLogger) -> None:
    """
    Download and decompress the UniProt Swiss-Prot FASTA file if not already present,
    or force re-download if specified.

    Parameters
    ----------
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
        uniprot_dir = os.path.join(BASE_PATH, "data/external/uniprot")
        fasta_gz_path = os.path.join(uniprot_dir, "uniprot_sprot.fasta.gz")
        fasta_path = os.path.join(uniprot_dir, "uniprot_sprot.fasta")

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
                if directory_exists(dir_path=uniprot_dir):
                    shutil.rmtree(uniprot_dir)

            # Re-create directory
            os.makedirs(uniprot_dir, exist_ok=True)

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

        # Final summary block
        summary_lines = get_task_completion_message(start_time)
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(summary_lines),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'download_uniprot_sprot()'")
        raise


def make_blast_database(force_rebuild: bool, logger: CustomLogger) -> None:
    """
    Create a BLAST protein database from the UniProt FASTA file using makeblastdb,
    with optional force rebuild.

    Parameters
    ----------
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
        uniprot_dir = os.path.join(BASE_PATH, "data/external/uniprot")
        fasta_path = os.path.join(uniprot_dir, "uniprot_sprot.fasta")
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

        # Final summary block
        summary_lines = get_task_completion_message(start_time)
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(summary_lines),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    except Exception:
        logger.exception("Unexpected error in 'make_blast_database()'")
        raise


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
    output_path: str,
    logger: CustomLogger,
    evalue: float = 1e-5,
    max_target_seqs: int = 1,
    outfmt: str = "6 qseqid sseqid pident length evalue stitle",
) -> None:
    """
    Run BLASTP alignment for the given query FASTA against the specified BLAST database.

    Parameters
    ----------
    output_path : str
        Path to save the BLASTP result output (.tsv).
    logger : CustomLogger
        Logger instance for logging.
    evalue : float, optional
        E-value threshold for BLASTP (default: 1e-5).
    max_target_seqs : int, optional
        Max number of target sequences to report (default: 1).
    outfmt : str, optional
        Output format string for BLAST (default: "6 qseqid sseqid pident length evalue stitle").

    Returns
    -------
    None
    """
    logger.info("/ Task: Run 'BLASTP' alignment")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Define variables
        interim_dir_path = os.path.join(BASE_PATH, "data/interim")
        uniprot_dir = os.path.join(BASE_PATH, "data/external/uniprot")
        query_fasta_path = os.path.join(interim_dir_path, "dbAMP_tax_NA.fasta")
        uniprot_fasta_path = os.path.join(uniprot_dir, "uniprot_sprot.fasta")
        blast_db_prefix = os.path.splitext(uniprot_fasta_path)[0]
        psq_file = blast_db_prefix + ".psq"

        # Validate files
        if not file_exists(file_path=query_fasta_path):
            raise FileNotFoundError(f"File not found: '{query_fasta_path}'")
        if not file_exists(file_path=psq_file):
            raise FileNotFoundError(f"File not found: '{psq_file}'")

        # Log parameters
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "Running BLASTP with the following parameters:\n"
                f"- 'Query': '{query_fasta_path}'\n"
                f"- 'Database': '{blast_db_prefix}'\n"
                f"- 'E-value': {evalue}\n"
                f"- 'Max target seqs': {max_target_seqs}"
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
            "-outfmt",
            outfmt,
            "-out",
            output_path,
        ]
        subprocess.run(cmd, check=True)

        logger.log_with_borders(
            level=logging.INFO,
            message="'BLASTP' alignment completed successfully.",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # -------------------- Parse BLAST TSV --------------------
        df = load_dataframe_by_columns(file_path=output_path, has_header=False)
        df.columns = [
            "dbAMP_ID",
            "subject_id",
            "pident",
            "align_len",
            "evalue",
            "stitle",
        ]

        parsed_df = df.join(df["stitle"].apply(extract_fields_from_stitle))
        parsed_df.to_csv(output_path, sep="\t", index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Parsed TSV saved:\n'{output_path}'",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Final summary block
        logger.log_with_borders(
            level=logging.INFO,
            message="\n".join(get_task_completion_message(start_time)),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

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


def integrate_amp_taxonomy_info(
    parsed_blast_path: str, output_path: str, logger: CustomLogger
) -> None:
    """
    Merge BLAST output (with TaxID) and AMP metadata, append full taxonomy lineage.

    Parameters
    ----------
    parsed_blast_path : str
        Path to the parsed BLAST result file with taxid.
    output_path : str
        Path to save the final merged output.
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
        # Define variables
        interim_dir_path = os.path.join(BASE_PATH, "data/interim")
        amp_source_path = os.path.join(interim_dir_path, "dbAMP_tax_NA.csv")

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

        # Final summary block
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
        input_path = kwargs.get("resolve_input_dbamp") or os.path.join(
            base_path, "data/raw/dbAMP/dbAMP3_pepinfo.xlsx"
        )
        force_download = kwargs.get("resolve_force_download", False)
        force_build = kwargs.get("resolve_force_build", False)
        blast_result_path = kwargs.get("resolve_blast_result_path") or os.path.join(
            base_path, "data/interim/blast_results.tsv"
        )
        resolve_output_path = kwargs.get("resolve_output_path") or os.path.join(
            base_path, "data/processed/dbAMP_with_blast_tax.csv"
        )

        # -------------------- Pipeline Execution --------------------
        # Step 1: Load dbAMP file, split by whether TaxID is present, and export the missing ones to FASTA.
        # Output: data/interim/dbAMP_tax_NA.fasta
        prepare_fasta_for_blast(input_path=input_path, logger=logger)
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 2: Download and extract UniProt Swiss-Prot FASTA file.
        # Skips download if the file already exists and --resolve_force_download is not set.
        download_uniprot_sprot(force_download=force_download, logger=logger)
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 3: Build BLAST database using makeblastdb.
        # Skips building if database already exists unless --resolve_force_build is set.
        make_blast_database(force_rebuild=force_build, logger=logger)
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 4: Run BLASTP using query FASTA and UniProt BLAST DB.
        # Output: BLAST results saved to .tsv file.
        blastp_alignment(output_path=blast_result_path, logger=logger)
        logger.add_spacer(level=logging.INFO, lines=1)

        # Step 5: Merge BLAST result with AMP metadata and resolve full NCBI taxonomy lineages.
        # Output: data/processed/dbAMP_with_blast_tax.csv
        integrate_amp_taxonomy_info(
            parsed_blast_path=blast_result_path,
            output_path=resolve_output_path,
            logger=logger,
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
