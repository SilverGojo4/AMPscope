# pylint: disable=import-error, wrong-import-position, broad-exception-caught, too-many-arguments, too-many-positional-arguments
"""
Run nf-core/ampliseq pipeline for microbiome 16S analysis
"""

# ============================== Standard Library Imports ==============================
import logging
import os
import subprocess
import sys
import time

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
from utils.log_utils import get_pipeline_completion_message, get_task_completion_message


# ============================== Custom Functions ==============================
def nfcore_ampliseq(
    input_csv: str,
    metadata_tsv: str,
    silva_train: str,
    silva_species: str,
    output_dir: str,
    logger: CustomLogger,
) -> None:
    """
    Execute nf-core/ampliseq pipeline using Nextflow for microbiome composition analysis.

    Parameters
    ----------
    input_csv : str
        Path to the sample sheet CSV input.
    metadata_tsv : str
        Path to the sample metadata TSV.
    silva_train : str
        Path to custom SILVA training set FASTA.
    silva_species : str
        Path to custom SILVA species FASTA.
    output_dir : str
        Output directory for results.
    logger : CustomLogger
        Logger instance for tracking.

    Returns
    -------
    None
    """
    logger.info("/ Task: Run 'nf-core/ampliseq' via Nextflow for microbiome profiling")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Log parameters
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "Running nf-core/ampliseq with the following parameters:\n"
                f"- 'Input samplesheet'     : '{input_csv}'\n"
                f"- 'Sample metadata'       : '{metadata_tsv}'\n"
                f"- 'SILVA train FASTA'     : '{silva_train}'\n"
                f"- 'SILVA species FASTA'   : '{silva_species}'\n"
                f"- 'Output directory'      : '{output_dir}'\n"
                f"- 'Profile'               : 'docker'\n"
                f"- 'Nextflow revision'     : 2.11.0\n"
                f"- 'Flags'                 : '--single_end, --skip_cutadapt'"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Run nf-core/ampliseq
        cmd = [
            "nextflow",
            "run",
            "nf-core/ampliseq",
            "-r",
            "2.11.0",
            "-profile",
            "docker",
            "--input",
            input_csv,
            "--single_end",
            "--skip_cutadapt",
            "--metadata",
            metadata_tsv,
            "--dada_ref_tax_custom",
            silva_train,
            "--dada_ref_tax_custom_sp",
            silva_species,
            "--outdir",
            output_dir,
        ]
        subprocess.run(cmd, check=True)

        logger.log_with_borders(
            level=logging.INFO,
            message="'nf-core/ampliseq' completed successfully.",
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Results saved to:\n'{output_dir}'",
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
        logger.exception("Unexpected error during 'nfcore_ampliseq()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_microbiome_composition(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Run microbiome composition analysis pipeline using nf-core/ampliseq.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance for progress tracking.
    kwargs : dict
        Optional keyword arguments (for CLI overrides).

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        # -------------------- Retrieve input parameters (CLI or defaults) --------------------
        input_csv = kwargs.get("microbiome_input_csv") or os.path.join(
            base_path, "data/raw/HMP2/biopsy_16S/samplesheet.csv"
        )
        metadata_tsv = kwargs.get("microbiome_metadata") or os.path.join(
            base_path, "data/interim/nfcore_metadata.tsv"
        )
        silva_train = kwargs.get("microbiome_silva_train") or os.path.join(
            base_path, "data/external/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
        )
        silva_species = kwargs.get("microbiome_silva_species") or os.path.join(
            base_path, "data/external/silva_species_assignment_v138.1.fa.gz"
        )
        output_dir = kwargs.get("microbiome_output_dir") or os.path.join(
            base_path, "experiments/HMP2"
        )

        # -------------------- Pipeline Execution --------------------
        # Step 1: Run nf-core/ampliseq pipeline
        nfcore_ampliseq(
            input_csv=input_csv,
            metadata_tsv=metadata_tsv,
            silva_train=silva_train,
            silva_species=silva_species,
            output_dir=output_dir,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_microbiome_composition()'")
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
