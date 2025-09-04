# pylint: disable=line-too-long, import-error, wrong-import-position, too-many-arguments, too-many-positional-arguments
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
from utils.io_utils import directory_exists, file_exists
from utils.log_utils import get_pipeline_completion_message, get_task_completion_message


# ============================== Custom Functions ==============================
def nfcore_ampliseq(
    input_csv: str,
    metadata_tsv: str,
    silva_train: str,
    silva_species: str,
    output_dir: str,
    logger: CustomLogger,
    dada2_threads: int = 4,
    fastp_threads: int = 4,
    cutadapt_threads: int = 4,
    max_cpus: int = 16,
    max_memory: str = "64.GB",
    max_time: str = "48.h",
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
    dada2_threads : int
        Threads for DADA2 denoising step.
    fastp_threads : int
        Threads for fastp preprocessing step.
    cutadapt_threads : int
        Threads for cutadapt trimming step.
    max_cpus : int
        Max CPUs to allocate for any task.
    max_memory : str
        Max memory to allocate (e.g., '64.GB').
    max_time : str
        Max time to allocate (e.g., '48.h').

    Returns
    -------
    None
    """
    logger.info("/ Task: Run 'nf-core/ampliseq' via Nextflow for microbiome profiling")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Validate files
        for path in [input_csv, metadata_tsv, silva_train, silva_species]:
            if not file_exists(file_path=path):
                raise FileNotFoundError(f"File not found: '{path}'")

        # Log parameters
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "Running nf-core/ampliseq with the following parameters:\n"
                f"- 'Input samplesheet' : '{input_csv}'\n"
                f"- 'Sample metadata' : '{metadata_tsv}'\n"
                f"- 'SILVA train FASTA' : '{silva_train}'\n"
                f"- 'SILVA species FASTA' : '{silva_species}'\n"
                f"- 'Output directory' : '{output_dir}'\n"
                f"- 'Profile' : 'docker'\n"
                f"- 'Nextflow revision' : 2.11.0\n"
                f"- 'Threads (dada2/fastp/cutadapt)' : {dada2_threads}/{fastp_threads}/{cutadapt_threads}\n"
                f"- 'Resources (cpus/mem/time)' : {max_cpus}/{max_memory}/{max_time}\n"
                f"- 'Flags' : '--single_end, --skip_cutadapt'"
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
            "--dada2_threads",
            str(dada2_threads),
            "--fastp_threads",
            str(fastp_threads),
            "--cutadapt_threads",
            str(cutadapt_threads),
            "--max_cpus",
            str(max_cpus),
            "--max_memory",
            str(max_memory),
            "--max_time",
            str(max_time),
        ]
        subprocess.run(cmd, check=True)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"Results saved to:\n'{output_dir}'",
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
        input_csv = kwargs.get(
            "microbiome_input_csv",
            os.path.join(base_path, "data/raw/HMP2/biopsy_16S/samplesheet.csv"),
        )
        metadata_tsv = kwargs.get(
            "microbiome_metadata",
            os.path.join(base_path, "data/interim/nfcore_metadata.tsv"),
        )
        silva_train = kwargs.get(
            "microbiome_silva_train",
            os.path.join(
                base_path, "data/external/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
            ),
        )
        silva_species = kwargs.get(
            "microbiome_silva_species",
            os.path.join(
                base_path, "data/external/silva_species_assignment_v138.1.fa.gz"
            ),
        )
        output_dir = kwargs.get(
            "microbiome_output_dir", os.path.join(base_path, "experiments/HMP2")
        )
        dada2_threads = int(kwargs.get("microbiome_dada2_threads", 4))
        fastp_threads = int(kwargs.get("microbiome_fastp_threads", 4))
        cutadapt_threads = int(kwargs.get("microbiome_cutadapt_threads", 4))
        max_cpus = int(kwargs.get("microbiome_max_cpus", 16))
        max_memory = kwargs.get("microbiome_max_memory", "64.GB")
        max_time = kwargs.get("microbiome_max_time", "48.h")

        # Check folders exist
        if not directory_exists(dir_path=output_dir):
            os.makedirs(output_dir)

        # -------------------- Pipeline Execution --------------------
        # Step 1: Run nf-core/ampliseq pipeline
        nfcore_ampliseq(
            input_csv=input_csv,
            metadata_tsv=metadata_tsv,
            silva_train=silva_train,
            silva_species=silva_species,
            output_dir=output_dir,
            logger=logger,
            dada2_threads=dada2_threads,
            fastp_threads=fastp_threads,
            cutadapt_threads=cutadapt_threads,
            max_cpus=max_cpus,
            max_memory=max_memory,
            max_time=max_time,
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
