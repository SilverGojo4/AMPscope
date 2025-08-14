# pylint: disable=import-error, wrong-import-position, too-many-function-args, broad-exception-caught
"""
AMPscope - Project Main Entry Point

This script serves as the centralized execution interface for the AMPscope pipeline.
"""
# ============================== Standard Library Imports ==============================
import argparse
import importlib
import logging
import os
import sys

# ============================== Project Root Path Setup ==============================
BASE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
LOGGING_PATH = os.path.join(BASE_PATH, "src/logging_toolkit/src/python")

# Fix Python Path
if BASE_PATH not in sys.path:
    sys.path.append(BASE_PATH)
if LOGGING_PATH not in sys.path:
    sys.path.append(LOGGING_PATH)

# ============================== Project-Specific Imports ==============================
# Logging configuration and custom logger
from setup_logging import setup_logging

# ============================== Stage Configuration ==============================
SUPPORTED_STAGES = {
    "resolve_taxid_by_blast": {
        "title": "Resolve missing AMP TaxIDs via BLAST alignment",
        "import_path": "src.data.dbAMP.resolve_taxid_by_blast.run_resolve_taxid_by_blast",
    },
    "microbiome_composition": {
        "title": "Analyze 16S microbiome composition with nf-core/ampliseq",
        "import_path": "src.data.HMP2.microbiome_composition.run_microbiome_composition",
    },
}


# ============================== Pipeline Dispatcher ==============================
def dispatch_stage(args: argparse.Namespace) -> None:
    """
    Dispatch execution to the appropriate pipeline stage using lazy import.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing stage-specific options.
    """
    # -------------------- Stage Validation --------------------
    # Ensure that the provided stage name is supported.
    stage = args.stage.lower()
    if stage not in SUPPORTED_STAGES:
        available = ", ".join(SUPPORTED_STAGES.keys())
        raise ValueError(f"Unknown stage '{stage}'. Available stages: {available}.")

    # -------------------- Log Path Validation --------------------
    # Ensure that the --log_path argument is provided.
    if not args.log_path:
        raise ValueError("Please specify '--log_path' to enable logging output.")

    # -------------------- Dynamic Stage Import --------------------
    # Perform a lazy import of the selected pipeline stage based on its import path.
    # This avoids loading all modules at startup and improves modularity.
    stage_info = SUPPORTED_STAGES[stage]
    module_path, func_name = stage_info["import_path"].rsplit(".", 1)
    module = importlib.import_module(module_path)
    stage_func = getattr(module, func_name)

    # -------------------- Logger Initialization --------------------
    # Set up the logging environment, loading the logging config and preparing output directory.
    log_config_file = os.path.join(BASE_PATH, "configs/logging.json")
    stage_log_path = os.path.abspath(args.log_path)
    stage_log_dir = os.path.dirname(stage_log_path)
    os.makedirs(name=stage_log_dir, exist_ok=True)

    logger = setup_logging(
        input_config_file=log_config_file,
        logger_name="combo_logger",
        handler_name="file",
        output_log_path=stage_log_path,
    )

    # -------------------- Log Pipeline Metadata --------------------
    logger.info("[ 'Pipeline Initialization Summary' ]")
    logger.log_pipeline_initialization(
        project_name="AMPscope",
        line_width=120,
    )
    logger.add_spacer(level=logging.INFO, lines=1)

    # -------------------- Execute Stage Function --------------------
    # Dynamically call the selected stage function, passing all stage-specific arguments.
    extra_args = vars(args)
    extra_args.pop("stage", None)
    extra_args.pop("log_path", None)
    stage_func(base_path=BASE_PATH, logger=logger, **extra_args)


# ============================== Main Entry ==============================
def main():
    """
    Main CLI entry point for the AMPscope pipeline.
    Parses CLI arguments and routes execution to the selected pipeline stages.
    """
    # -------------------- Argument Parser --------------------
    available_stage_lines = [
        f"  - {stage:<15} {info['title']}" for stage, info in SUPPORTED_STAGES.items()
    ]
    available_stages_text = "\n".join(available_stage_lines)
    example_stage = list(SUPPORTED_STAGES.keys())[0]
    example_command = (
        f"  python main.py --stage {example_stage} --log_path logs/{example_stage}.log"
    )

    parser = argparse.ArgumentParser(
        description=(
            "Min - Phosphoproteome-derived Pan-cancer Kinase Activity Profiling Pipeline\n\n"
            "Available stages:\n"
            f"{available_stages_text}\n\n"
            "Example:\n"
            f"{example_command}\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # -------------------- General Options --------------------
    parser.add_argument(
        "--stage",
        type=str,
        choices=list(SUPPORTED_STAGES.keys()),
        required=True,
        help="Pipeline stage to run. Choose one of: "
        + ", ".join(SUPPORTED_STAGES.keys()),
    )
    parser.add_argument(
        "--log_path",
        type=str,
        required=True,
        help="Path for log file output.",
    )

    # -------------------- Resolve TaxID via BLAST Parameters --------------------
    parser.add_argument(
        "--resolve_input_dbamp",
        type=str,
        help="Path to the input dbAMP file (.xlsx) for resolving TaxID",
    )
    parser.add_argument(
        "--resolve_force_download",
        action="store_true",
        help="Force re-download of UniProt FASTA file even if it exists",
    )
    parser.add_argument(
        "--resolve_force_build",
        action="store_true",
        help="Force rebuild of the UniProt BLAST database even if it already exists",
    )
    parser.add_argument(
        "--resolve_blast_result_path",
        type=str,
        help="Path to BLAST result file (.tsv) used for resolving TaxID",
    )
    parser.add_argument(
        "--resolve_output_path",
        type=str,
        help="Path to save the TaxID-resolved dbAMP output (.csv)",
    )

    # -------------------- Microbiome Composition Parameters --------------------
    parser.add_argument(
        "--microbiome_input_csv",
        type=str,
        help="Path to the input sample sheet CSV for nf-core/ampliseq",
    )
    parser.add_argument(
        "--microbiome_metadata",
        type=str,
        help="Path to the sample metadata TSV file",
    )
    parser.add_argument(
        "--microbiome_silva_train",
        type=str,
        help="Path to custom SILVA training set FASTA file",
    )
    parser.add_argument(
        "--microbiome_silva_species",
        type=str,
        help="Path to custom SILVA species FASTA file",
    )
    parser.add_argument(
        "--microbiome_output_dir",
        type=str,
        help="Output directory for nf-core/ampliseq results",
    )

    args = parser.parse_args()

    # -------------------- Stage Execution --------------------
    try:
        dispatch_stage(args)
        print("Pipeline execution completed successfully.")
        sys.exit(0)

    except Exception as e:
        print(f"[Pipeline Error] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
