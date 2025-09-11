# pylint: disable=line-too-long, import-error, wrong-import-position, broad-exception-caught, too-many-statements
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
LOGGING_PATH = os.path.join(BASE_PATH, "src/utils/logging_toolkit/src/python")

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
        "import_path": "src.preprocess.dbAMP.resolve_taxid_by_blast.run_resolve_taxid_by_blast",
    },
    "taxonomy_nomatch_mapping": {
        "title": "Fill unmatched taxonomy via manual mapping + NCBI",
        "import_path": "src.preprocess.dbAMP.prepare_nomatch_mapping.run_prepare_nomatch_mapping",
    },
    "merge_checked_taxonomies": {
        "title": "Merge checked taxonomy tables and append BLAST-derived results",
        "import_path": "src.preprocess.dbAMP.merge_checked_taxonomies.run_merge_checked_taxonomies",
    },
    "clean_targets": {
        "title": "Clean AMP Targets text to long table",
        "import_path": "src.preprocess.dbAMP.clean_targets.run_clean_targets",
    },
    "targets_mapping": {
        "title": "Resolve AMP targets via manual mapping + BacDive",
        "import_path": "src.preprocess.dbAMP.resolve_targets_by_bacdive.run_prepare_targets_mapping",
    },
    "microbiome_composition": {
        "title": "Analyze 16S microbiome composition with nf-core/ampliseq",
        "import_path": "src.analysis.HMP2.microbiome_composition.run_microbiome_composition",
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
    logger.add_spacer(level=logging.INFO, lines=1)


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
            "AMPscope - Accelerating the Discovery of Antimicrobial Peptides through Computational Intelligence Pipeline\n\n"
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
        "--resolve_output_dir",
        type=str,
        help=(
            "Directory to save all output files for the resolve_taxid_by_blast stage. "
            "If not provided, a timestamped folder will be created inside data/interim/blast"
        ),
    )
    parser.add_argument(
        "--resource_dir",
        type=str,
        default=os.path.join(BASE_PATH, "data/external/uniprot"),
        help=(
            "Directory to store and manage all external resources "
            "(e.g., UniProt FASTA, BLAST databases, SILVA references). "
            "Default: data/external/"
        ),
    )
    parser.add_argument(
        "--resolve_blast_evalue",
        type=float,
        default=1e-5,
        help="E-value threshold for BLASTP. Default: 1e-5",
    )
    parser.add_argument(
        "--resolve_blast_max_hits",
        type=int,
        default=25,
        help="Max number of BLAST hits per query. Default: 25",
    )
    parser.add_argument(
        "--resolve_blast_max_hsps",
        type=int,
        default=1,
        help="Max HSPs per subject to report. Default: 1",
    )
    parser.add_argument(
        "--resolve_blast_pident_min",
        type=float,
        default=30.0,
        help="Minimum percent identity threshold for BLASTP filtering. Default: 30.0",
    )
    parser.add_argument(
        "--resolve_blast_qcovs_min",
        type=float,
        default=60.0,
        help="Minimum query coverage (qcovs) threshold for BLASTP filtering. Default: 60.0",
    )
    parser.add_argument(
        "--resolve_blast_threads",
        type=int,
        default=4,
        help="Number of CPU threads for BLASTP. Default: 4",
    )
    parser.add_argument(
        "--resolve_blast_pe_max",
        type=int,
        default=2,
        help="Maximum UniProt Protein Existence (PE) level to keep (1-5). Use 1-2 for high confidence; set to 3-5 to relax; set to -1 to disable.",
    )

    # -------------------- Taxonomy Nomatch Mapping Parameters --------------------
    parser.add_argument(
        "--nomatch_input_yaml",
        type=str,
        help="Path to YAML of manual taxonomy mappings (default: data/manual/taxid_mapping/notna_nomatch_mapping.yml)",
    )
    parser.add_argument(
        "--nomatch_output_csv",
        type=str,
        help="Path to save resolved mapping table CSV (default: data/processed/dbAMP/resolved_manual_taxonomy.csv)",
    )
    parser.add_argument(
        "--nomatch_amp_input_csv",
        type=str,
        help="Path to AMP table to apply mapping (default: data/interim/resolve_taxid/02_tax_check/notna_nomatch.csv)",
    )
    parser.add_argument(
        "--nomatch_amp_output_csv",
        type=str,
        help="Path to save AMP table with mapped taxonomy (default: data/interim/resolve_taxid/02_tax_check/notna_nomatch_filled.csv)",
    )

    # -------------------- Merge Checked Taxonomies Parameters --------------------
    parser.add_argument(
        "--merge_multi_path",
        type=str,
        help="Path to notna_multi.csv (default: data/interim/resolve_taxid/02_tax_check/notna_multi.csv)",
    )
    parser.add_argument(
        "--merge_single_path",
        type=str,
        help="Path to notna_single.csv (default: data/interim/resolve_taxid/02_tax_check/notna_single.csv)",
    )
    parser.add_argument(
        "--merge_filled_path",
        type=str,
        help="Path to notna_nomatch_filled.csv (default: data/interim/resolve_taxid/02_tax_check/notna_nomatch_filled.csv)",
    )
    parser.add_argument(
        "--merge_na_with_blast_path",
        type=str,
        help="Path to na_with_blast.csv (default: data/interim/blast/06_merge/na_with_blast.csv)",
    )
    parser.add_argument(
        "--merge_output_path",
        type=str,
        help="Path to save merged CSV (checked-only). Final '+_with_blast' is auto-saved next to it. (default: data/processed/dbAMP/checked_taxonomies.csv)",
    )

    # -------------------- Clean Targets Parameters --------------------
    parser.add_argument(
        "--targets_input_path",
        type=str,
        help="Path to the AMP table (Excel/CSV/TSV) with columns [dbAMP_ID, Targets] "
        "(default: data/raw/dbAMP/dbAMP3_pepinfo.xlsx)",
    )
    parser.add_argument(
        "--targets_output_csv",
        type=str,
        help="Path to save the cleaned long-format targets CSV "
        "(default: data/processed/dbAMP/targets_clean.csv)",
    )

    # -------------------- Prepare Targets Mapping Parameters --------------------
    parser.add_argument(
        "--targetsmap_input_yaml",
        type=str,
        help="Path to YAML of manual target mappings "
        "(default: data/manual/targets_mapping/targets_mapping.yml)",
    )
    parser.add_argument(
        "--targetsmap_bacdive_config",
        type=str,
        help="Path to BacDive JSON config with credentials "
        "(default: configs/bacdive.json)",
    )
    parser.add_argument(
        "--targetsmap_output_csv",
        type=str,
        help="Path to save resolved target mapping table "
        "(default: data/processed/dbAMP/resolved_manual_targets.csv)",
    )
    parser.add_argument(
        "--targetsmap_amp_input_csv",
        type=str,
        help="Path to AMP table to apply target mapping "
        "(default: data/interim/resolve_targets/targets_clean.csv)",
    )
    parser.add_argument(
        "--targetsmap_amp_output_csv",
        type=str,
        help="Path to save AMP table with resolved targets "
        "(default: data/processed/dbAMP/targets_resolved.csv)",
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
    parser.add_argument(
        "--microbiome_dada2_threads",
        type=int,
        default=4,
        help="Number of threads for DADA2 step in nf-core/ampliseq (default: 4)",
    )
    parser.add_argument(
        "--microbiome_fastp_threads",
        type=int,
        default=4,
        help="Number of threads for Fastp pre-processing (default: 4)",
    )
    parser.add_argument(
        "--microbiome_cutadapt_threads",
        type=int,
        default=4,
        help="Number of threads for Cutadapt trimming (default: 4)",
    )
    parser.add_argument(
        "--microbiome_max_cpus",
        type=int,
        default=8,
        help="Maximum number of CPUs Nextflow can use for the entire workflow (default: 8)",
    )
    parser.add_argument(
        "--microbiome_max_memory",
        type=str,
        default="16.GB",
        help="Maximum memory available to the pipeline (e.g., '16.GB') (default: 16.GB)",
    )
    parser.add_argument(
        "--microbiome_max_time",
        type=str,
        default="24.h",
        help="Maximum execution time for the pipeline (e.g., '24.h') (default: 24.h)",
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
