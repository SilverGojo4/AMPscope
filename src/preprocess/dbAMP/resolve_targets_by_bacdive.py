# pylint: disable=line-too-long, import-error, wrong-import-position, broad-exception-caught, too-many-locals, too-many-statements, too-many-arguments, too-many-positional-arguments
"""
Resolve AMP targets via manual mapping + BacDive lookup
"""
# ============================== Standard Library Imports ==============================
import logging
import os
import sys
import time

# ============================== Third-Party Library Imports ==============================
import pandas as pd
import yaml
from bacdive import BacdiveClient

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
    load_json_config,
)
from utils.log_utils import get_pipeline_completion_message, get_task_completion_message


# ============================== Custom Functions ==============================
def extract_gram_stains(morph: dict) -> dict[str, int]:
    """
    Extract and count Gram stain information from a BacDive strain's morphology record.

    Parameters
    ----------
    morph : dict
        The 'Morphology' section of a BacDive strain record.

    Returns
    -------
    dict of str -> int
        Dictionary with counts of Gram stain results:
        {
            "positive": <count>,
            "negative": <count>
        }
        Returns {"positive": 0, "negative": 0} if no Gram stain information is present.
    """
    counts = {"positive": 0, "negative": 0}
    if not isinstance(morph, dict):
        return counts

    cell_morph = morph.get("cell morphology")

    # Case 1: single dict
    if isinstance(cell_morph, dict):
        gram = cell_morph.get("gram stain")
        if gram:
            gram_norm = str(gram).strip().lower()
            if gram_norm in counts:
                counts[gram_norm] += 1

    # Case 2: list of dicts
    elif isinstance(cell_morph, list):
        for item in cell_morph:
            if isinstance(item, dict):
                gram = item.get("gram stain")
                if gram:
                    gram_norm = str(gram).strip().lower()
                    if gram_norm in counts:
                        counts[gram_norm] += 1

    return counts


def cached_gram_lookup(
    name: str,
    client: BacdiveClient,
    gram_cache: dict[str, dict],
    logger: CustomLogger,
) -> dict:
    """
    Query BacDive for Gram stain information of a given target name,
    with caching to avoid redundant queries.

    Parameters
    ----------
    name : str
        Target name (normalized, after manual mapping).
    client : BacdiveClient
        An authenticated BacDive API client.
    gram_cache : dict
        Cache dictionary {name -> result dict}.
    logger : CustomLogger
        Logger instance for tracking progress.

    Returns
    -------
    dict
        Result dictionary with Gram stain information:
        {
            "positive": int | None,
            "negative": int | None,
            "final": str
        }

        Possible values of "final":
        - "positive"     → majority of Gram stain results are positive
        - "negative"     → majority of Gram stain results are negative
        - "unclassified" → strains found but no Gram info, or tie (0/0, equal counts)
        - "not_found"    → no strains found in BacDive
        - "error"        → API or query error occurred
    """
    # Cache hit
    if name in gram_cache:
        return gram_cache[name]

    try:
        # Query BacDive
        count = client.search(taxonomy=name)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"[Querying] '{name}' → {count} strains found",
            border="|",
            length=120,
        )

        # Case 1: No strains found
        if count == 0:
            result = {"positive": None, "negative": None, "final": "not_found"}
            gram_cache[name] = result
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
            return result

        positive_count = 0
        negative_count = 0

        # Case 2: Iterate over retrieved strains
        for strain in client.retrieve():
            morph = strain.get("Morphology", {})
            if isinstance(morph, dict):
                gram_counts = extract_gram_stains(morph)
                positive_count += gram_counts["positive"]
                negative_count += gram_counts["negative"]

        # Case 3: Decide final classification
        if positive_count > negative_count:
            final = "positive"
        elif negative_count > positive_count:
            final = "negative"
        else:
            final = "unclassified"

        result = {
            "positive": positive_count,
            "negative": negative_count,
            "final": final,
        }
        gram_cache[name] = result

        logger.log_with_borders(
            level=logging.INFO,
            message=(
                f"▸ Final result: '{final}'\n"
                f"  [+] Positive: {positive_count}\n"
                f"  [-] Negative: {negative_count}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")
        return result

    except Exception:
        # Case 4: Error during query
        logger.error(f"Error occurred while querying '{name}'")
        result = {"positive": None, "negative": None, "final": "error"}
        gram_cache[name] = result
        return result


def resolve_manual_targets(
    input_yaml: str,
    bacdive_config: str,
    output_csv: str,
    amp_input_csv: str,
    amp_output_csv: str,
    logger: CustomLogger,
) -> None:
    """
    Resolve AMP target using manually mapped names and BacDive.
    Then apply this mapping back to a given AMP table and export updated results.

    Parameters
    ----------
    input_yaml : str
        Path to the YAML file with manual mappings.
    bacdive_config : str
        Path to the JSON config file containing BacDive credentials (email, password).
    output_csv : str
        Path to save resolved target results as CSV (mapping table).
    amp_input_csv : str
        Path to the original long-format AMP targets table [dbAMP_ID, Targets].
    amp_output_csv : str
        Path to save the updated AMP targets table with resolved targets.
    logger : CustomLogger
        Logger instance to track progress.

    Returns
    -------
    None
    """
    logger.info("/ Task: Resolve targets via manual mapping + BacDive")
    logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

    # Start timing
    start_time = time.time()

    try:
        # Load manual mappings from YAML into DataFrame
        with open(input_yaml, "r", encoding="utf-8") as f:
            manual_mapping = yaml.safe_load(f)
        df_map = pd.DataFrame(
            manual_mapping.items(), columns=["original_name", "mapped_name"]
        )

        # Prepare query stats (before BacDive)
        skip_keywords = {"unknown", "unclassified", "na", "not available", "n/a", ""}

        def classify_target(name: str) -> str:
            """
            Classify target into queryable or skip categories.

            Parameters
            ----------
            name : str
                Target name after manual mapping and normalization.

            Returns
            -------
            str
                "queryable" if eligible for BacDive,
                otherwise one of {"placeholder", "homo_sapiens", "virus"}.
            """
            name_clean = str(name).strip()
            name_lower = name_clean.lower()

            # 1. Placeholders
            if name_lower in skip_keywords:
                return "placeholder"

            # 2. Human (cell lines or disease-specific)
            if name_clean.startswith("Homo sapiens ("):
                return "homo_sapiens"

            # 3. Virus
            virus_keywords = [
                "virus",
                "hcov",
                "influenza",
                "hepatitis",
                "herpes",
                "hiv",
                "sars",
                "coronavirus",
            ]
            if any(vk in name_lower for vk in virus_keywords):
                return "virus"

            return "queryable"

        # Apply classification
        df_map["skip_reason"] = df_map["mapped_name"].apply(classify_target)
        df_map["is_queryable"] = df_map["skip_reason"] == "queryable"
        total_rows = df_map.shape[0]
        num_to_query = int(df_map["is_queryable"].sum())
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Manual mapping query plan ]\n"
                f"▸ total rows             : {total_rows}\n"
                f"  - to query (BacDive)   : {num_to_query}\n"
                f"  - skipped total        : {total_rows - num_to_query}\n"
                f"    • placeholders       : {int((df_map['skip_reason'] == 'placeholder').sum())}\n"
                f"    • Homo sapiens       : {int((df_map['skip_reason'] == 'homo_sapiens').sum())}\n"
                f"    • virus              : {int((df_map['skip_reason'] == 'virus').sum())}"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # BacDive login
        try:
            config = load_json_config(bacdive_config)
            email, password = config.get("email"), config.get("password")
            if not email or not password:
                raise ValueError("BacDive credentials missing in config file")
            client = BacdiveClient(email, password)
            client.setSearchType("exact")
            logger.log_with_borders(
                level=logging.INFO,
                message="Successfully logged into BacDive API ('Authentication successful')",
                border="|",
                length=120,
            )
            logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        except Exception as exc:
            raise SystemExit(
                "Failed to log into BacDive API (Authentication failed)"
            ) from exc

        # Query BacDive for Gram info
        gram_cache: dict[str, dict] = {}
        results = []
        for _, row in df_map.iterrows():
            name = row["mapped_name"]
            if not row["is_queryable"]:
                results.append(
                    {"positive": None, "negative": None, "final": row["skip_reason"]}
                )
                continue

            result = cached_gram_lookup(name, client, gram_cache, logger)
            results.append(result)

        # Expand results into DataFrame
        df_results = pd.DataFrame(results)
        df_map = pd.concat([df_map, df_results], axis=1)
        df_map["positive"] = pd.to_numeric(df_map["positive"], errors="coerce").astype(
            "Int64"
        )
        df_map["negative"] = pd.to_numeric(df_map["negative"], errors="coerce").astype(
            "Int64"
        )
        summary_counts = df_map["final"].value_counts(dropna=False).to_dict()
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ BacDive query results summary ]\n"
                + "\n".join(f"▸ {k:12s}: {v}" for k, v in summary_counts.items())
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Apply mapping to AMP table
        required_columns = ["dbAMP_ID", "Targets"]
        df_targets = load_dataframe_by_columns(
            amp_input_csv, required_columns=required_columns, has_header=True
        )
        df_merged = df_targets.merge(
            df_map, how="left", left_on="Targets", right_on="original_name"
        )
        df_merged["Targets"] = df_merged["mapped_name"].fillna(df_merged["Targets"])
        df_merged.drop(
            columns=["original_name", "mapped_name"], inplace=True, errors="ignore"
        )

        # Target-level classification stats
        target_class_counts = df_merged["final"].value_counts(dropna=False).to_dict()
        total_targets = len(df_merged)
        max_class_len = max(
            (len(str(cls)) for cls in target_class_counts.keys()), default=10
        )
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Target-level Gram stain classification (after merge) ]\n"
                f"▸ Total targets  : {total_targets}\n"
                + "\n".join(
                    f"  - {str(cls):<{max_class_len}s} : {count} ('{count/total_targets:.2%}')"
                    for cls, count in target_class_counts.items()
                )
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Genus-level statistics
        genus_series = (
            df_merged.loc[
                df_merged["final"].isin(["positive", "negative", "unclassified"]),
                "Targets",
            ]
            .dropna()
            .apply(lambda x: str(x).split()[0])
        )
        genus_counts = genus_series.value_counts().to_dict()
        unique_genera = len(genus_counts)
        max_genus_len = max(
            (len(str(genus)) for genus in genus_counts.keys()), default=10
        )
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ Genus-level statistics (resolved only) ]\n"
                f"▸ Unique genera : {unique_genera}\n"
                + "\n".join(
                    f"  - {genus:<{max_genus_len}s} : {count}"
                    for genus, count in genus_counts.items()
                )
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # AMP–Target linkage statistics
        amp_target_counts = df_merged.groupby("dbAMP_ID")["Targets"].nunique()
        total_amps = amp_target_counts.shape[0]
        total_unique_targets = df_merged["Targets"].nunique()
        avg_targets_per_amp = amp_target_counts.mean()
        max_targets = amp_target_counts.max()
        min_targets = amp_target_counts.min()
        amp_with_max = amp_target_counts.idxmax()
        amp_with_min = amp_target_counts.idxmin()
        logger.log_with_borders(
            level=logging.INFO,
            message=(
                "[ AMP-Target linkage statistics ]\n"
                f"▸ total AMP entries       : {total_amps}\n"
                f"▸ total unique targets    : {total_unique_targets}\n"
                f"▸ average targets per AMP : {avg_targets_per_amp:.2f}\n"
                f"▸ max targets per AMP     : '{amp_with_max}' ({max_targets} targets)\n"
                f"▸ min targets per AMP     : '{amp_with_min}' ({min_targets} target)"
            ),
            border="|",
            length=120,
        )
        logger.add_divider(level=logging.INFO, length=120, border="+", fill="-")

        # Save result to CSV
        df_map.to_csv(output_csv, index=False)
        df_merged.to_csv(amp_output_csv, index=False)
        logger.log_with_borders(
            level=logging.INFO,
            message=f"CSV saved:\n'{output_csv}'\n{amp_output_csv}",
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
        logger.exception("Unexpected error in 'resolve_manual_targets()'")
        raise


# ============================== Pipeline Entry Point ==============================
def run_prepare_targets_mapping(base_path: str, logger: CustomLogger, **kwargs) -> None:
    """
    Resolve AMP targets using manual mappings (YAML)
    combined with BacDive lookup. Saves both the resolved
    mapping table and the updated AMP table.

    Parameters
    ----------
    base_path : str
        Project root path.
    logger : CustomLogger
        Logger instance for tracking progress.
    **kwargs : dict
        Optional overrides for input/output file paths.

    Returns
    -------
    None
    """
    # Start timing
    start_time = time.time()

    try:
        # -------------------- Retrieve input parameters (CLI or default) --------------------
        input_yaml = kwargs.get(
            "targetsmap_input_yaml",
            os.path.join(base_path, "data/manual/targets_mapping/targets_mapping.yml"),
        )
        bacdive_config = kwargs.get(
            "targetsmap_bacdive_config",
            os.path.join(base_path, "configs/bacdive.json"),
        )
        output_csv = kwargs.get(
            "targetsmap_output_csv",
            os.path.join(base_path, "data/processed/dbAMP/resolved_manual_targets.csv"),
        )
        amp_input_csv = kwargs.get(
            "targetsmap_amp_input_csv",
            os.path.join(base_path, "data/interim/resolve_targets/targets_clean.csv"),
        )
        amp_output_csv = kwargs.get(
            "targetsmap_amp_output_csv",
            os.path.join(base_path, "data/processed/dbAMP/targets_resolved.csv"),
        )

        # Check input files exist
        for path in [input_yaml, bacdive_config, amp_input_csv]:
            if not file_exists(file_path=path):
                raise FileNotFoundError(f"File not found: '{path}'")

        # Ensure output directories exist
        for path in [output_csv, amp_output_csv]:
            output_dir = os.path.dirname(path)
            if not directory_exists(dir_path=output_dir):
                os.makedirs(name=output_dir)

        # -------------------- Pipeline Execution --------------------
        resolve_manual_targets(
            input_yaml=input_yaml,
            bacdive_config=bacdive_config,
            output_csv=output_csv,
            amp_input_csv=amp_input_csv,
            amp_output_csv=amp_output_csv,
            logger=logger,
        )
        logger.add_spacer(level=logging.INFO, lines=1)

    except Exception:
        logger.exception("Critical failure in 'run_prepare_targets_mapping()'")
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
