"""
Pre-flight validation of a SQANTI-evidence configuration file.

Called by the `sqanti_evidence` wrapper before Snakemake starts so that mistakes are reported
with a clear message and a non-zero exit code instead of a Snakemake traceback. The same rules
are enforced (with defaults filled in) by `rules/setup/functions.smk:validate_and_fill_config`
when Snakemake parses the workflow; this module only validates what the user wrote.

Relative paths are interpreted relative to the current working directory, exactly as Snakemake
does when it builds the DAG.
"""
import os
import sys
import logging

# Use the pipeline logger if available, otherwise create a basic one
logger = logging.getLogger('pipeline')
if not logger.handlers:
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)

ALLOWED_VALUES = {
    ("training", "mode"): ["mixed", "busco_only", "sqanti_only"],
    ("prediction", "mode"): ["split", "full"],
    ("prediction", "filter_mode"): ["strict", "medium", "none"],
    ("curation", "mode"): ["placebo", "user", "ab_initio"],
    ("curation", "data_type"): ["pacbio", "pacbio_ccs", "nanopore", "ont", "assembly", "transcripts"],
}

# Optional file paths: validated only when the user set a non-empty value.
OPTIONAL_FILES = [
    ("project", "prediction_genome"),
    ("prediction", "hint_config"),
    ("prediction", "hint_weights"),
    ("curation", "filter_rules"),
    ("evaluation", "reference_gtf"),
]

# Options run_isoquant sets itself; passing them again through curation.isoquant_args is rejected.
ISOQUANT_RESERVED_OPTIONS = ["--reference", "-r", "--data_type", "-d", "--prefix", "-p", "--threads", "-t",
                             "-o", "--output", "--fastq", "--bam", "--unmapped_bam", "--fastq_list", "--bam_list"]


def validate_isoquant_args(extra):
    """Returns an error message if curation.isoquant_args repeats an option the rule already sets, else None."""
    if not extra:
        return None
    if not isinstance(extra, str):
        return f"ERROR: 'curation.isoquant_args' must be a string (current: {extra!r})"
    tokens = [t.split("=")[0] for t in extra.split()]
    clashes = [t for t in tokens if t in ISOQUANT_RESERVED_OPTIONS]
    if clashes:
        return (f"ERROR: 'curation.isoquant_args' must not contain {clashes}: the pipeline already sets "
                "the reference, input, data type, prefix, threads and output directory for IsoQuant.")
    return None


AB_INITIO_TRAINING_ERROR = (
    "ERROR: 'curation.mode' is 'ab_initio' but 'training.mode' is not 'busco_only'. "
    "To generate an ab initio annotation for SQANTI3 to classify against, the Augustus gene model "
    "must be trained using only BUSCO genes: the SQANTI3-derived training genes ('mixed' or "
    "'sqanti_only') do not exist yet at that point, which makes the workflow circular. "
    "Set training.mode: busco_only, or use curation.mode: placebo / user."
)


def validate_modes(training_mode, curation_mode):
    """
    Single source of truth for constraints between modes. Returns an error message or None.
    Used by the pre-flight check (wrapper) and by rules/setup/functions.smk (Snakemake parse time).
    """
    if curation_mode == "ab_initio" and training_mode != "busco_only":
        return AB_INITIO_TRAINING_ERROR
    return None


def die(msg):
    logger.error(msg)
    sys.exit(1)


def warn(msg):
    logger.warning(msg)


def validate_path(path, name, is_file=True):
    """Checks if a path exists and is the correct type."""
    if not isinstance(path, str) or not path:
        die(f"ERROR: '{name}' must be a non-empty path string (current: {path!r})")
    if not os.path.exists(path):
        die(f"ERROR: Missing required {name} at: {path}")
    if is_file and not os.path.isfile(path):
        die(f"ERROR: Expected {name} to be a file, but it is not: {path}")
    if not is_file and not os.path.isdir(path):
        die(f"ERROR: Expected {name} to be a directory, but it is not: {path}")
    logger.info(f"Validated {name}: {path}")


def validate_enum(config, section, key):
    """Dies if config[section][key] is set to a value outside ALLOWED_VALUES."""
    value = config.get(section, {}).get(key)
    if value is None:
        return  # functions.smk fills the default
    allowed = ALLOWED_VALUES[(section, key)]
    if value not in allowed:
        die(f"ERROR: '{section}.{key}' must be one of {allowed} (current: {value!r})")


def check_augustus_config_path(toolsdir):
    """
    Augustus writes the trained species model into $AUGUSTUS_CONFIG_PATH/species/<name>.
    The bioconda Augustus package sets that variable when its conda environment is activated,
    pointing inside <toolsdir>/conda_envs/<env>/config/, so what must be writable is the tools
    directory that will host the environment. If the user exported AUGUSTUS_CONFIG_PATH
    themselves it is checked too, with a warning that conda activation overrides it.
    """
    env_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
    if env_path:
        species_dir = os.path.join(env_path, "species")
        target = species_dir if os.path.isdir(species_dir) else env_path
        if not os.access(target, os.W_OK):
            die(f"ERROR: AUGUSTUS_CONFIG_PATH is set to '{env_path}' but it is not writable; "
                "Augustus training (new_species.pl / etraining) needs to create the species directory there.")
        warn(f"AUGUSTUS_CONFIG_PATH is set in your shell ({env_path}). Note that activating the Augustus "
             "conda environment overrides it with the environment's own config directory.")

    conda_envs = os.path.join(toolsdir, "conda_envs")
    target = conda_envs if os.path.isdir(conda_envs) else toolsdir
    if not os.access(target, os.W_OK):
        die(f"ERROR: '{target}' is not writable. Conda environments (including the Augustus config "
            "directory that receives the trained species model) are installed under project.toolsdir.")
    logger.info(f"Validated writable tools directory for Augustus species models: {target}")


def check_inputs(config):
    """
    Performs hard validation of the configuration before launching Snakemake.
    Focuses on file existence and parameter sanity.
    """
    logger.info("Starting comprehensive input validation...")

    # 1. Check for main sections
    required_sections = ["project", "training", "prediction", "curation"]
    for section in required_sections:
        if section not in config or not isinstance(config[section], dict):
            die(f"ERROR: Missing required configuration section: '{section}'")

    # 2. Project Section (Hard requirements)
    prj = config["project"]
    for key in ["genome", "input", "toolsdir"]:
        if key not in prj or not prj[key]:
            die(f"ERROR: Missing required variable 'project.{key}'")

    validate_path(prj["genome"], "project.genome", is_file=True)
    validate_path(prj["input"], "project.input", is_file=True)
    validate_path(prj["toolsdir"], "project.toolsdir", is_file=False)

    if not config["prediction"].get("species"):
        die("ERROR: Missing required variable 'prediction.species' (Augustus species name).")

    # 3. Output directory check (Warning only)
    outdir = prj.get("outdir", "SQANTI_evidence_results")
    if os.path.exists(outdir):
        warn(f"Output directory '{outdir}' already exists. Files may be overwritten.")

    # 4. Enumerated parameters
    for section, key in ALLOWED_VALUES:
        validate_enum(config, section, key)

    # 5. Mode combinations
    trn = config["training"]
    cur = config["curation"]
    training_mode = trn.get("mode", "mixed")
    curation_mode = cur.get("mode", "placebo")

    mode_error = validate_modes(training_mode, curation_mode)
    if mode_error:
        die(mode_error)

    iq_error = validate_isoquant_args(cur.get("isoquant_args", ""))
    if iq_error:
        die(iq_error)

    if training_mode in ["mixed", "busco_only"] and not trn.get("lineage"):
        die(f"ERROR: 'training.lineage' (BUSCO lineage) is required when training.mode is '{training_mode}'.")

    # 6. Numerical Parameter Validation
    if "miniprot_threshold" in trn:
        try:
            val = float(trn["miniprot_threshold"])
        except (ValueError, TypeError):
            die(f"ERROR: 'training.miniprot_threshold' must be a number (current: {trn['miniprot_threshold']!r}).")
        if not 0 <= val <= 1:
            die(f"ERROR: 'training.miniprot_threshold' must be between 0 and 1 (current: {val})")

    for key in ["flanking_region", "test_size"]:
        if key in trn:
            try:
                val = int(trn[key])
            except (ValueError, TypeError):
                die(f"ERROR: 'training.{key}' must be an integer (current: {trn[key]!r}).")
            if val < 0:
                die(f"ERROR: 'training.{key}' cannot be negative.")

    # 7. Conditional and optional files
    if curation_mode == "user":
        user_gtf = cur.get("user_gtf")
        if not user_gtf:
            die("ERROR: 'curation.mode' is set to 'user', but 'curation.user_gtf' is missing.")
        validate_path(user_gtf, "curation.user_gtf", is_file=True)

    for section, key in OPTIONAL_FILES:
        value = config.get(section, {}).get(key)
        if value:
            validate_path(value, f"{section}.{key}", is_file=True)

    # 8. Augustus species model location
    check_augustus_config_path(prj["toolsdir"])

    logger.info("Input validation completed successfully!")


if __name__ == "__main__":
    # This block is for direct testing. Usually called via sqanti_evidence.
    import yaml
    if len(sys.argv) != 2:
        die("Usage: python input_check.py <config.yaml>")

    with open(sys.argv[1], 'r') as f:
        cfg = yaml.safe_load(f)
    check_inputs(cfg)
