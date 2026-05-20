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

def die(msg):
    logger.error(msg)
    sys.exit(1)

def warn(msg):
    logger.warning(msg)

def validate_path(path, name, is_file=True):
    """Checks if a path exists and is the correct type."""
    if not os.path.exists(path):
        die(f"ERROR: Missing required {name} at: {path}")
    if is_file and not os.path.isfile(path):
        die(f"ERROR: Expected {name} to be a file, but it is not: {path}")
    if not is_file and not os.path.isdir(path):
        die(f"ERROR: Expected {name} to be a directory, but it is not: {path}")
    logger.info(f"Validated {name}: {path}")

def check_inputs(config):
    """
    Performs hard validation of the configuration before launching Snakemake.
    Focuses on file existence and parameter sanity.
    """
    logger.info("Starting comprehensive input validation...")
    
    # 1. Check for main sections
    required_sections = ["project", "training", "prediction", "curation"]
    for section in required_sections:
        if section not in config:
            die(f"ERROR: Missing required configuration section: '{section}'")

    # 2. Project Section (Hard requirements)
    prj = config["project"]
    for key in ["genome", "input", "toolsdir"]:
        if key not in prj or not prj[key]:
            die(f"ERROR: Missing required variable 'project.{key}'")
    
    validate_path(prj["genome"], "project.genome", is_file=True)
    validate_path(prj["input"], "project.input", is_file=True)
    validate_path(prj["toolsdir"], "project.toolsdir", is_file=False)

    # 3. Output directory check (Warning only)
    outdir = prj.get("outdir", "results")
    if os.path.exists(outdir):
        warn(f"Output directory '{outdir}' already exists. Files may be overwritten.")

    # 4. Numerical Parameter Validation
    # Training
    trn = config["training"]
    if "miniprot_threshold" in trn:
        try:
            val = float(trn["miniprot_threshold"])
            if not 0 <= val <= 1:
                die(f"ERROR: 'training.miniprot_threshold' must be between 0 and 1 (current: {val})")
        except ValueError:
            die(f"ERROR: 'training.miniprot_threshold' must be a number.")

    for key in ["flanking_region", "test_size"]:
        if key in trn:
            try:
                val = int(trn[key])
                if val < 0:
                    die(f"ERROR: 'training.{key}' cannot be negative.")
            except ValueError:
                die(f"ERROR: 'training.{key}' must be an integer.")

    # 5. Conditional File Validation
    # User-defined GTF in curation
    cur = config["curation"]
    if cur.get("mode") == "user":
        user_gtf = cur.get("user_gtf")
        if not user_gtf:
            die("ERROR: 'curation.mode' is set to 'user', but 'curation.user_gtf' is missing.")
        validate_path(user_gtf, "curation.user_gtf", is_file=True)

    # Optional prediction genome
    if prj.get("prediction_genome"):
        validate_path(prj["prediction_genome"], "project.prediction_genome", is_file=True)

    logger.info("Input validation completed successfully!")

if __name__ == "__main__":
    # This block is for direct testing. Usually called via sqanti_evidence.
    import yaml
    if len(sys.argv) != 2:
        die("Usage: python input_check.py <config.yaml>")
    
    with open(sys.argv[1], 'r') as f:
        cfg = yaml.safe_load(f)
    check_inputs(cfg)
