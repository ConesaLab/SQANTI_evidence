
import os
import sys
import logging

# Use the pipeline logger if available, otherwise create a basic one
logger = logging.getLogger('pipeline')
if not logger.handlers:
    # Fallback for standalone execution
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

def validate_file(path, name):
    if not os.path.exists(path):
        die(f"Missing required variable  \"{name}\" file or directory: {path}")
    else:
        logger.info(f"Validated {name}: {path}")


def check_inputs(config):
    logger.info("Starting input validation...")
    
    # Check required files and directories
    for key in ["genome", "input", "toolsdir"]:
        validate_file(config["required"][key], key)


    if os.path.exists(config["required"]["outdir"]):
        warn(f"Output directory {config['required']['outdir']} already exists and may be overwritten.")
    
    if not os.path.exists(config["augustus"].get("reference_gtf", "")):
        warn("Reference GTF file for the organism not found.")
        warn("Running ab initio prediction first for SQANTI3")
        config["augustus"]["prediction"] = "ab_initio"
    else:
        logger.info(f"Validated reference GTF: {config['augustus']['reference_gtf']}")
        config["augustus"]["prediction"] = "evidence_driven"
        
    # Check numerical variables
    if not 0 <= config["ab_initio"]["miniprot_threshold"] <= 1:
        die("miniprot_threshold must be between 0 and 1.")
    else:
        logger.info(f"Validated miniprot_threshold: {config['ab_initio']['miniprot_threshold']}")

    if not isinstance(config["ab_initio"]["flanking_region"], int):
        die("flanking_region must be an integer.")
    else:
        logger.info(f"Validated flanking_region: {config['ab_initio']['flanking_region']}")

    if not isinstance(config["ab_initio"]["test_size"], int):
        die("test_size must be an integer.")
    else:
        logger.info(f"Validated test_size: {config['ab_initio']['test_size']}")
    
    logger.info("Input validation completed successfully!")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        die("Usage: python input_check.py <config_file>")
    check_inputs(sys.argv[1])
