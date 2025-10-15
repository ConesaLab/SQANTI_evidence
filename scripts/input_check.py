
import os
import sys

def die(msg):
    sys.stderr.write(msg + "\n")
    sys.exit(1)

def warn(msg):
    sys.stderr.write("Warning: " + msg + "\n")

def validate_file(path, name):
    if not os.path.exists(path):
        die(f"Missing required variable  \"{name}\" file or directory: {path}")


def check_inputs(config):
    # Check required files and directories
    for key in ["genome", "input", "toolsdir"]:
        validate_file(config["required"][key], key)


    if os.path.exists(config["required"]["outdir"]):
        warn(f"Output directory {config['required']['outdir']} already exists and may be overwritten.")
    
    if not os.path.exists(config["augustus"].get("reference_gtf", "")):
        warn("Reference GTF file for the organism not found.\nRunning and ab initio prediciton first for SQANTI3")
        config["augustus"]["prediction"] = "ab_initio"
    else:
        config["augustus"]["prediction"] = "evidence_driven"
        
    # Check numerical variables
    if not 0 <= config["ab_initio"]["miniprot_threshold"] <= 1:
        die("miniprot_threshold must be between 0 and 1.")

    if not isinstance(config["ab_initio"]["flanking_region"], int):
        die("flanking_region must be an integer.")

    if not isinstance(config["ab_initio"]["test_size"], int):
        die("test_size must be an integer.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        die("Usage: python input_check.py <config_file>")
    check_inputs(sys.argv[1])
