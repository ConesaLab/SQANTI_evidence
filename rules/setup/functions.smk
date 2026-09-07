import os
import logging

AB_INITIO_TRAINING_ERROR = (
    "ERROR: 'curation.mode' is 'ab_initio' but 'training.mode' is not 'busco_only'. "
    "To generate an ab initio annotation for SQANTI3 to classify against, the Augustus gene model "
    "must be trained using only BUSCO genes: the SQANTI3-derived training genes ('mixed' or "
    "'sqanti_only') do not exist yet at that point, which makes the workflow circular. "
    "Set training.mode: busco_only, or use curation.mode: placebo / user."
)

def validate_and_fill_config(config_dict):
    """
    Validates hierarchical parameters and fills in default values.
    """
    # 1. project
    if "project" not in config_dict:
        config_dict["project"] = {}
    prj = config_dict["project"]
    if not prj.get("genome"):
        raise ValueError("ERROR: 'project.genome' (reference genome FASTA) is required.")
    if not prj.get("input"):
        raise ValueError("ERROR: 'project.input' (input RNA-Seq data) is required.")
    if not prj.get("toolsdir"):
        raise ValueError("ERROR: 'project.toolsdir' (tools directory) is required.")
        
    prj.setdefault("outdir", "SQANTI_evidence_results")
    prj.setdefault("log_level", "INFO")
    
    if not prj.get("prediction_genome") or not os.path.exists(prj.get("prediction_genome")):
        prj["prediction_genome"] = prj.get("genome")
    
    # 2. training
    if "training" not in config_dict:
        config_dict["training"] = {}
    trn = config_dict["training"]
    trn.setdefault("mode", "mixed")
    if trn["mode"] not in ["mixed", "busco_only", "sqanti_only"]:
        raise ValueError("ERROR: 'training.mode' must be one of 'mixed', 'busco_only', or 'sqanti_only'.")

    if trn["mode"] in ["mixed", "busco_only"]:
        if not trn.get("lineage"):
            raise ValueError("ERROR: 'training.lineage' (BUSCO lineage) is required when training.mode is 'mixed' or 'busco_only'.")
    else:
        trn.setdefault("lineage", "none")
    
    trn.setdefault("skip", False)
    trn.setdefault("miniprot_threshold", 0.95)
    trn.setdefault("flanking_region", 1000)
    trn.setdefault("test_size", 5000)
    trn.setdefault("busco_downloads", "busco_downloads")
    
    # 3. prediction
    if "prediction" not in config_dict:
        config_dict["prediction"] = {}
    prd = config_dict["prediction"]
    if not prd.get("species"):
        raise ValueError("ERROR: 'prediction.species' (Augustus species name) is required.")
        
    prd.setdefault("mode", "full")
    if prd["mode"] not in ["split", "full"]:
        raise ValueError("ERROR: 'prediction.mode' must be either 'split' or 'full'.")
    prd.setdefault("utr", True)
    # Tier 2 (Augustus) single-exon noise filter applied by resolve_transcript_tiers.py
    prd.setdefault("filter_mode", "medium")
    if prd["filter_mode"] not in ["strict", "medium", "none"]:
        raise ValueError("ERROR: 'prediction.filter_mode' must be one of 'strict', 'medium', or 'none'.")

    envs_dir = os.path.abspath(os.path.join(workflow.basedir, "envs"))

    if not prd.get("hint_config") or not os.path.isfile(prd.get("hint_config")):
        prd["hint_config"] = os.path.join(envs_dir, "hint_config.tsv")

    if not prd.get("hint_weights") or not os.path.isfile(prd.get("hint_weights")):
        prd["hint_weights"] = os.path.join(envs_dir, "extrinsic.hints_weights_default.cfg")

    # 4. curation
    if "curation" not in config_dict:
        config_dict["curation"] = {}
    cur = config_dict["curation"]
    cur.setdefault("data_type","pacbio")
    if not cur.get("filter_rules") or not os.path.isfile(cur.get("filter_rules")):
        cur["filter_rules"] = os.path.join(envs_dir, "filter_rules.json")
        
    cur.setdefault("mode", "placebo")
    if cur["mode"] not in ["placebo", "user", "ab_initio"]:
        raise ValueError("ERROR: 'curation.mode' must be one of 'placebo', 'user', or 'ab_initio'.")
    cur.setdefault("user_gtf", "")
    if cur["mode"] == "user" and not cur["user_gtf"]:
        raise ValueError("ERROR: 'curation.mode' is 'user' but 'curation.user_gtf' is empty.")
    # The ab initio reference for SQANTI3 is predicted with the trained Augustus model, and
    # mixed/sqanti_only training needs SQANTI3's output first -> cyclic DAG.
    if cur["mode"] == "ab_initio" and trn["mode"] != "busco_only":
        raise ValueError(AB_INITIO_TRAINING_ERROR)
    
    # 5. evaluation
    if "evaluation" not in config_dict:
        config_dict["evaluation"] = {}
    eva = config_dict["evaluation"]
    eva.setdefault("omark_db", "LUCA")
    eva.setdefault("omark_taxid", 10090)
    eva.setdefault("reference_gtf", "")
 
    # 6. resources
    # Convention: Slurm-style human-readable values, e.g. mem: "8GB", time: "2h" (Snakemake parses
    # both via humanfriendly for the `mem` and `runtime` resources). Missing keys are filled per tier.
    if "resources" not in config_dict:
        config_dict["resources"] = {}
    res = config_dict["resources"]
    resource_defaults = {
        "small":  {"cpus": 2,  "mem": "8GB",  "time": "2h",  "qos": "short"},
        "medium": {"cpus": 8,  "mem": "12GB", "time": "10h", "qos": "short"},
        "big":    {"cpus": 10, "mem": "20GB", "time": "24h", "qos": "long"},
        "busco":  {"cpus": 30, "mem": "60GB", "time": "72h", "qos": "long"},
    }
    for tier, defaults in resource_defaults.items():
        tier_cfg = res.setdefault(tier, {})
        for key, value in defaults.items():
            tier_cfg.setdefault(key, value)

    return config_dict

def get_sqanti_gtf(config):
    logger = logging.getLogger('pipeline')
    if config.curation.mode == "ab_initio":
        result = os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
        logger.debug(f"Using ab_initio prediction GTF: {result}")
        return result
    elif config.curation.mode == "user":
        result = config.curation.user_gtf
        logger.debug(f"Using user-defined reference GTF: {result}")
        return result
    elif config.curation.mode == "placebo":
        results = os.path.join(dir.out.ab_initio,"placebo.gtf")
        logger.debug(f"Using placebo GTF: {results}")
        return results

def get_sample_name(file):
    logger = logging.getLogger('pipeline')
    sample = os.path.splitext(os.path.basename(file))[0]
    filetype = os.path.splitext(file)[1]
    logger.debug(f"Detected sample name: {sample}, file type: {filetype}")
    return sample, filetype

def get_pbmm2_input(filetype, config, sample):
    logger = logging.getLogger('pipeline')
    if filetype == ".bam":
        logger.debug(f"Input is BAM, using direct input: {config.project.input}")
        return config.project.input
    else:
        result = os.path.join(dir.out.isoquant, f"{sample}.bam")
        logger.debug(f"Input is not BAM, will use converted file: {result}")
        return result

def get_chromosomes(file):
    logger = logging.getLogger('pipeline')
    chromosomes = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chromosomes.append(line[1:].split()[0])
    logger.debug(f"Found {len(chromosomes)} chromosomes in {file}")
    return chromosomes

def get_genome_name(file):
    logger = logging.getLogger('pipeline')
    genome_name = os.path.splitext(os.path.basename(file))[0]
    logger.debug(f"Extracted genome name: {genome_name}")
    return genome_name

def check_augustus_species(species):
    """
    Checks if the augustus species exists.
    """
    import subprocess
    try:
        config_path = os.environ.get('AUGUSTUS_CONFIG_PATH')
        if config_path:
            species_dir = os.path.join(config_path, "species", species)
            if os.path.isdir(species_dir):
                return True
        res = subprocess.run(["augustus", f"--species={species}"], capture_output=True, text=True)
        if "unknown species" in res.stderr.lower() or "not found" in res.stderr.lower():
            return False
        return True
    except Exception:
        return True

def get_isoquant_input_flag(reads_path):
    path_lower = str(reads_path).lower()
    if path_lower.endswith('.gz'):
        path_lower = path_lower[:-3]
        
    if path_lower.endswith('.bam'):
        return f"--unmapped_bam {reads_path}"
    elif any(path_lower.endswith(ext) for ext in ['.fastq', '.fq', '.fasta', '.fa']):
        return f"--fastq {reads_path}"
    else:
        raise ValueError(f"Unknown input format for IsoQuant: {reads_path}")
