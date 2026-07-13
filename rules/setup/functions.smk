import os
import logging

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
    if not trn.get("lineage"):
        raise ValueError("ERROR: 'training.lineage' (BUSCO lineage) is required.")
    
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
    prd.setdefault("utr", True)
    prd.setdefault("filter_mode", "monoexon")

    if prd["filter_mode"] not in ["monoexon", "all"]:
        raise ValueError("ERROR: 'prediction.filter_mode' must be either 'monoexon' or 'all'.")

    envs_dir = os.path.abspath(os.path.join(workflow.basedir, "envs"))

    if not prd.get("hint_config") or not os.path.isfile(prd.get("hint_config")):
        prd["hint_config"] = os.path.join(envs_dir, "hint_config.tsv")
        
    if not prd.get("hint_weights") or not os.path.isfile(prd.get("hint_weights")):
        prd["hint_weights"] = os.path.join(envs_dir, "extrinsic.hints_weights_default.cfg")
    
    # 4. curation
    if "curation" not in config_dict:
        config_dict["curation"] = {}
    cur = config_dict["curation"]

    if not cur.get("filter_rules") or not os.path.isfile(cur.get("filter_rules")):
        cur["filter_rules"] = os.path.join(envs_dir, "filter_rules.json")
        
    cur.setdefault("mode", "placebo")
    cur.setdefault("user_gtf", "")
    
    # 5. evaluation
    if "evaluation" not in config_dict:
        config_dict["evaluation"] = {}
    eva = config_dict["evaluation"]
    eva.setdefault("omark_db", "LUCA")
    eva.setdefault("omark_taxid", 10090)
    eva.setdefault("reference_gtf", "")
    
    # 6. resources
    if "resources" not in config_dict:
        config_dict["resources"] = {}
    res = config_dict["resources"]
    res.setdefault("small", {"cpus": 2, "mem_mb": 8192, "time_min": 120, "qos": "short"})
    res.setdefault("medium", {"cpus": 8, "mem_mb": 12288, "time_min": 600, "qos": "short"})
    res.setdefault("big", {"cpus": 10, "mem_mb": 20480, "time_min": 1440, "qos": "long"})
    res.setdefault("small_bigMem", {"mem_mb": 20480})
    res.setdefault("busco", {"cpus": 30, "mem_mb": 61440, "time_min": 7200, "qos": "long"})

    # Backward compatibility for mem and time keys
    for rname, rdict in res.items():
        if "mem_mb" in rdict:
            rdict["mem"] = f"{rdict['mem_mb']}M"
        if "time_min" in rdict:
            # Convert minutes to HH:MM:SS if needed, or just keep as minutes string
            # Most Slurm executors accept minutes, but Snakemake's 'runtime' often expects minutes
            rdict["time"] = f"{rdict['time_min']}"

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
        result = os.path.join(dir.out.isoseq, f"{sample}.bam")
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
