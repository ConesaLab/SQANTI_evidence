import os
import logging

def validate_and_fill_config(config_dict):
    """
    Validates required parameters and fills in default values for optional parameters.
    Returns the modified config dictionary.
    """
    # 1. required
    if "required" not in config_dict:
        config_dict["required"] = {}
    req = config_dict["required"]
    if not req.get("genome"):
        raise ValueError("ERROR: 'required.genome' (path to the reference genome FASTA) must be provided in the config.")
    if not req.get("input"):
        raise ValueError("ERROR: 'required.input' (path to the input lrRNA data) must be provided in the config.")
    if not req.get("toolsdir"):
        raise ValueError("ERROR: 'required.toolsdir' (path to the tools directory) must be provided in the config.")
        
    req.setdefault("outdir", "results")
    req.setdefault("log_level", "INFO")
    
    # 2. augustus
    if "augustus" not in config_dict:
        config_dict["augustus"] = {}
    aug = config_dict["augustus"]
    if not aug.get("species_name"):
        raise ValueError("ERROR: 'augustus.species_name' must be provided in the config.")
        
    aug.setdefault("utr", True)
    aug.setdefault("mode", "full")
    aug.setdefault("filter_mode", "monoexon")
    aug.setdefault("config", None)
    
    # 3. sqanti
    if "sqanti" not in config_dict:
        config_dict["sqanti"] = {}
    sq = config_dict["sqanti"]
    if not sq.get("json_rules"):
        raise ValueError("ERROR: 'sqanti.json_rules' must be provided in the config.")
        
    sq.setdefault("gtf_type", "placebo")
    sq.setdefault("reference_gtf", "")
    
    # 4. ab_initio
    if "ab_initio" not in config_dict:
        config_dict["ab_initio"] = {}
    ab = config_dict["ab_initio"]
    if not ab.get("lineage"):
        raise ValueError("ERROR: 'ab_initio.lineage' (BUSCO lineage) must be provided in the config.")
        
    ab.setdefault("miniprot_threshold", 0.95)
    ab.setdefault("flanking_region", 1000)
    ab.setdefault("test_size", 5000)
    
    if "busco_downloads" in ab:
        del ab["busco_downloads"]

    # 5. qc
    if "qc" not in config_dict:
        config_dict["qc"] = {}
    qc = config_dict["qc"]
    qc.setdefault("omark_db", "LUCA")
    qc.setdefault("omark_taxid", 10090)
    
    # 6. resources
    if "resources" not in config_dict:
        config_dict["resources"] = {}
    res = config_dict["resources"]
    res.setdefault("small", {"cpus": 2, "mem_mb": 8192, "time_min": 120, "qos": "short"})
    res.setdefault("medium", {"cpus": 8, "mem_mb": 12288, "time_min": 600, "qos": "short"})
    res.setdefault("big", {"cpus": 10, "mem_mb": 20480, "time_min": 1440, "qos": "long"})
    res.setdefault("small_bigMem", {"mem_mb": 20480})
    res.setdefault("busco", {"cpus": 30, "mem_mb": 61440, "time_min": 7200, "qos": "long"})

    return config_dict

def get_sqanti_gtf(config):
    logger = logging.getLogger('pipeline')
    if config.sqanti.gtf_type == "ab_initio":
        result = os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
        logger.debug(f"Using ab_initio prediction GTF: {result}")
        return result
    elif config.sqanti.gtf_type == "user_defined":
        result = config.sqanti.reference_gtf
        logger.debug(f"Using evidence-driven reference GTF: {result}")
        return result
    elif config.sqanti.gtf_type == "placebo":
        results = os.path.join(dir.out.ab_initio,"placebo.gtf")
        logger.debug(f"Using placebo GTF: {results}")
        return results

def get_sample_name(file):
    logger = logging.getLogger('pipeline')
    sample = os.path.splitext(os.path.basename(file))[0]
    filetype = os.path.splitext(file)[1]
    logger.debug(f"Detected sample name: {sample}, file type: {filetype}")
    return sample, filetype

#TODO: Change this so bam is predefined input
def get_pbmm2_input(filetype,config,sample):
    logger = logging.getLogger('pipeline')
    if filetype == ".bam":
        logger.debug(f"Input is BAM, using direct input: {config.required.input}")
        return config.required.input
    else:
        result = os.path.join(dir.out.isoseq,f"{sample}.bam")
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
