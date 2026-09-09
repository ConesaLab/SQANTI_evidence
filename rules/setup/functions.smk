import os
import sys
import logging

# Shared validation tables/messages live in scripts/input_check.py (plain Python, no Snakemake
# dependency) so the wrapper pre-flight and this parse-time validation cannot drift apart.
sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from input_check import ALLOWED_VALUES, validate_modes, validate_isoquant_args

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
    if trn["mode"] not in ALLOWED_VALUES[("training", "mode")]:
        raise ValueError(f"ERROR: 'training.mode' must be one of {ALLOWED_VALUES[('training', 'mode')]}.")

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
    if prd["mode"] not in ALLOWED_VALUES[("prediction", "mode")]:
        raise ValueError(f"ERROR: 'prediction.mode' must be one of {ALLOWED_VALUES[('prediction', 'mode')]}.")
    prd.setdefault("utr", True)
    # Tier 2 (Augustus) single-exon noise filter applied by resolve_transcript_tiers.py
    prd.setdefault("filter_mode", "medium")
    if prd["filter_mode"] not in ALLOWED_VALUES[("prediction", "filter_mode")]:
        raise ValueError(f"ERROR: 'prediction.filter_mode' must be one of {ALLOWED_VALUES[('prediction', 'filter_mode')]}.")

    envs_dir = os.path.abspath(os.path.join(workflow.basedir, "envs"))

    if not prd.get("hint_config") or not os.path.isfile(prd.get("hint_config")):
        prd["hint_config"] = os.path.join(envs_dir, "hint_config.tsv")

    if not prd.get("hint_weights") or not os.path.isfile(prd.get("hint_weights")):
        prd["hint_weights"] = os.path.join(envs_dir, "extrinsic.hints_weights_default.cfg")

    # Miniprot options for the self-protein alignment that feeds src=P hints. Miniprot only *writes*
    # secondary alignments scoring >= --outs * best (default 0.99); --outs=0.5 lets paralog loci receive
    # protein hints (validated on Arabidopsis 2026-09-07: locus Sn/Sp 55.9/74.2 -> 57.7/76.2).
    prd.setdefault("miniprot_args", "-N 30 -p 0.6 --outs=0.5")

    # 4. curation
    if "curation" not in config_dict:
        config_dict["curation"] = {}
    cur = config_dict["curation"]
    cur.setdefault("data_type","pacbio")
    # Extra IsoQuant options appended verbatim (e.g. "--polya_trimmed all --stranded forward" for
    # Iso-Seq FLNC reads, whose poly(A) tails were removed by `isoseq refine`; without a tail IsoQuant
    # cannot build novel mono-exon transcripts). Options the rule sets itself are rejected by input_check.
    cur.setdefault("isoquant_args", "")
    iq_error = validate_isoquant_args(cur["isoquant_args"])
    if iq_error:
        raise ValueError(iq_error)
    if cur["data_type"] not in ALLOWED_VALUES[("curation", "data_type")]:
        raise ValueError(f"ERROR: 'curation.data_type' must be one of {ALLOWED_VALUES[('curation', 'data_type')]}.")
    if not cur.get("filter_rules") or not os.path.isfile(cur.get("filter_rules")):
        cur["filter_rules"] = os.path.join(envs_dir, "filter_rules.json")
        
    cur.setdefault("mode", "placebo")
    if cur["mode"] not in ALLOWED_VALUES[("curation", "mode")]:
        raise ValueError(f"ERROR: 'curation.mode' must be one of {ALLOWED_VALUES[('curation', 'mode')]}.")
    cur.setdefault("user_gtf", "")
    if cur["mode"] == "user" and not cur["user_gtf"]:
        raise ValueError("ERROR: 'curation.mode' is 'user' but 'curation.user_gtf' is empty.")
    mode_error = validate_modes(trn["mode"], cur["mode"])
    if mode_error:
        raise ValueError(mode_error)
    
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
