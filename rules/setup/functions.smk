import os
import logging

def get_sqanti_gtf(config):
    logger = logging.getLogger('pipeline')
    if config.augustus.sqanti_gtf_type == "ab_initio":
        result = os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gtf")
        logger.debug(f"Using ab_initio prediction GTF: {result}")
        return result
    elif config.augustus.sqanti_gtf_type == "user_defined":
        result = config.augustus.reference_gtf
        logger.debug(f"Using evidence-driven reference GTF: {result}")
        return result
    elif config.augustus.sqanti_gtf_type == "placebo":
        results = os.path.join(dir.out.ab_initio,"placebo.gtf")
        logger.debug(f"Using placebo GTF: {results}")
        return results

def get_sample_name(file):
    logger = logging.getLogger('pipeline')
    sample = os.path.splitext(os.path.basename(file))[0]
    filetype = os.path.splitext(file)[1]
    logger.debug(f"Detected sample name: {sample}, file type: {filetype}")
    return sample, filetype

def get_pbmm2_input(filetype,config,sample):
    logger = logging.getLogger('pipeline')
    if filetype == ".fastq":
        logger.debug(f"Input is FASTQ, using direct input: {config.required.input}")
        return config.required.input
    else:
        result = os.path.join(dir.out.isoseq,f"{sample}.fastq")
        logger.debug(f"Input is not FASTQ, will use converted file: {result}")
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
