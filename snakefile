
import glob
import attrmap as ap

config = ap.AttrMap(config)

localrules: all, install_tama, install_sqanti

# Setup rules
include: os.path.join("rules","setup","directories.smk")
include: os.path.join("rules","setup","installations.smk")
include: os.path.join("rules","setup","logging_setup.smk")
include: os.path.join("rules","setup","functions.smk")

sample,filetype = get_sample_name(str(config.required.input))
genome_name = get_genome_name(str(config.required.genome))

# Snakemake hooks for logging
onstart:
    global pipeline_logger
    log_level = config.required.log_level 
    if log_level is None:
        log_level = "INFO"
    pipeline_logger = setup_pipeline_logger(log_level=log_level, log_dir=dir.logs)
    
    # Compress old logs
    compress_old_logs(log_dir=dir.logs)
    
    # Log pipeline start
    pipeline_logger.info("=" * 70)
    pipeline_logger.info("GENOME ANNOTATION PIPELINE STARTED")
    pipeline_logger.info("=" * 70)
    
    # Log sample and genome information
    pipeline_logger.info(f"Sample name: {sample}")
    pipeline_logger.info(f"Input file type: {filetype}")
    pipeline_logger.info(f"Genome name: {genome_name}")
    
    # Log configuration summary
    log_config_summary(pipeline_logger, config)
    
    pipeline_logger.info("")
    pipeline_logger.info("Starting workflow execution...")


onsuccess:
    if pipeline_logger:
        pipeline_logger.info("")
        pipeline_logger.info("=" * 70)
        pipeline_logger.info("PIPELINE COMPLETED SUCCESSFULLY")
        pipeline_logger.info("=" * 70)
        pipeline_logger.info(f"Final annotation file: {os.path.join(dir.out.evidence_driven, 'Final_clean_prediction.gff')}")
        pipeline_logger.info(f"GAQET plot: {os.path.join(dir.out.qc_gaqet2, f'{sample}_GAQET.plot.png')}")
        pipeline_logger.info("All output files have been generated.")


onerror:
    if pipeline_logger:
        pipeline_logger.error("")
        pipeline_logger.error("=" * 70)
        pipeline_logger.error("PIPELINE FAILED")
        pipeline_logger.error("=" * 70)
        pipeline_logger.error("An error occurred during pipeline execution.")
        pipeline_logger.error("Check the individual rule logs in the logs/rules/ directory for details.")


include: os.path.join("rules","transcript_modelling.smk")

include: os.path.join("rules","ab_initio.smk")

include: os.path.join("rules","evidence_driven.smk")

include: os.path.join("rules","quality_control.smk")

rule all:
    input:
        os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff"),
        os.path.join(dir.out.qc_gaqet2,f"{sample}_GAQET.plot.png")
