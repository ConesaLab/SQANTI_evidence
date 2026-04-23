import logging
import os
from datetime import datetime
import gzip
import shutil
import glob

# Hardcoded logging configuration
LOG_CONSOLE = True
LOG_FILE = True
LOG_ROTATION = 10  # Keep last 10 log files
LOG_CONFIG_AT_START = True
TRACK_PERFORMANCE = False

def setup_pipeline_logger(log_level="INFO", log_dir="logs"):
    """
    Initialize the pipeline logging system.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_dir: Directory to store log files
    
    Returns:
        logger: Configured logger instance
    """
    # Create logs directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)
    
    # Create logger
    logger = logging.getLogger('pipeline')
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Prevent duplicate handlers if logger already exists
    if logger.handlers:
        return logger
    
    # Create formatters
    detailed_formatter = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler (now with timestamps)
    if LOG_CONSOLE:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(getattr(logging, log_level.upper()))
        console_handler.setFormatter(detailed_formatter)
        logger.addHandler(console_handler)
    
    # File handler
    if LOG_FILE:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_filename = os.path.join(log_dir, f'pipeline_{timestamp}.log')
        
        file_handler = logging.FileHandler(log_filename)
        file_handler.setLevel(logging.DEBUG)  # Always log everything to file
        file_handler.setFormatter(detailed_formatter)
        logger.addHandler(file_handler)
        
        logger.info(f"Log file created: {log_filename}")
    
    return logger


def compress_old_logs(log_dir="logs", keep_last=LOG_ROTATION):
    """
    Compress old log files and keep only the most recent ones.
    
    Args:
        log_dir: Directory containing log files
        keep_last: Number of most recent logs to keep (0 = keep all)
    """
    if keep_last == 0:
        return
    
    # Find all uncompressed log files (exclude current session)
    log_pattern = os.path.join(log_dir, 'pipeline_*.log')
    log_files = sorted(glob.glob(log_pattern), key=os.path.getmtime, reverse=True)
    
    # Skip the most recent log (current session)
    if len(log_files) > 1:
        old_logs = log_files[1:]
        
        for log_file in old_logs:
            gz_file = log_file + '.gz'
            
            # Compress if not already compressed
            if not os.path.exists(gz_file):
                try:
                    with open(log_file, 'rb') as f_in:
                        with gzip.open(gz_file, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    os.remove(log_file)
                except Exception as e:
                    print(f"Warning: Could not compress {log_file}: {e}")
    
    # Clean up old compressed logs
    gz_pattern = os.path.join(log_dir, 'pipeline_*.log.gz')
    gz_files = sorted(glob.glob(gz_pattern), key=os.path.getmtime, reverse=True)
    
    if len(gz_files) > keep_last:
        for old_gz in gz_files[keep_last:]:
            try:
                os.remove(old_gz)
            except Exception as e:
                print(f"Warning: Could not remove {old_gz}: {e}")


def log_config_summary(logger, config):
    """
    Log a summary of the pipeline configuration.
    
    Args:
        logger: Logger instance
        config: Configuration object
    """
    if not LOG_CONFIG_AT_START:
        return
    
    logger.info("=" * 60)
    logger.info("PIPELINE CONFIGURATION SUMMARY")
    logger.info("=" * 60)
    
    # Required parameters
    logger.info("Required Parameters:")
    logger.info(f"  Genome: {config.required.genome}")
    logger.info(f"  Input: {config.required.input}")
    logger.info(f"  Output directory: {config.required.outdir}")
    logger.info(f"  Tools directory: {config.required.toolsdir}")
    
    # Augustus parameters
    logger.info("Augustus Configuration:")
    logger.info(f"  Prediction mode: {config.augustus.prediction}")
    if config.augustus.prediction == "evidence_driven":
        logger.info(f"  Reference GTF: {config.augustus.reference_gtf}")
    logger.info(f"  Species name: {config.augustus.species_name}")
    logger.info(f"  UTR: {config.augustus.utr}")
    logger.info(f"  Mode: {config.augustus.mode}")
    
    # Ab initio parameters
    if config.augustus.prediction == "ab_initio":
        logger.info("Ab Initio Configuration:")
        logger.info(f"  BUSCO lineage: {config.ab_initio.lineage}")
        logger.info(f"  Miniprot threshold: {config.ab_initio.miniprot_threshold}")
        logger.info(f"  Flanking region: {config.ab_initio.flanking_region}")
    
    # QC parameters
    logger.info("Quality Control Configuration:")
    logger.info(f"  OMARK database: {config.qc.omark_db}")
    logger.info(f"  OMARK taxid: {config.qc.omark_taxid}")
    
    logger.info("=" * 60)


def log_phase_transition(logger, phase_name):
    """
    Log the start of a new pipeline phase.
    
    Args:
        logger: Logger instance
        phase_name: Name of the phase
    """
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"PHASE: {phase_name}")
    logger.info("=" * 60)


# Initialize global logger (will be configured in onstart)
pipeline_logger = None
