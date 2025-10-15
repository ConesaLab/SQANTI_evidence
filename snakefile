import glob
import attrmap as ap

config = ap.AttrMap(config)

localrules: all, install_tama, install_sqanti
# Setup rules
include: os.path.join("rules","setup","directories.smk")
include: os.path.join("rules","setup","installations.smk")
include: os.path.join("rules","setup","functions.smk")

sample,filetype = get_sample_name(config.required.input)
genome_name = get_genome_name(config.required.genome)


include: os.path.join("rules","transcript_modelling.smk")

include: os.path.join("rules","ab_initio.smk")

include: os.path.join("rules","evidence_driven.smk")

include: os.path.join("rules","quality_control.smk")

rule all:
    input:
        os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff"),
        os.path.join(dir.out.qc_omark,"Final.pdf"),
        os.path.join(dir.out.qc_busco),
        os.path.join(dir.out.qc_agat,"Final_stats.txt"),
