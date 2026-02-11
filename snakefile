import glob
import attrmap as ap

config = ap.AttrMap(config)

localrules: all, install_tama, install_sqanti
# Setup rules
include: os.path.join("rules","setup","directories.smk")
include: os.path.join("rules","setup","installations.smk")
include: os.path.join("rules","setup","functions.smk")

sample,filetype = get_sample_name(str(config.required.input))
genome_name = get_genome_name(str(config.required.genome))


include: os.path.join("rules","transcript_modelling.smk")

include: os.path.join("rules","ab_initio.smk")

include: os.path.join("rules","evidence_driven.smk")

include: os.path.join("rules","quality_control.smk")

# This rule eliminates more than 50% of the size of the final directory.
# All of the files are non-essential, and snakemake will not re-run any of the previous steps if they are deleted, so this is a safe way to save space.
rule cleanup:
    input:
        os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff")
    output:
        temp(touch(os.path.join(dir.out.base,"cleanup.done")))
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "cleanup.log")
    shell:
        """
        rm -rf {dir.out.ab_busco}/logs \
               {dir.out.ab_busco}/run_*/miniprot_output/ref.mpi \
               {dir.out.ed_sqanti}/TD2 \
               {dir.out.isoseq_collapsed}/{sample}.collapsed.fastq 
        """

rule all:
    input:
        os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff"),
        os.path.join(dir.out.base,"cleanup.done"),
        # os.path.join(dir.out.qc_omark,"Final.pdf"),
        # os.path.join(dir.out.qc_busco),
        # os.path.join(dir.out.qc_agat,"Final_stats.txt"),
        os.path.join(dir.out.qc_gaqet2,f"{sample}_GAQET.stats.tsv")
