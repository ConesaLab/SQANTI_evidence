# Pure ab initio Augustus prediction (no hints). Only reached when curation.mode == "ab_initio",
# where ab_initio_prediction.gtf serves as the reference annotation for SQANTI3 (get_sqanti_gtf).
# In the default "placebo" mode this file contributes nothing to the DAG.
# In split mode the per-chromosome rules in rules/split_augustus.smk (included from the
# snakefile) produce ab_initio_prediction.gff instead.
if config.prediction.mode != "split":

    rule run_augustus:
        input:
            genome=config.project.genome,
            mod=os.path.join(dir.out.ab_augustus_training, "SC_freq_mod.done"),
        output:
            os.path.join(dir.out.ab_augustus, "ab_initio_prediction.gff"),
        log:
            os.path.join(dir.logs, "run_augustus.log"),
        conda:
            os.path.join(dir.envs, "augustus.yaml")
        resources:
            slurm_extra=f"'--qos={config.resources.big.qos}'",
            cpus_per_task=config.resources.big.cpus,
            mem=config.resources.big.mem,
            runtime=config.resources.big.time,
        params:
            name=config.prediction.species,
        shell:
            "augustus --species={params.name} {input.genome} --protein=on --codingseq=on > {output} 2> {log}"

rule gff2gtf:
    input:
        os.path.join(dir.out.ab_augustus, "ab_initio_prediction.gff"),
    output:
        os.path.join(dir.out.ab_augustus, "ab_initio_prediction.gtf"),
    log:
        os.path.join(dir.logs, "gff2gtf.log"),
    conda:
        os.path.join(dir.envs, "sqanti3.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
        cpus_per_task=config.resources.small.cpus,
    shell:
        "gffread {input} -T -o {output} &> {log}"