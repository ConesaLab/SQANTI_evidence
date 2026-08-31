sp_name = config.prediction.species



rule extract_hints:
    input:
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.gtf"),
        classification=os.path.join(dir.out.ed_sqanti, f"{sp_name}_classification.txt"),
        hint_config=config.prediction.hint_config,
    output:
        os.path.join(dir.out.ed_hints, f"{sp_name}.hints.gff"),
    conda:
        os.path.join(dir.envs, "busco.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        utr=config.prediction.utr,
    script:
        os.path.join(dir.scripts, "generate_hints.py")


if config.prediction.mode == "split":

    include: "split_augustus.smk"

else:

    rule augustus_hints:
        input:
            genome=config.project.prediction_genome,
            mod=os.path.join(dir.out.ab_augustus_training, "SC_freq_mod.done"),
            gff=os.path.join(dir.out.ed_hints, f"{sp_name}.hints.gff"),
        output:
            os.path.join(dir.out.ed_augustus, "Augustus_prediction.gff"),
        log:
            os.path.join(dir.logs, "run_augustus_ed.log"),
        conda:
            os.path.join(dir.envs, "augustus.yaml")
        resources:
            slurm_extra=f"'--qos={config.resources.big.qos}'",
            cpus_per_task=config.resources.big.cpus,
            mem=config.resources.big.mem,
            runtime=config.resources.big.time,
        params:
            name=config.prediction.species,
            extcfg=config.prediction.hint_weights,
        shell:
            """
            augustus --species={params.name} {input.genome} --hintsfile={input.gff} \
            --extrinsicCfgFile={params.extcfg} --protein=on --codingseq=on \
            --alternatives-from-evidence=true > {output} 2> {log}
            """


rule filter_monoexons:
    input:
        gff=os.path.join(dir.out.ed_augustus, "Augustus_prediction.gff"),
        hints=os.path.join(dir.out.ed_hints, f"{sp_name}.hints.gff"),
    output:
        gff=os.path.join(dir.out.ed_augustus, "Augustus_prediction.filtered.gff"),
    log:
        os.path.join(dir.logs, "filter_monoexons.log"),
    conda:
        os.path.join(dir.envs, "busco.yaml")
    threads: config.resources.small.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        filter_mode=config.prediction.filter_mode,
    script:
        os.path.join(dir.scripts, "filter_monoexons.py")
