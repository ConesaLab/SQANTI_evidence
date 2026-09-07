sp_name = config.prediction.species



rule extract_rna_hints:
    input:
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.gtf"),
        classification=os.path.join(dir.out.ed_sqanti, f"{sp_name}_classification.txt"),
        hint_config=config.prediction.hint_config,
    output:
        os.path.join(dir.out.ed_hints, f"{sp_name}.rna.hints.gff"),
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


rule align_proteins_miniprot:
    input:
        genome=config.project.genome,
        proteins=os.path.join(dir.out.ab_augustus_model, "sqanti_dominant.faa"),
    output:
        miniprot_gff=os.path.join(dir.out.ed_hints, "sqanti_proteins.miniprot.gff"),
    log:
        os.path.join(dir.logs, "align_proteins_miniprot.log"),
    conda:
        os.path.join(dir.envs, "busco.yaml")
    threads: config.resources.small.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        miniprot_args=config.prediction.miniprot_args,
    shell:
        """
        miniprot --gff -t {threads} {params.miniprot_args} {input.genome} {input.proteins} > {output.miniprot_gff} 2> {log}
        """


rule convert_miniprot_hints:
    input:
        miniprot_gff=os.path.join(dir.out.ed_hints, "sqanti_proteins.miniprot.gff"),
    output:
        protein_hints=os.path.join(dir.out.ed_hints, f"{sp_name}.protein.hints.gff"),
    log:
        os.path.join(dir.logs, "convert_miniprot_hints.log"),
    conda:
        os.path.join(dir.envs, "basic.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        src="P",
        priority=2,
    script:
        os.path.join(dir.scripts, "miniprot_to_hints.py")


rule combine_evidence_hints:
    input:
        rna_hints=os.path.join(dir.out.ed_hints, f"{sp_name}.rna.hints.gff"),
        protein_hints=os.path.join(dir.out.ed_hints, f"{sp_name}.protein.hints.gff"),
    output:
        combined_hints=os.path.join(dir.out.ed_hints, f"{sp_name}.hints.gff"),
    log:
        os.path.join(dir.logs, "combine_evidence_hints.log"),
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    shell:
        """
        cat {input.rna_hints} {input.protein_hints} | sort -k1,1 -k4,4n > {output.combined_hints} 2> {log}
        """


# In split mode the per-chromosome rules in rules/split_augustus.smk (included from the
# snakefile) produce Augustus_prediction.gff instead.
if config.prediction.mode != "split":

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


rule resolve_transcript_tiers:
    input:
        sqanti_gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.gtf"),
        augustus_gff=os.path.join(dir.out.ed_augustus, "Augustus_prediction.gff"),
        hints=os.path.join(dir.out.ed_hints, f"{sp_name}.hints.gff"),
    output:
        resolved_gtf=os.path.join(dir.out.ed_augustus, "resolved_prediction.gtf"),
    log:
        os.path.join(dir.logs, "resolve_transcript_tiers.log"),
    conda:
        os.path.join(dir.envs, "basic.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        min_monoexon_len=300,
        filter_mode=config.prediction.filter_mode,
    script:
        os.path.join(dir.scripts, "resolve_transcript_tiers.py")

