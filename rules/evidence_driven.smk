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
    shell:
        """
        miniprot --gff -t {threads} -N 5 -p 0.6 {input.genome} {input.proteins} > {output.miniprot_gff} 2> {log}
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


rule run_tsebra:
    input:
        sqanti_gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.gtf"),
        augustus_gff=os.path.join(dir.out.ed_augustus, "Augustus_prediction.gff"),
        hints=os.path.join(dir.out.ed_hints, f"{sp_name}.hints.gff"),
        cfg=config.prediction.tsebra_config,
    output:
        tsebra_gtf=os.path.join(dir.out.ed_augustus, "tsebra_prediction.gtf"),
    log:
        os.path.join(dir.logs, "run_tsebra.log"),
    conda:
        os.path.join(dir.envs, "tsebra.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    shell:
        """
        tsebra.py -g {input.augustus_gff} \
                  -k {input.sqanti_gtf} \
                  -e {input.hints} \
                  -c {input.cfg} \
                  -o {output.tsebra_gtf} &> {log}
        """
