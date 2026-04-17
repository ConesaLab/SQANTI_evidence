sp_name=config.augustus.species_name

rule run_sqanti:
    input:
        isoforms = os.path.join(dir.out.isoseq_collapsed,f"{sample}.collapsed.gff"),
        ref_gff = get_sqanti_gtf(config),
        ref_genome = config.required.genome,
        sqanti = os.path.join(dir.tools_sqanti,"sqanti_installed.done")
    output:
        classification = os.path.join(dir.out.ed_sqanti,f"{sp_name}_classification.txt"),
        gtf = os.path.join(dir.out.ed_sqanti,f"{sp_name}_corrected.cds.gtf"),
        td2 = temp(directory(os.path.join(dir.out.ed_sqanti,"TD2")))
    threads:
        config.resources.medium.cpus
    conda:
        f"{dir.envs}/sqanti3.yaml"
    params:
        sp_name = sp_name
    log:
        os.path.join(dir.logs,"run_sqanti.log")
    resources:
        slurm_extra = f"\'--qos={config.resources.medium.qos}\'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.big.mem,
        runtime = config.resources.medium.time
    shell:
        #TODO: Eliminate this for the final release, as it is only used in Garnatxa
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        python {dir.tools_sqanti}/sqanti3_qc.py --isoforms {input.isoforms} --refGTF {input.ref_gff} --refFasta {input.ref_genome} \
            --dir {dir.out.ed_sqanti} --output {params.sp_name} -t {threads} --include_ORF &> {log}
        mv {dir.out.ed_sqanti}/{params.sp_name}_corrected.cds.gff3 {output.gtf}
        """

rule filter_isoforms:
    input:
        classification = os.path.join(dir.out.ed_sqanti,f"{sp_name}_classification.txt"),
        gtf = os.path.join(dir.out.ed_sqanti,f"{sp_name}_corrected.cds.gtf")
    output:
        gtf = os.path.join(dir.out.ed_sqanti,f"{sp_name}.filtered.gtf")
    conda:
        os.path.join(dir.envs,"sqanti3.yaml")
    log:
        os.path.join(dir.logs,"filter_sqanti.log")
    threads:
        config.resources.small.cpus
    params:
        json_rules = config.sqanti.json_rules,
        sp_name = sp_name
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime = config.resources.small.time
    shell:
        """
        export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        python {dir.tools_sqanti}/sqanti3_filter.py rules --sqanti_class {input.classification} --filter_gtf {input.gtf} \
            -j {params.json_rules} --dir {dir.out.ed_sqanti} \
            --output {params.sp_name} &> {log}
        """

rule extract_hints:
    input:
        gtf = os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.gtf"),
        classification = os.path.join(dir.out.ed_sqanti,f"{sp_name}_classification.txt")
    output:
        os.path.join(dir.out.ed_hints, f"{sp_name}.hints.gff")
    conda:
        os.path.join(dir.envs,"busco.yaml")
    params:
        utr = config.augustus.utr
        #TODO: Perhaps add techonolgy and priority options
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime = config.resources.small.time
    script:
        os.path.join(dir.scripts,"generate_hints.py")

if config.augustus.mode == "split":
    include: "split_augustus.smk"

else:
    rule augustus_hints:
        input:
            genome = config.required.genome,
            mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done"),
            gff = os.path.join(dir.out.ed_hints,f"{sp_name}.hints.gff")
        output:
            os.path.join(dir.out.ed_augustus,"Augustus_prediction.gff")
        conda:
            os.path.join(dir.envs,"augustus.yaml")
        params:
            name = config.augustus.species_name,
            extcfg = config.augustus.config if config.evidence_driven.config else f"{dir.envs}/extrinsic.M.RM.PB.cfg"
        log:
            os.path.join(dir.logs,"run_augustus_ed.log")
        resources:
            slurm_extra = f"\'--qos={config.resources.big.qos}\'",
            cpus_per_task = config.resources.big.cpus,
            mem = config.resources.big.mem,
            runtime = config.resources.big.time
        shell:
            """
            augustus --species={params.name} {input.genome} --hintsfile={input.gff} \
            --extrinsicCfgFile={params.extcfg} --protein=on --codingseq=on > {output} 2> {log}
            """

rule filter_monoexons:
    input:
        gff = os.path.join(dir.out.ed_augustus,"Augustus_prediction.gff"),
        hints = os.path.join(dir.out.ed_hints,f"{sp_name}.hints.gff")
    output:
        gff = os.path.join(dir.out.ed_augustus,"Augustus_prediction.filtered.gff")
    conda:
        os.path.join(dir.envs,"busco.yaml")
    log:
        os.path.join(dir.logs,"filter_monoexons.log")
    threads:
        config.resources.small.cpus
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime = config.resources.small.time
    script:
        os.path.join(dir.scripts,"filter_monoexons.py")
