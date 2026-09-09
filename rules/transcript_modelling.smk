sp_name = config.prediction.species


rule run_isoquant:
    input:
        reads = config.project.input,
        ref = config.project.genome,
    output:
        gtf = os.path.join(dir.out.isoquant,sample, f"{sample}.transcript_models.gtf")
    conda:
        f"{dir.envs}/isoquant.yaml"
    threads:
        config.resources.big.cpus
    resources:
        cpus_per_task = config.resources.big.cpus,
        slurm_extra = f"\'--qos={config.resources.big.qos}\'",
        mem = config.resources.big.mem,
        runtime = config.resources.big.time
    params:
        input_flag = lambda wildcards, input: get_isoquant_input_flag(input.reads),
        outdir = dir.out.isoquant,
        data_type = config.curation.data_type,
        prefix = sample,
        extra = config.curation.isoquant_args
    log:
        os.path.join(dir.logs, "isoquant.log")
    shell:
        """
        isoquant \
            --reference {input.ref} \
            {params.input_flag} \
            --data_type {params.data_type} \
            --prefix {params.prefix} \
            --threads {threads} \
            {params.extra} \
            -o {params.outdir} &> {log}
        """

rule run_sqanti:
    input:
        isoforms=os.path.join(dir.out.isoquant, sample, f"{sample}.transcript_models.gtf"),
        ref_gff=get_sqanti_gtf(config),
        ref_genome=config.project.genome,
    output:
        classification=os.path.join(dir.out.ed_sqanti, f"{sp_name}_classification.txt"),
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}_corrected.cds.gtf"),
    log:
        os.path.join(dir.logs, "run_sqanti.log"),
    conda:
        f"{dir.envs}/sqanti3.yaml"
    threads: config.resources.medium.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.medium.qos}'",
        cpus_per_task=config.resources.medium.cpus,
        mem=config.resources.big.mem,
        runtime=config.resources.medium.time,
    params:
        sp_name=sp_name,
        fl_matrix=os.path.join(dir.out.isoquant, sample, f"{sample}.discovered_transcript_counts.tsv"),
    shell:
        """
        sqanti3_qc.py --isoforms {input.isoforms} --refGTF {input.ref_gff} --refFasta {input.ref_genome} \
            --dir {dir.out.ed_sqanti} --output {params.sp_name} -t {threads} --include_ORF -fl {params.fl_matrix} --report skip &> {log}
        mv {dir.out.ed_sqanti}/{params.sp_name}_corrected.cds.gff3 {output.gtf}
        """


rule filter_isoforms:
    input:
        classification=os.path.join(dir.out.ed_sqanti, f"{sp_name}_classification.txt"),
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}_corrected.cds.gtf"),
    output:
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.gtf"),
        classif=os.path.join(dir.out.ed_sqanti, f"{sp_name}_RulesFilter_classification.txt")
    log:
        os.path.join(dir.logs, "filter_sqanti.log"),
    conda:
        os.path.join(dir.envs, "sqanti3.yaml")
    threads: config.resources.small.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        json_rules=config.curation.filter_rules,
        sp_name=sp_name,
    shell:
        """
        #export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
        sqanti3_filter.py rules --sqanti_class {input.classification} --filter_gtf {input.gtf} \
            -j {params.json_rules} --dir {dir.out.ed_sqanti} --skip_report \
            --output {params.sp_name} &> {log}
        """


rule restore_gene_ids:
    """Put the transcript-model gene grouping back after SQANTI3.

    SQANTI3 gives every intergenic isoform its own novelGene id; against the placebo reference
    that is every isoform, so the filtered GTF has one gene per transcript. The IsoQuant
    gene_id is restored per transcript_id (GTF and RulesFilter classification), so gene
    counts and the dominant-isoform selection see real genes again.
    """
    input:
        models=os.path.join(dir.out.isoquant, sample, f"{sample}.transcript_models.gtf"),
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.gtf"),
        classification=os.path.join(dir.out.ed_sqanti, f"{sp_name}_RulesFilter_classification.txt"),
    output:
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.regrouped.gtf"),
        classification=os.path.join(dir.out.ed_sqanti, f"{sp_name}_RulesFilter_classification.regrouped.txt"),
    log:
        os.path.join(dir.logs, "restore_gene_ids.log"),
    conda:
        os.path.join(dir.envs, "basic.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    script:
        os.path.join(dir.scripts, "restore_gene_ids.py")
