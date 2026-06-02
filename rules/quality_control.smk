localrules: gaqet2_setup

rule agat_cleaning:
    input:
        os.path.join(dir.out.ed_augustus,"Augustus_prediction.filtered.gff")
    output:
        os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff")
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.medium.mem,
        runtime = config.resources.small.time
    log:
        os.path.join(dir.logs, "agat_cleaning.log")
    conda:
        os.path.join(dir.envs, "gaqet2.yaml")
    threads:
        config.resources.small.cpus
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input} -o {output} &> {log}
        """
rule gaqet2_setup:
    input:
        genome = config.project.prediction_genome,
        annotation = os.path.join(dir.out.ed_augustus,"Augustus_prediction.filtered.gff"),
    output:
        os.path.join(dir.out.qc_gaqet2,"gaqet2_config.yaml")
    conda:
        os.path.join(dir.envs, "gaqet2.yaml")
    params:
        config = os.path.join(dir.envs, "gaqet2_conf.yaml"),
        id = sample,
        outdir = dir.out.qc_gaqet2,
        busco_lineage = config.training.lineage,
        taxid = config.evaluation.omark_taxid,
        omark_db = os.path.join(dir.tools_omark,f"{config.evaluation.omark_db}.h5"),
        threads = config.resources.medium.cpus
    log:
        os.path.join(dir.logs, "gaqet2_setup.log")
    script:
        f"{dir.scripts}/gaqet2_setup.py"

rule gaqet2:
    input:
        config = os.path.join(dir.out.qc_gaqet2,"gaqet2_config.yaml"),
    output:
        os.path.join(dir.out.qc_gaqet2,f"{sample}_GAQET.stats.tsv")
    threads:
        config.resources.medium.cpus
    conda:
        os.path.join(dir.envs, "gaqet2.yaml")
    log:
        os.path.join(dir.logs, "gaqet2.log")
    resources:
        slurm_extra = f"\'--qos={config.resources.medium.qos}\'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.big.mem,
        runtime = config.resources.medium.time
    shell:
        """
        GAQET --yaml {input.config} &> {log}
        """

rule gaqet2_plot:
    input:
        os.path.join(dir.out.qc_gaqet2,f"{sample}_GAQET.stats.tsv")
    output:
        os.path.join(dir.out.qc_gaqet2,f"{sample}_GAQET.plot.png")
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime = config.resources.small.time
    threads:
        config.resources.small.cpus
    conda:
        os.path.join(dir.envs, "gaqet2.yaml")
    log:
        os.path.join(dir.logs, "gaqet2_plot.log")
    shell:
        """
        GAQET_PLOT -i {input} -o {output} &> {log}
        """

rule gffcompare_eval:
    input:
        ref = config.evaluation.reference_gtf,
        anno = os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff")
    output:
        stats = os.path.join(dir.out.qc_gffcompare, f"{sample}.stats") 
    params:
        out_prefix = os.path.join(dir.out.qc_gffcompare, f"{sample}")
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime = config.resources.small.time
    conda:
        os.path.join(dir.envs, "gffcompare.yaml")
    log:
        os.path.join(dir.logs, "gffcompare.log")
    shell:
        "gffcompare -r {input.ref} -o {params.out_prefix} {input.anno} -T &> {log}"

rule subset_reference_cds:
    input:
        ref = config.evaluation.reference_gtf
    output:
        ref_cds = os.path.join(dir.out.qc_gffcompare, "reference.cds.gtf")
    log:
        os.path.join(dir.logs, "subset_reference_cds.log")
    conda:
        os.path.join(dir.envs, "busco.yaml")
    shell:
        "python3 {dir.scripts}/subset_cds.py {input.ref} {output.ref_cds} &> {log}"

rule subset_prediction_cds:
    input:
        anno = os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff")
    output:
        anno_cds = os.path.join(dir.out.qc_gffcompare, f"{sample}.cds.gtf")
    log:
        os.path.join(dir.logs, "subset_prediction_cds.log")
    conda:
        os.path.join(dir.envs, "busco.yaml")
    shell:
        "python3 {dir.scripts}/subset_cds.py {input.anno} {output.anno_cds} &> {log}"

rule gffcompare_cds_eval:
    input:
        ref = os.path.join(dir.out.qc_gffcompare, "reference.cds.gtf"),
        anno = os.path.join(dir.out.qc_gffcompare, f"{sample}.cds.gtf")
    output:
        stats = os.path.join(dir.out.qc_gffcompare, f"{sample}_cds.stats")
    params:
        out_prefix = os.path.join(dir.out.qc_gffcompare, f"{sample}_cds")
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime = config.resources.small.time
    conda:
        os.path.join(dir.envs, "gffcompare.yaml")
    log:
        os.path.join(dir.logs, "gffcompare_cds.log")
    shell:
        "gffcompare -r {input.ref} -o {params.out_prefix} {input.anno} -T &> {log}"
    
