
rule generate_proteome:
    input:
        annot = os.path.join(dir.out.ed_augustus,"Augustus_prediction.gff"), 
        reference = config.required.genome,
    output:
        os.path.join(dir.out.ed_augustus,"Augustus_prediction.aa")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "generate_proteome.log")
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    shell:
        """
        getAnnoFasta.pl {input.annot} --seqfile={input.reference}
        """

# OMARk results

rule omamer:
    input:
        proteome = os.path.join(dir.out.ed_augustus,"Augustus_prediction.aa"),
        omark_db = os.path.join(dir.tools_omark,f"{config.qc.omark_db}.h5")
    output:
        os.path.join(dir.out.qc_omark,"Final.omamer")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.medium.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "omamer.log")
    conda:
        os.path.join(dir.envs, "omark.yaml")
    shell:
        """
        omamer search --db {input.omark_db} --query {input.proteome} --out {output} &> {log}
        """

rule omark:
    input:
        omamer = os.path.join(dir.out.qc_omark,"Final.omamer"),
        omark_db = os.path.join(dir.tools_omark,f"{config.qc.omark_db}.h5")
    output:
        os.path.join(dir.out.qc_omark,"Final.pdf")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "omark.log")
    conda:
        os.path.join(dir.envs, "omark.yaml")
    threads:
        config.resources.small.cpus
    shell:
        """
        omark -f {input.omamer} -d {input.omark_db} -o $(dirname {output}) &> {log}
        """

rule busco_qc:
    input:
        proteome = os.path.join(dir.out.ed_augustus,"Augustus_prediction.aa"),        
    output:
        directory(os.path.join(dir.out.qc_busco)),
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.busco.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.medium.time
    params:
        busco_dir = dir.tools_busco,
        lineage = config.ab_initio.lineage,
    log: 
        os.path.join(dir.logs, "busco_qc.log")
    conda:
        os.path.join(dir.envs, "busco.yaml")
    threads:
        config.resources.busco.cpus
    shell:
        """
        busco -i {input.proteome} -o {output} -l {params.lineage} \
            -m proteins -c {threads} --force --download_path {params.busco_dir} &> {log}
        """

rule agat_cleaning:
    input:
        os.path.join(dir.out.ed_augustus,"Augustus_prediction.gff")
    output:
        os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.medium.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "agat_cleaning.log")
    conda:
        os.path.join(dir.envs, "agat.yaml")
    threads:
        config.resources.small.cpus
    shell:
        """
        agat_convert_sp_gxf2gxf.pl -g {input} -o {output} &> {log}
        """

rule agat_stats:
    input:
        os.path.join(dir.out.evidence_driven,"Final_clean_prediction.gff")
    output:
        os.path.join(dir.out.qc_agat,"Final_stats.txt")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "agat_stats.log")
    conda:
        os.path.join(dir.envs, "agat.yaml")
    threads:
        config.resources.small.cpus
    shell:
        """
        agat_sp_statistics.pl --gff {input} -o {output} &> {log}
        """

rule gaqet2_setup:
    input:
        os.path.join(dir.tools_gaqet2,"gaqet2_installed.done")
    output:
        touch(os.path.join(dir.tools_gaqet2,"gaqet2_setup.done"))
    conda:
        os.path.join(dir.envs, "gaqet2.yaml")
    params:
        config = os.path.join(dir.tools_gaqet2, "gaqet2_config.yaml"),
        lineage = config.busco.lineage,
    log:
        os.path.join(dir.logs, "gaqet2_setup.log")
    script:
        f"{dir.scripts}/gaqet2_setup.py"

    
#rule gaqet2:
#    input:
#        touch(os.path.join(dir.tools_gaqet2,"gaqet2_installed.done")),
#        config = os.path.join(dir.envs, "gaqet2_config.yaml"),
#    output:
#    
#    conda:
#        os.path.join(dir.envs, "gaqet2.yaml")
