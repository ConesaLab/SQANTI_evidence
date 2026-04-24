
if filetype == ".bam":
    rule fastq2bam:
        input:
            config.required.input
        output:
            os.path.join(dir.out.isoseq,f"{sample}.bam")
        conda:
            f"{dir.envs}/sqanti3.yaml"
        params:
            bam = os.path.join(dir.envs,"pacbio_mock.bam")
        threads:
            config.resources.small.cpus
        resources:
            slurm_extra = f"\'--qos={config.resources.small.qos}\'",
            cpus_per_task = config.resources.small.cpus,
            mem = config.resources.small.mem,
            runtime = config.resources.small.time
        script:
            os.path.join(dir.scripts,"fastq2bam.py")

rule index_genome:
    input:
        genome = config.required.genome
    output:
        os.path.join(dir.tools_index,genome_name,"index.mmi")
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads:
        config.resources.medium.cpus
    resources:
        slurm_extra = f"\'--qos={config.resources.medium.qos}\'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.big.mem,
        runtime = config.resources.medium.time
    shell:
        """
        mkdir -p {dir.tools_index}
        pbmm2 index {input.genome} {output}
        """

rule mapping_reads_pbmm2:
    input:
        reads = get_pbmm2_input(filetype,config,sample),
        index = os.path.join(dir.tools_index,genome_name,"index.mmi")
    output:
        os.path.join(dir.out.isoseq_mapping,f"{sample}.mapping_pbmm2.bam"),
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads:
        config.resources.big.cpus
    resources:
        cpus_per_task = config.resources.big.cpus,
        slurm_extra = f"\'--qos={config.resources.big.qos}\'",
        mem = config.resources.big.mem,
        runtime = config.resources.big.time
    log:
        os.path.join(dir.logs,"isoseq_mapping.log")
    shell:
        """
        pbmm2 align --preset ISOSEQ --sort {input.index} {input.reads}  {output} &> {log}
        """

# TODO: see how to use the FLNC BAM if we decide to use IsoSeq3
rule collapse_isoforms:
    input:
        mapped = os.path.join(dir.out.isoseq_mapping,f"{sample}.mapping_pbmm2.bam"),
        #flnc = os.path.join(config.required.flnc_dir,"{sample}","{sample}.flnc.bam")
    output:
        gff = os.path.join(dir.out.isoseq_collapsed,f"{sample}.collapsed.gff"),
        fastq = temp(os.path.join(dir.out.isoseq_collapsed,f"{sample}.collapsed.fastq"))
    conda:
        f"{dir.envs}/isoseq.yaml"
    log:
        os.path.join(dir.logs,"isoseq_collapse.log")
    threads:
        config.resources.small.cpus
    resources:
        slurm_extra = f"\'--qos={config.resources.small.qos}\'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.medium.mem,
        runtime = config.resources.small.time
    shell:
        """
        isoseq collapse --do-not-collapse-extra-5exons {input.mapped} {output.gff} -j {threads} &> {log}
        """
