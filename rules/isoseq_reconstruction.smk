# IsoSeq3 transcriptome reconstruction: (fastq2bam) -> cluster -> align -> collapse.
#
# Alternative producer of the transcriptome contract consumed by run_sqanti and restore_gene_ids:
# a GTF with gene_id/transcript_id attributes plus a two-column FL-count matrix (see
# get_transcriptome() in rules/setup/functions.smk). Included only when
# curation.reconstruction == "isoseq"; everything downstream of the transcriptome is unchanged.
#
# Commands follow the pre-refactor dev-IsoSeq branch (rules/isoseq.smk and
# rules/transcript_modelling.smk at dc36bb7^), rewritten against the current configuration schema.
# Design notes: docs/isoseq_reconstruction_branch_design.md.


# `isoseq cluster2` reads PacBio BAM only, so a FASTA/FASTQ input is wrapped in PacBio records using a
# template. The run metadata of the result is synthetic; input_check.py warns about it at configuration
# time and the script repeats the warning in its log.
if filetype != ".bam":

    rule fastq2bam:
        input:
            config.project.input,
        output:
            os.path.join(dir.out.isoseq, f"{sample}.bam"),
        log:
            os.path.join(dir.logs, "fastq2bam.log"),
        conda:
            f"{dir.envs}/isoseq.yaml"
        threads: config.resources.small.cpus
        params:
            template=os.path.join(dir.envs, "pacbio_mock.bam"),
        resources:
            slurm_extra=f"'--qos={config.resources.small.qos}'",
            cpus_per_task=config.resources.small.cpus,
            mem=config.resources.medium.mem,
            runtime=config.resources.medium.time,
        script:
            os.path.join(dir.scripts, "fastq2bam.py")


# Clustering is a mandatory step of the PacBio Iso-Seq workflow, not a tuning choice: it merges the
# full-length reads into non-redundant transcript sequences before they are aligned.
rule isoseq_cluster:
    input:
        flnc=get_isoseq_reads(filetype, config, sample),
    output:
        bam=os.path.join(dir.out.isoseq_cluster, f"{sample}.cluster.bam"),
    log:
        os.path.join(dir.logs, "isoseq_cluster.log"),
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads: config.resources.big.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.big.qos}'",
        cpus_per_task=config.resources.big.cpus,
        mem=config.resources.big.mem,
        runtime=config.resources.big.time,
    shell:
        """
        isoseq cluster2 {input.flnc} {output.bam} -j {threads} &> {log}
        """


rule pbmm2_index:
    input:
        genome=config.project.genome,
    output:
        os.path.join(dir.tools_pbmm2, genome_name, "isoseq_index.mmi"),
    log:
        os.path.join(dir.logs, "pbmm2_index.log"),
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads: config.resources.medium.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.medium.qos}'",
        cpus_per_task=config.resources.medium.cpus,
        mem=config.resources.big.mem,
        runtime=config.resources.medium.time,
    shell:
        """
        mkdir -p $(dirname {output})
        pbmm2 index {input.genome} {output} &> {log}
        """


rule pbmm2_align:
    input:
        reads=os.path.join(dir.out.isoseq_cluster, f"{sample}.cluster.bam"),
        index=os.path.join(dir.tools_pbmm2, genome_name, "isoseq_index.mmi"),
    output:
        os.path.join(dir.out.isoseq_mapping, f"{sample}.aligned.bam"),
    log:
        os.path.join(dir.logs, "pbmm2_align.log"),
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads: config.resources.big.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.big.qos}'",
        cpus_per_task=config.resources.big.cpus,
        mem=config.resources.big.mem,
        runtime=config.resources.big.time,
    shell:
        """
        pbmm2 align --preset ISOSEQ --sort -j {threads} {input.index} {input.reads} {output} &> {log}
        """


# --do-not-collapse-extra-5exons keeps 5'-degraded variants as separate models. That permissiveness is
# the property under comparison against IsoQuant, so the flag is always set and is rejected if a user
# repeats it through curation.isoseq_args.
rule isoseq_collapse:
    input:
        mapped=os.path.join(dir.out.isoseq_mapping, f"{sample}.aligned.bam"),
    output:
        gff=os.path.join(dir.out.isoseq_collapsed, f"{sample}.collapsed.gff"),
        abundance=os.path.join(dir.out.isoseq_collapsed, f"{sample}.collapsed.abundance.txt"),
    log:
        os.path.join(dir.logs, "isoseq_collapse.log"),
    conda:
        f"{dir.envs}/isoseq.yaml"
    threads: config.resources.medium.cpus
    params:
        extra=config.curation.isoseq_args,
    resources:
        slurm_extra=f"'--qos={config.resources.medium.qos}'",
        cpus_per_task=config.resources.medium.cpus,
        mem=config.resources.big.mem,
        runtime=config.resources.medium.time,
    shell:
        """
        isoseq collapse --do-not-collapse-extra-5exons {params.extra} \
            {input.mapped} {output.gff} -j {threads} &> {log}
        """


rule collapse_fl_counts:
    input:
        abundance=os.path.join(dir.out.isoseq_collapsed, f"{sample}.collapsed.abundance.txt"),
    output:
        os.path.join(dir.out.isoseq_collapsed, f"{sample}.fl_counts.tsv"),
    log:
        os.path.join(dir.logs, "collapse_fl_counts.log"),
    conda:
        f"{dir.envs}/basic.yaml"
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    script:
        os.path.join(dir.scripts, "collapse_counts_to_fl.py")
