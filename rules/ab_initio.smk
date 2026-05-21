# Snakefile for ab initio gene prediction
import os


# Setup local rules (do not require much resources)
localrules:
    new_species,
    identify_bad_genes,
    extract_stop_codon_freq,


rule busco_run:
    input:
        genome=config.project.genome,
    output:
        dir=directory(dir.out.ab_busco),
    log:
        os.path.join(dir.logs, "busco_run.log"),
    conda:
        f"{dir.envs}/busco.yaml"
    threads: config.resources.busco.cpus
    resources:
        slurm_extra=f"'--qos={config.resources.busco.qos}'",
        cpus_per_task=config.resources.busco.cpus,
        mem=config.resources.busco.mem,
        runtime=config.resources.busco.time,
    params:
        busco_dir=dir.tools_busco,
        lineage=config.training.lineage,
        out_name=os.path.basename(dir.out.ab_busco),
        out_path=os.path.abspath(os.path.dirname(dir.out.ab_busco)),
    shell:
        """
        busco -i {input} -o {params.out_name} --out_path {params.out_path} \
            -l {params.lineage} -m genome --miniprot \
            -c {threads} --download_path {params.busco_dir} -f &> {log}
        """


rule busco_gather:
    input:
        dir.out.ab_busco,
    output:
        genes=os.path.join(dir.out.ab_augustus_model, "busco_genes.faa"),
    log:
        os.path.join(dir.logs, "busco_gather.log"),
    conda:
        os.path.join(dir.envs, "busco.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        lineage=config.training.lineage,
        gene_type="single",
    script:
        os.path.join(dir.scripts, "busco_complete_aa.py")


rule clustering_busco_genes:
    input:
        os.path.join(dir.out.ab_augustus_model, "busco_genes.faa"),
    output:
        os.path.join(dir.out.ab_augustus_model, "cdhit.lst"),
    log:
        os.path.join(dir.logs, "clustering_busco_genes.log"),
    conda:
        os.path.join(dir.envs, "busco.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    shell:
        """
        dir=$(dirname {output})
        cd-hit -o $dir/complete_buscos.cdhit -c 0.8 -i {input} -p 1 -d 0 -T 4 -M 48000 &> {log}
        grep ">" $dir/complete_buscos.cdhit | cut -f2 -d">" | cut -f1 > {output}
        """


rule concatenate_gff:
    input:
        gene_list=os.path.join(dir.out.ab_augustus_model, "cdhit.lst"),
        busco_path=dir.out.ab_busco,
    output:
        os.path.join(dir.out.ab_augustus_model, "busco_genes.gff"),
    log:
        os.path.join(dir.logs, "concatenate_gff.log"),
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        lineage=config.training.lineage,
        gene_type="single",
    run:
        with open(input.gene_list) as f:
            gene_names = [line.strip() for line in f if line.strip()]
        gff_path = os.path.join(
            input.busco_path,
            f"run_{params.lineage}",
            "busco_sequences",
            f"{params.gene_type}_copy_busco_sequences",
        )
        gff_files = [f"{gff_path}/{name}.gff" for name in gene_names]
        shell("cat {files} > {output}", files=" ".join(gff_files))


rule filter_miniprot_genes:
    input:
        os.path.join(dir.out.ab_augustus_model, "busco_genes.gff"),
    output:
        os.path.join(dir.out.ab_augustus_model, "busco_genes.filtered.gff"),
    log:
        os.path.join(dir.logs, "filter_miniprot_genes.log"),
    conda:
        os.path.join(dir.envs, "basic.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        threshold=config.training.miniprot_threshold,
    shell:
        """
        Rscript {dir.scripts}/filter_miniprot_genes.R {input} {output} {params.threshold} &> {log}
        """


rule gff2genbank:
    input:
        genome=config.project.genome,
        gff=os.path.join(dir.out.ab_augustus_model, "busco_genes.filtered.gff"),
    output:
        gen_bank=temp(os.path.join(dir.out.ab_augustus_model, "busco_genes.gb")),
    log:
        os.path.join(dir.logs, "gff2genbank.log"),
    conda:
        os.path.join(dir.envs, "busco.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        flanking_region=config.training.flanking_region,
    shell:
        """
        gff2gbSmallDNA.pl {input.gff} {input.genome} {params.flanking_region} {output} &> {log}
        """


rule generate_subsets:
    input:
        gen_bank_in=os.path.join(dir.out.ab_augustus_model, "busco_genes.gb"),
    output:
        gen_bank_out=os.path.join(dir.out.ab_augustus_model, "busco_genes.subset.gb"),
    log:
        os.path.join(dir.logs, "generate_subset.log"),
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        size=config.training.test_size,
        seed=123,
    script:
        os.path.join(dir.scripts, "generate_subset.py")


# TODO: Skip this rule if the directory of the new species exist
rule new_species:
    input:
        gen_bank=os.path.join(dir.out.ab_augustus_model, "busco_genes.subset.gb"),
    output:
        touch(
            os.path.join(
                dir.out.ab_augustus_model, f"{config.prediction.species}.done"
            )
        ),
    log:
        os.path.join(dir.logs, "new_species.log"),
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        name=config.prediction.species,
        augustus_dir=os.environ.get("AUGUSTUS_CONFIG_PATH"),
    shell:
        """
        rm -rf $AUGUSTUS_CONFIG_PATH/species/{params.name}
        new_species.pl --species={params.name} &> {log}
        """


rule initial_etraining:
    input:
        gb=os.path.join(dir.out.ab_augustus_model, "busco_genes.subset.gb"),
        new_species=os.path.join(
            dir.out.ab_augustus_model, f"{config.prediction.species}.done"
        ),
    output:
        training=os.path.join(dir.out.ab_augustus_training, "etrain.out"),
    log:
        os.path.join(dir.logs, "initial_etraining.log"),
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        name=config.prediction.species,
    shell:
        "etraining --species={params.name} {input.gb} &> {output}"


rule identify_bad_genes:
    input:
        training=os.path.join(dir.out.ab_augustus_training, "etrain.out"),
    output:
        bad=os.path.join(dir.out.ab_augustus_training, "bad.lst"),
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    shell:
        "grep 'in sequence' {input} | cut -f7 -d' ' | sed s/://g | sort -u > {output}"


rule filter_genes:
    input:
        bad_list=os.path.join(dir.out.ab_augustus_training, "bad.lst"),
        gb=os.path.join(dir.out.ab_augustus_model, "busco_genes.subset.gb"),
    output:
        filt=os.path.join(dir.out.ab_augustus_training, "filtered.gb"),
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    shell:
        "filterGenes.pl {input.bad_list} {input.gb} > {output}"


rule retrain:
    input:
        bad=os.path.join(dir.out.ab_augustus_training, "filtered.gb"),
    output:
        train=os.path.join(dir.out.ab_augustus_training, "etrain_filtered.out"),
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        name=config.prediction.species,
    shell:
        "etraining --species={params.name} {input} > {output}"


rule extract_stop_codon_freq:
    input:
        train=os.path.join(dir.out.ab_augustus_training, "etrain_filtered.out"),
    output:
        train=os.path.join(dir.out.ab_augustus_training, "SC_freq.txt"),
    shell:
        "tail -6 {input} | head -3 > {output}"


# TODO: Rewrite this rule to be adapted to snakemake nature. It creates a new file and then substitutes the frequency one.
# TODO: Find a way to give the shell variable AUGUSUTUS_CONFIG_PATH directly to the file
rule modify_stop_codon_freq:
    input:
        train=os.path.join(dir.out.ab_augustus_training, "SC_freq.txt"),
    output:
        mod=os.path.join(dir.out.ab_augustus_training, "SC_freq_mod.done"),
    log:
        os.path.join(dir.logs, "modify_stop_codon_freq.log"),
    conda:
        os.path.join(dir.envs, "augustus.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        name=config.prediction.species,
    script:
        os.path.join(dir.scripts, "modify_SC_freq.py")


# TODO: Is there any way to increase augustus usage to >1 core?
if config.prediction.mode == "split":

    include: "split_augustus.smk"

else:

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


rule placebo_gtf:
    input:
        config.project.genome,
    output:
        os.path.join(dir.out.ab_initio, "placebo.gtf"),
    log:
        os.path.join(dir.logs, "placebo_gtf.log"),
    conda:
        os.path.join(dir.envs, "sqanti3.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    script:
        os.path.join(dir.scripts, "generate_placebo_gtf.py")
