# Snakefile for ab initio gene prediction
import os

sp_name = config.prediction.species

# Setup local rules (do not require much resources)
localrules:
    new_species,
    identify_bad_genes,
    extract_stop_codon_freq,


rule select_dominant_sqanti_genes:
    input:
        gtf=os.path.join(dir.out.ed_sqanti, f"{sp_name}.filtered.regrouped.gtf"),
        classification=os.path.join(dir.out.ed_sqanti, f"{sp_name}_RulesFilter_classification.regrouped.txt"),
        genome=config.project.genome,
    output:
        gff=os.path.join(dir.out.ab_augustus_model, "sqanti_dominant.gff"),
        faa=os.path.join(dir.out.ab_augustus_model, "sqanti_dominant.faa"),
    log:
        os.path.join(dir.logs, "select_dominant_sqanti_genes.log"),
    conda:
        os.path.join(dir.envs, "basic.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    script:
        os.path.join(dir.scripts, "select_dominant_isoforms.py")


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


rule concatenate_gff:
    input:
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
    shell:
        """
        gff_path="{input.busco_path}/run_{params.lineage}/busco_sequences/{params.gene_type}_copy_busco_sequences"
        # sorted concatenation so the training set is identical between runs
        cat $(ls $gff_path/*.gff | sort) > {output} 2> {log}
        """


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


rule filter_busco_sqanti_overlaps:
    input:
        busco_gff=os.path.join(dir.out.ab_augustus_model, "busco_genes.filtered.gff"),
        busco_faa=os.path.join(dir.out.ab_augustus_model, "busco_genes.faa"),
        sqanti_gff=os.path.join(dir.out.ab_augustus_model, "sqanti_dominant.gff"),
    output:
        busco_clean_gff=os.path.join(dir.out.ab_augustus_model, "busco_non_overlapping.gff"),
        busco_clean_faa=os.path.join(dir.out.ab_augustus_model, "busco_non_overlapping.faa"),
    log:
        os.path.join(dir.logs, "filter_busco_sqanti_overlaps.log"),
    conda:
        os.path.join(dir.envs, "basic.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        flanking_size=config.training.flanking_region,
    script:
        os.path.join(dir.scripts, "filter_busco_overlaps.py")


def get_cluster_faa_inputs(wildcards):
    mode = config.training.mode
    if mode == "sqanti_only":
        return [os.path.join(dir.out.ab_augustus_model, "sqanti_dominant.faa")]
    elif mode == "busco_only":
        return [os.path.join(dir.out.ab_augustus_model, "busco_genes.faa")]
    else:  # mixed
        return [
            os.path.join(dir.out.ab_augustus_model, "sqanti_dominant.faa"),
            os.path.join(dir.out.ab_augustus_model, "busco_non_overlapping.faa"),
        ]


def get_assemble_gff_inputs(wildcards):
    mode = config.training.mode
    inputs = {
        "cdhit_lst": os.path.join(dir.out.ab_augustus_model, "cdhit.lst"),
    }
    if mode in ["mixed", "sqanti_only"]:
        inputs["sqanti_gff"] = os.path.join(dir.out.ab_augustus_model, "sqanti_dominant.gff")
    if mode in ["mixed", "busco_only"]:
        inputs["busco_gff"] = os.path.join(
            dir.out.ab_augustus_model,
            "busco_non_overlapping.gff" if mode == "mixed" else "busco_genes.filtered.gff",
        )
    return inputs


rule cluster_training_proteins:
    input:
        get_cluster_faa_inputs,
    output:
        combined_faa=os.path.join(dir.out.ab_augustus_model, "training_proteins.faa"),
        cdhit_lst=os.path.join(dir.out.ab_augustus_model, "cdhit.lst"),
    log:
        os.path.join(dir.logs, "cluster_training_proteins.log"),
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
        cat {input} > {output.combined_faa}
        dir=$(dirname {output.cdhit_lst})
        cd-hit -i {output.combined_faa} -o $dir/training_proteins.cdhit \
               -c 0.8 -p 1 -d 0 -T {threads} -M 48000 &> {log}
        grep ">" $dir/training_proteins.cdhit | cut -f2 -d">" | cut -f1 > {output.cdhit_lst}
        """


rule assemble_training_gff:
    input:
        unpack(get_assemble_gff_inputs),
    output:
        training_gff=os.path.join(dir.out.ab_augustus_model, "training_genes.gff"),
    log:
        os.path.join(dir.logs, "assemble_training_gff.log"),
    conda:
        os.path.join(dir.envs, "basic.yaml")
    resources:
        slurm_extra=f"'--qos={config.resources.small.qos}'",
        cpus_per_task=config.resources.small.cpus,
        mem=config.resources.small.mem,
        runtime=config.resources.small.time,
    params:
        max_genes=config.training.test_size,
        strategy=lambda wildcards: "busco_core" if config.training.mode == "mixed" else config.training.mode,
    script:
        os.path.join(dir.scripts, "assemble_training_gff.py")


rule gff2genbank:
    input:
        genome=config.project.genome,
        gff=os.path.join(dir.out.ab_augustus_model, "training_genes.gff"),
    output:
        gen_bank=temp(os.path.join(dir.out.ab_augustus_model, "training_genes.gb")),
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
        # PERL_HASH_SEED=0: Perl randomises hash iteration order per process, which changed which of
        # two overlapping/flanking-adjacent loci gff2gbSmallDNA.pl kept between runs (same input, 9/4,696
        # loci swapped on Arabidopsis) and thus the trained model. Fixing the seed makes it reproducible.
        export PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0
        gff2gbSmallDNA.pl {input.gff} {input.genome} {params.flanking_region} {output} &> {log}
        """


# Creates the Augustus species from the generic template. The species directory lives in
# $AUGUSTUS_CONFIG_PATH/species/<name> (set by the Augustus conda env, i.e. inside toolsdir), not in
# outdir. Any previous copy is removed so retraining always starts from the template; the .done
# sentinel keeps Snakemake from redoing this unless the training GenBank changed.
rule new_species:
    input:
        gen_bank=os.path.join(dir.out.ab_augustus_model, "training_genes.gb"),
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
    shell:
        """
        if [ -z "$AUGUSTUS_CONFIG_PATH" ] || [ ! -w "$AUGUSTUS_CONFIG_PATH/species" ]; then
            echo "ERROR: AUGUSTUS_CONFIG_PATH/species ('${{AUGUSTUS_CONFIG_PATH:-unset}}/species') is not writable;" \
                 "cannot create the Augustus species '{params.name}'." | tee {log} >&2
            exit 1
        fi
        echo "Creating Augustus species '{params.name}' in $AUGUSTUS_CONFIG_PATH/species/{params.name}" > {log}
        rm -rf "$AUGUSTUS_CONFIG_PATH/species/{params.name}"
        export PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0
        new_species.pl --species={params.name} &>> {log}
        """


rule initial_etraining:
    input:
        gb=os.path.join(dir.out.ab_augustus_model, "training_genes.gb"),
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
        gb=os.path.join(dir.out.ab_augustus_model, "training_genes.gb"),
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
        "export PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0; filterGenes.pl {input.bad_list} {input.gb} > {output}"


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


# Writes the observed stop-codon frequencies into the species parameters file (in
# $AUGUSTUS_CONFIG_PATH, outside outdir) and keeps a copy of the modified file in the training
# directory so the run is inspectable. The .done sentinel is what downstream rules depend on.
rule modify_stop_codon_freq:
    input:
        train=os.path.join(dir.out.ab_augustus_training, "SC_freq.txt"),
    output:
        mod=os.path.join(dir.out.ab_augustus_training, "SC_freq_mod.done"),
        params_cfg=os.path.join(dir.out.ab_augustus_training, f"{config.prediction.species}_parameters.cfg"),
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
