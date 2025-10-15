chromosomes=get_chromosomes(config.required.genome)

rule split_fasta:
    input:
        fasta = config.required.genome
    output:
        touch(os.path.join(dir.tools_reference,genome_name,f"{genome_name}_split.done"))
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    conda:
        f"{dir.envs}/sqanti3.yaml"
    threads:
        config.resources.small.cpus
    log:
        os.path.join(dir.logs,"split_fasta.log")
    script:
        os.path.join(dir.scripts,"splitfasta.py")

rule ed_augusuts_per_chromosome:
    input:
        os.path.join(dir.tools_reference,genome_name,f"{genome_name}_split.done"),
        mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done"),
        gff = os.path.join(dir.out.ed_hints,"IsoSeq.hints.gff")
    output:
        os.path.join(dir.out.ed_augustus,"split","{chromosome}.prediction.gff")
    conda:
        f"{dir.envs}/augustus.yaml"
    params:
        name = config.augustus.species_name,
        extcfg = f"{dir.envs}/extrinsic.M.RM.PB.cfg"
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    log:
        os.path.join(dir.logs,"ed_augustus_{chromosome}.log")
    threads:
        config.resources.small.cpus
    shell:
        """
        chromosome={dir.tools_reference}/{genome_name}/{wildcards.chromosome}.fasta
        augustus --species={params.name} $chromosome --hintsfile={input.gff} \
        --extrinsicCfgFile={params.extcfg} --protein=on --codingseq=on \
        --alternatives-from-evidence=true \
         > {output}  2>{log}
        """

rule merge_ed_predictions:
    input:
        expand(os.path.join(dir.out.ed_augustus,"split","{chromosome}.prediction.gff"),chromosome=chromosomes)
    output:
        temp(os.path.join(dir.out.ed_augustus,"Naive_prediction.gff"))
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    threads:
        config.resources.small.cpus
    shell:
        """
        for file in {input} ; do
            if [[ $(grep -v -c "#" $file) -gt 0 ]]; then
                cat $file >> {output}
            fi
        done
        """

rule rename_ed_augustus:
    input:
        os.path.join(dir.out.ed_augustus,"Naive_prediction.gff")
    output:
        os.path.join(dir.out.ed_augustus,"Augustus_prediction.gff")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "rename_augustus.log")
    script:
        f"{dir.scripts}/rename_augustus_genes.py"


rule ab_augustus_per_chromosome:
    input:
        os.path.join(dir.tools_reference,genome_name,f"{genome_name}_split.done"),
        mod = os.path.join(dir.out.ab_augustus_training,"SC_freq_mod.done")
    output:
        os.path.join(dir.out.ab_augustus,"split","{chromosome}.prediction.gff")
    conda:
        f"{dir.envs}/augustus.yaml"
    params:
        name = config.augustus.species_name,
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    log:
        os.path.join(dir.logs,"ab_augustus_{chromosome}.log")
    threads:
        config.resources.small.cpus
    shell:
        """
        chromosome={dir.tools_reference}/{genome_name}/{wildcards.chromosome}.fasta
        augustus --species={params.name} $chromosome --protein=on --codingseq=on > {output} 2>{log} 
        """

rule merge_ab_predictions:
    input:
        expand(os.path.join(dir.out.ab_augustus,"split","{chromosome}.prediction.gff"),chromosome=chromosomes)
    output:
        os.path.join(dir.out.ab_augustus,"split","ab_initio_prediction.gff")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    threads:
        config.resources.small.cpus
    shell:
        """
        cat {input} | grep -v "#" > {output}
        """

# Renaming predictions in the split mode
rule rename_ab_augustus:
    input:
        os.path.join(dir.out.ab_augustus,"split","ab_initio_prediction.gff")
    output:
        os.path.join(dir.out.ab_augustus,"ab_initio_prediction.gff")
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    log:
        os.path.join(dir.logs, "rename_augustus_ab.log")
    script:
        f"{dir.scripts}/rename_augustus_genes.py"

