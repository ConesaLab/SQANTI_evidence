import glob
import attrmap as ap

config = ap.AttrMap(config)

# Functions

def get_samples(filename):
    """
    Get samples from the config file.
    """
    with open(filename, 'r') as f:
        files = [line.strip() for line in f]
        samples = [os.path.basename(file).split('.')[0] for file in files]
    return samples, files


samples, files = get_samples_and_files("samples.txt")

rule lima:
    input:
        lambda wildcards: files[samples.index(wildcards.sample)]
    output:
        os.path.join(dir.out.isoseq_lima,"{sample}","fl.lima.bam")
    conda:
        f"{dir.envs}/isoseq.yaml"
    params:
        output = os.path.join(dir.out.isoseq_lima,"{sample}","fl.bam"),
        primers = config.isoseq.primers,
        samples = config.isoseq.primers_to_samples
    log:
        os.path.join(dir.logs,"lima_demultiplexing_{sample}.log")
    threads:
        config.resources.big.get("cpus", 4)
    resources:
        slurm_extra = f"'--qos={config.resources.big.qos}'",
        cpus_per_task = config.resources.big.cpus,
        mem = config.resources.big.mem,
        runtime =  config.resources.big.time
    shell:
        """
        lima {input} {params.primers} {params.output} \
            --isoseq --peek-guess --overwrite-biosample-names \
            --split-named --biosample-csv {params.samples} --log-level INFO --log-file {log}
        """

rule refine:
    input:
        lima = os.path.join(dir.out.isoseq_lima,"{sample}","fl.lima.bam")
    output:
        os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
    conda:
        f"{dir.envs}/isoseq.yaml"
    log:
        os.path.join(dir.logs,"isoseq_refine_{sample}.log")
    params:
        primers = config.isoseq.primers
    threads:
        config.resources.small.get("cpus", 4)
    resources:
        slurm_extra = f"'--qos={config.resources.small.qos}'",
        cpus_per_task = config.resources.small.cpus,
        mem = config.resources.small.mem,
        runtime =  config.resources.small.time
    shell:
        """
        isoseq refine --require-polya {input.lima} {params.primers} {output} &> {log}
        """

rule create_new_fofn:
    input:
        expand(os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam"),sample=samples)
    output:
        "results/isoannot/input.fofn"
    conda:
        f"{dir.envs}/isoseq.yaml"
    shell:
        """
        ls {input} > {output}
        """

rule cluster:
    input:
        os.path.join(dir.out.isoseq_refine,"{sample}","{sample}.flnc.bam")
    output:
        bam=os.path.join(dir.out.isoseq_cluster,"{sample}.cluster.bam"),
    conda:
        f"{dir.envs}/isoseq.yaml"
    log:
        os.path.join(dir.logs,"isoseq_cluster_{sample}.log")
    threads:
        config.resources.small.get("cpus", 4)
    resources:
        slurm_extra = f"'--qos={config.resources.medium.qos}'",
        cpus_per_task = config.resources.medium.cpus,
        mem = config.resources.small_bigMem.mem,
        runtime =  config.resources.medium.time
    shell:
        """
        isoseq cluster2 {input} {output.bam} &> {log}
        """
