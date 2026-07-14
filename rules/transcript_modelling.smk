
rule run_isoquant:
    input:
        reads = config.project.input,
        ref = config.project.genome,
    output:
        gtf = os.path.join(dir.out.isoquant, f"{sample}.transcript_models.gtf")
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
        data_type = config.isoquant.data_type,
        prefix = sample
    log:
        os.path.join(dir.logs, "isoquant.log")
    shell:
        """
        isoquant.py \
            --reference {input.ref} \
            {params.input_flag} \
            --data_type {params.data_type} \
            --prefix {params.prefix} \
            --threads {threads} \
            -o {params.outdir} &> {log}
        """
