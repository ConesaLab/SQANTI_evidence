
import yaml

def main():
    with open(snakemake.params.config, "r") as f_in:
        base_config = yaml.safe_load(f_in)

    base_config["ID"] = snakemake.params.id
    base_config["Assembly"] = snakemake.input.genome
    base_config["Annotation"] = snakemake.input.annotation
    base_config["Basedir"] = snakemake.params.outdir
    base_config["Threads"] = snakemake.params.threads
    base_config["Analysis"] = ["AGAT", "BUSCO", "OMARK"]
    base_config["OMARK_db"] = snakemake.params.omark_db
    base_config["OMARK_taxid"] = snakemake.params.taxid
    base_config["BUSCO_lineages"] = [snakemake.params.busco_lineage]

    with open(snakemake.output[0], "w") as f_out:
        yaml_dump = yaml.dump(base_config, default_flow_style=False)
        f_out.write(yaml_dump)

if __name__ == "__main__":
    main()
