#!/usr/bin/env python3
"""
Script to assemble final Augustus training GFF by combining non-redundant
BUSCO and SQANTI3 gene models based on CD-HIT clustering results.

Strategy: BUSCO Core + SQANTI Majority (Option 1)
- Guaranteed inclusion of all clean non-overlapping BUSCO representative genes.
- Fills remaining capacity up to max_genes with empirical SQANTI3 complete genes.

Author: Pablo Atienza & Gemini
"""
import sys
import os


def parse_cdhit_list(cdhit_lst_path: str) -> set:
    """Read representative gene/protein IDs from CD-HIT output list."""
    rep_ids = set()
    with open(cdhit_lst_path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                # First whitespace token is the sequence ID
                rep_ids.add(line.split()[0])
    return rep_ids


def parse_gff_by_gene(gff_path: str) -> tuple:
    """
    Parse a GFF/GTF file and group lines by gene ID.

    Returns:
        tuple: (gene_order_list, dict_of_gene_id_to_lines)
    """
    gene_order = []
    gene_lines = {}

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            attrs = fields[8]
            gene_id = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id"):
                    gene_id = attr.split("gene_id")[-1].strip().strip('"').strip("'").strip()
                    break
                elif attr.startswith("ID="):
                    gene_id = attr.split("ID=")[-1].strip().split(".")[0]
                    break
                elif attr.startswith("Parent="):
                    gene_id = attr.split("Parent=")[-1].strip().split(".")[0]
                    break

            if not gene_id:
                first_attr = attrs.split(";")[0].strip()
                gene_id = first_attr.split("=")[-1].strip()

            if gene_id not in gene_lines:
                gene_order.append(gene_id)
                gene_lines[gene_id] = []

            gene_lines[gene_id].append(line)

    return gene_order, gene_lines


def assemble_training_genes(
    cdhit_ids: set,
    busco_order: list,
    busco_dict: dict,
    sqanti_order: list,
    sqanti_dict: dict,
    max_genes: int,
    strategy: str = "busco_core",
) -> tuple:
    """
    Select training genes based on CD-HIT representatives and the chosen strategy.

    Strategies:
    - "busco_core" (Option 1): Include all BUSCO representatives first, fill remainder with SQANTI.
    - "sqanti_priority" (Option 3): Include all SQANTI representatives first, fill remainder with BUSCO.
    - "proportional": Split capacity evenly between available BUSCO and SQANTI models.

    Returns:
        tuple: (selected_busco_ids, selected_sqanti_ids)
    """
    # Filter for IDs present in CD-HIT representatives
    valid_busco_reps = [gid for gid in busco_order if gid in cdhit_ids or f"{gid}.t1" in cdhit_ids]
    valid_sqanti_reps = [gid for gid in sqanti_order if gid in cdhit_ids]

    sys.stderr.write(f"### Total CD-HIT representative IDs: {len(cdhit_ids)}\n")
    sys.stderr.write(f"### Available non-redundant BUSCO representatives: {len(valid_busco_reps)}\n")
    sys.stderr.write(f"### Available non-redundant SQANTI representatives: {len(valid_sqanti_reps)}\n")
    sys.stderr.write(f"### Target maximum training genes: {max_genes}\n")
    sys.stderr.write(f"### Selected Assembly Strategy: {strategy.upper()}\n")

    selected_busco = []
    selected_sqanti = []

    if strategy in ["busco_only"]:
        selected_busco = valid_busco_reps[:max_genes]

    elif strategy in ["sqanti_only"]:
        selected_sqanti = valid_sqanti_reps[:max_genes]

    elif strategy == "busco_core":
        # 1. Take all available BUSCO representatives
        selected_busco = valid_busco_reps[:max_genes]

        # 2. Fill remaining capacity with SQANTI representatives
        remaining_slots = max(0, max_genes - len(selected_busco))
        selected_sqanti = valid_sqanti_reps[:remaining_slots]

    elif strategy == "sqanti_priority":
        selected_sqanti = valid_sqanti_reps[:max_genes]
        remaining_slots = max(0, max_genes - len(selected_sqanti))
        selected_busco = valid_busco_reps[:remaining_slots]

    elif strategy == "proportional":
        half_cap = max_genes // 2
        selected_busco = valid_busco_reps[:half_cap]
        rem_for_sq = max_genes - len(selected_busco)
        selected_sqanti = valid_sqanti_reps[:rem_for_sq]

        # If SQANTI didn't use all its slots, give back to BUSCO
        if len(selected_sqanti) < rem_for_sq and len(valid_busco_reps) > len(selected_busco):
            extra_slots = max_genes - (len(selected_busco) + len(selected_sqanti))
            selected_busco.extend(valid_busco_reps[len(selected_busco) : len(selected_busco) + extra_slots])

    total_selected = len(selected_busco) + len(selected_sqanti)
    sys.stderr.write(f"### Assembled BUSCO core genes: {len(selected_busco)}\n")
    sys.stderr.write(f"### Assembled SQANTI empirical genes: {len(selected_sqanti)}\n")
    sys.stderr.write(f"### Final total training genes: {total_selected}\n")

    return selected_busco, selected_sqanti


def main():
    cdhit_lst_file = snakemake.input.cdhit_lst
    output_gff_file = snakemake.output.training_gff

    max_genes = int(getattr(snakemake.params, "max_genes", 5000))
    strategy = getattr(snakemake.params, "strategy", "busco_core")

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_gff_file), exist_ok=True)

    sys.stderr.write(f"### Reading CD-HIT list from: [{cdhit_lst_file}]\n")
    cdhit_ids = parse_cdhit_list(cdhit_lst_file)

    # Read BUSCO GFF if provided
    busco_order, busco_dict = [], {}
    busco_gff_file = getattr(snakemake.input, "busco_gff", None)
    if busco_gff_file and os.path.exists(busco_gff_file):
        sys.stderr.write(f"### Reading BUSCO GFF from: [{busco_gff_file}]\n")
        busco_order, busco_dict = parse_gff_by_gene(busco_gff_file)

    # Read SQANTI GFF if provided
    sqanti_order, sqanti_dict = [], {}
    sqanti_gff_file = getattr(snakemake.input, "sqanti_gff", None)
    if sqanti_gff_file and os.path.exists(sqanti_gff_file):
        sys.stderr.write(f"### Reading SQANTI GFF from: [{sqanti_gff_file}]\n")
        sqanti_order, sqanti_dict = parse_gff_by_gene(sqanti_gff_file)

    selected_busco, selected_sqanti = assemble_training_genes(
        cdhit_ids,
        busco_order,
        busco_dict,
        sqanti_order,
        sqanti_dict,
        max_genes,
        strategy=strategy,
    )

    sys.stderr.write(f"### Writing training GFF to: [{output_gff_file}]\n")
    with open(output_gff_file, "w") as out_f:
        # Write BUSCO records
        for gid in selected_busco:
            for line in busco_dict.get(gid, []):
                out_f.write(line)

        # Write SQANTI records
        for gid in selected_sqanti:
            for line in sqanti_dict.get(gid, []):
                out_f.write(line)


if __name__ == "__main__":
    main()
