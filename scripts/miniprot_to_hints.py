#!/usr/bin/env python3
"""
Convert Miniprot GFF spliced protein alignments into structured Augustus protein hints.

Features generated:
- CDSpart: Extracted from Miniprot CDS features with phase/frame (0, 1, 2) and src=P.
- intron: Inferred between consecutive CDS features with src=P.
- start/stop: Terminal boundary hints with src=P.

Author: Pablo Atienza & Gemini
"""
import sys
import os


def parse_miniprot_gff(miniprot_gff_path: str) -> list:
    """
    Parses a Miniprot GFF file and groups CDS segments by transcript/mRNA alignment.

    Returns:
        list of dict: Each entry contains {
            'id': str,
            'target': str,
            'chrom': str,
            'strand': str,
            'cds': list of (start, end, phase),
            'score': float or str,
            'identity': float
        }
    """
    alignments = {}
    current_mrna_id = None

    if not os.path.exists(miniprot_gff_path) or os.path.getsize(miniprot_gff_path) == 0:
        return []

    with open(miniprot_gff_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue

            chrom, source, ftype, start_s, end_s, score, strand, phase, attrs_s = fields
            start, end = int(start_s), int(end_s)

            # Parse attributes
            attrs = {}
            for item in attrs_s.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k.strip()] = v.strip()

            if ftype in ["mRNA", "transcript"]:
                mrna_id = attrs.get("ID", f"MP_{len(alignments)+1}")
                target = attrs.get("Target", mrna_id).split()[0]
                identity = float(attrs.get("Identity", 1.0))
                alignments[mrna_id] = {
                    "id": mrna_id,
                    "target": target,
                    "chrom": chrom,
                    "strand": strand,
                    "cds": [],
                    "score": score,
                    "identity": identity,
                }
                current_mrna_id = mrna_id

            elif ftype == "CDS":
                parent_id = attrs.get("Parent", current_mrna_id)
                if not parent_id:
                    parent_id = f"MP_unassigned_{len(alignments)+1}"
                    alignments[parent_id] = {
                        "id": parent_id,
                        "target": attrs.get("Target", parent_id).split()[0],
                        "chrom": chrom,
                        "strand": strand,
                        "cds": [],
                        "score": score,
                        "identity": 1.0,
                    }

                if parent_id not in alignments:
                    alignments[parent_id] = {
                        "id": parent_id,
                        "target": attrs.get("Target", parent_id).split()[0],
                        "chrom": chrom,
                        "strand": strand,
                        "cds": [],
                        "score": score,
                        "identity": 1.0,
                    }

                phase_val = phase if phase in ["0", "1", "2"] else "0"
                alignments[parent_id]["cds"].append((start, end, phase_val))

    return list(alignments.values())


def write_augustus_protein_hints(alignments: list, out_hints_path: str, src: str = "P", priority: int = 2):
    """
    Writes structured Augustus protein hints (CDSpart, intron, start, stop) from Miniprot alignments.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_hints_path)), exist_ok=True)
    count_cds = 0
    count_intron = 0
    count_start_stop = 0

    with open(out_hints_path, "w") as out_f:
        for aln in alignments:
            chrom = aln["chrom"]
            strand = aln["strand"]
            grp_id = aln["target"]
            cds_list = sorted(aln["cds"], key=lambda x: x[0])

            if not cds_list:
                continue

            hint_attr = f"grp={grp_id};pri={priority};src={src}"
            term_attr = f"grp={grp_id};pri={priority+1};src={src}"

            # 1. CDSpart hints
            for start, end, phase in cds_list:
                out_f.write(f"{chrom}\tminiprot\tCDSpart\t{start}\t{end}\t.\t{strand}\t{phase}\t{hint_attr}\n")
                count_cds += 1

            # 2. Intron hints
            for i in range(len(cds_list) - 1):
                intron_start = cds_list[i][1] + 1
                intron_end = cds_list[i + 1][0] - 1
                if intron_start <= intron_end:
                    out_f.write(f"{chrom}\tminiprot\tintron\t{intron_start}\t{intron_end}\t.\t{strand}\t.\t{hint_attr}\n")
                    count_intron += 1

            # 3. Start and Stop hints
            first_start = cds_list[0][0]
            last_end = cds_list[-1][1]
            if strand == "+":
                out_f.write(f"{chrom}\tminiprot\tstart\t{first_start}\t{first_start+2}\t.\t+\t0\t{term_attr}\n")
                out_f.write(f"{chrom}\tminiprot\tstop\t{last_end-2}\t{last_end}\t.\t+\t0\t{term_attr}\n")
                count_start_stop += 2
            elif strand == "-":
                out_f.write(f"{chrom}\tminiprot\tstart\t{last_end-2}\t{last_end}\t.\t-\t0\t{term_attr}\n")
                out_f.write(f"{chrom}\tminiprot\tstop\t{first_start}\t{first_start+2}\t.\t-\t0\t{term_attr}\n")
                count_start_stop += 2

    sys.stderr.write(f"### Converted {len(alignments)} Miniprot alignments into Augustus protein hints ({out_hints_path}):\n")
    sys.stderr.write(f"    - CDSpart hints: {count_cds}\n")
    sys.stderr.write(f"    - Intron hints:  {count_intron}\n")
    sys.stderr.write(f"    - Start/Stop:    {count_start_stop}\n")


def main():
    miniprot_gff_file = snakemake.input.miniprot_gff
    out_hints_file = snakemake.output.protein_hints

    src = getattr(snakemake.params, "src", "P")
    priority = int(getattr(snakemake.params, "priority", 2))

    sys.stderr.write(f"### Reading Miniprot GFF from: [{miniprot_gff_file}]\n")
    alignments = parse_miniprot_gff(miniprot_gff_file)

    sys.stderr.write(f"### Generating Augustus protein hints to: [{out_hints_file}]\n")
    write_augustus_protein_hints(alignments, out_hints_file, src=src, priority=priority)


if __name__ == "__main__":
    main()
