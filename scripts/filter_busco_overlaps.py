#!/usr/bin/env python3
"""
Script to filter BUSCO gene models against SQANTI dominant genes:
1. Direct gene locus overlap -> drop BUSCO gene.
2. Flanking safety check -> drop BUSCO gene if within flanking_size of a SQANTI gene.

Author: Pablo Atienza & Gemini
"""
import sys
import os
import re
from Bio import SeqIO


BUSCO_ID_RE = re.compile(r"^(\d+at\d+)")


def busco_gene_key(attrs: str):
    """
    Gene key for a BUSCO/Miniprot GFF record.

    BUSCO's miniprot GFF identifies alignments by an arbitrary `ID=MP######` / `Parent=MP######`,
    while the protein FASTA written by busco_complete_aa.py (and therefore the CD-HIT list) uses the
    BUSCO id, e.g. `10052at3699`. The BUSCO id is the prefix of the `Target=` attribute
    (`Target=10052at3699_29727_0:004f4a 1 586`). Keying BUSCO genes by it makes GFF, FAA and CD-HIT
    ids agree (roadmap 1.7b; before this, busco_non_overlapping.faa was always empty and no BUSCO gene
    ever reached the training set). Returns None if the record has no BUSCO-style Target.
    """
    for attr in attrs.split(";"):
        attr = attr.strip()
        if attr.startswith("Target="):
            m = BUSCO_ID_RE.match(attr[len("Target="):].strip())
            return m.group(1) if m else None
    return None


def resolve_gene_key(attrs: str, id_map: dict):
    """
    Gene key for one GFF record, remembering ID -> key so that child records that carry only
    `Parent=` (e.g. Miniprot `stop_codon` lines, which have no Target=) join their parent's gene.
    """
    key = busco_gene_key(attrs) or generic_gene_key(attrs)
    rec_id = parent = None
    for attr in attrs.split(";"):
        attr = attr.strip()
        if attr.startswith("ID="):
            rec_id = attr[3:].strip()
        elif attr.startswith("Parent="):
            parent = attr[7:].strip()
    if parent and parent in id_map:
        key = id_map[parent]
    if rec_id and key:
        id_map[rec_id] = key
    return key


def generic_gene_key(attrs: str):
    """gene_id "x" (GTF) or ID=/Parent= (GFF3, isoform suffix after '.' stripped); None if absent."""
    for attr in attrs.split(";"):
        attr = attr.strip()
        if attr.startswith("gene_id"):
            return attr.split("gene_id")[-1].strip().strip('"').strip("'").strip()
        if attr.startswith("ID="):
            return attr.split("ID=")[-1].strip().split(".")[0]
        if attr.startswith("Parent="):
            return attr.split("Parent=")[-1].strip().split(".")[0]
    return None


def extract_gene_intervals_from_gff(gff_path: str) -> dict:
    """
    Parse a GFF/GTF file and compute the bounding coordinates (min_start, max_end)
    for each gene locus, grouped by chromosome.

    Returns:
        dict: {chrom: [(start, end, gene_id), ...]}
    """
    genes = {}  # gene_id -> (chrom, min_start, max_end)
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attrs = fields[8]

            # Extract gene_id or ID
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
                gene_id = f"{chrom}_{start}_{end}"

            if gene_id not in genes:
                genes[gene_id] = [chrom, start, end]
            else:
                genes[gene_id][1] = min(genes[gene_id][1], start)
                genes[gene_id][2] = max(genes[gene_id][2], end)

    chrom_genes = {}
    for gid, (chrom, gstart, gend) in genes.items():
        if chrom not in chrom_genes:
            chrom_genes[chrom] = []
        chrom_genes[chrom].append((gstart, gend, gid))

    # Sort by start coordinate on each chromosome
    for chrom in chrom_genes:
        chrom_genes[chrom].sort(key=lambda x: x[0])

    return chrom_genes


def parse_busco_gff_records(busco_gff_path: str) -> tuple:
    """
    Parse BUSCO GFF file and group GFF lines by gene ID.

    Returns:
        tuple: (busco_genes_dict, busco_lines_by_gene)
    """
    busco_genes = {}  # gene_id -> [chrom, min_start, max_end]
    busco_lines = {}  # gene_id -> list of raw lines
    id_map = {}       # record ID -> gene key (lets Parent-only children follow their mRNA)

    with open(busco_gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attrs = fields[8]

            # BUSCO id from Target= (matches the FAA headers); fall back to ID=/Parent=/gene_id
            gene_id = resolve_gene_key(attrs, id_map)
            if not gene_id:
                first_attr = attrs.split(";")[0].strip()
                gene_id = first_attr.split("=")[-1].strip()

            if gene_id not in busco_genes:
                busco_genes[gene_id] = [chrom, start, end]
                busco_lines[gene_id] = []
            else:
                busco_genes[gene_id][1] = min(busco_genes[gene_id][1], start)
                busco_genes[gene_id][2] = max(busco_genes[gene_id][2], end)

            busco_lines[gene_id].append(line)

    return busco_genes, busco_lines


def filter_busco_genes(
    busco_genes: dict,
    sqanti_chrom_intervals: dict,
    flanking_size: int,
) -> tuple:
    """
    Filter BUSCO genes against SQANTI intervals for:
    1. Direct locus overlap
    2. Proximity within flanking_size
    """
    retained_ids = []
    dropped_overlap_count = 0
    dropped_flanking_count = 0

    for bg_id, (b_chrom, b_start, b_end) in busco_genes.items():
        if b_chrom not in sqanti_chrom_intervals:
            retained_ids.append(bg_id)
            continue

        sq_intervals = sqanti_chrom_intervals[b_chrom]
        has_direct_overlap = False
        too_close_flanking = False

        for s_start, s_end, _ in sq_intervals:
            # 1. Direct overlap check
            if b_start <= s_end and b_end >= s_start:
                has_direct_overlap = True
                break

            # 2. Distance check
            dist = max(0, max(s_start - b_end, b_start - s_end))
            if dist < flanking_size:
                too_close_flanking = True
                break

        if has_direct_overlap:
            dropped_overlap_count += 1
        elif too_close_flanking:
            dropped_flanking_count += 1
        else:
            retained_ids.append(bg_id)

    sys.stderr.write(f"### Initial BUSCO genes: {len(busco_genes)}\n")
    sys.stderr.write(f"### Dropped due to direct SQANTI overlap: {dropped_overlap_count}\n")
    sys.stderr.write(f"### Dropped due to flanking proximity (< {flanking_size} bp): {dropped_flanking_count}\n")
    sys.stderr.write(f"### Retained clean BUSCO genes: {len(retained_ids)}\n")

    return retained_ids


def filter_and_write_faa(busco_faa_path: str, retained_ids: set, out_faa_path: str) -> int:
    """Filter BUSCO FAA to write only retained gene sequences. Returns the number written."""
    written_count = 0
    with open(out_faa_path, "w") as out_f:
        for record in SeqIO.parse(busco_faa_path, "fasta"):
            # FAA headers are BUSCO ids (see busco_complete_aa.py); tolerate an isoform suffix
            rec_id = record.id.split(".")[0]
            if record.id in retained_ids or rec_id in retained_ids:
                SeqIO.write(record, out_f, "fasta")
                written_count += 1
    sys.stderr.write(f"### Wrote {written_count} sequences to {out_faa_path}\n")
    if retained_ids and written_count == 0:
        raise ValueError(
            f"No sequence of {busco_faa_path} matches any of the {len(retained_ids)} retained BUSCO gene ids "
            f"(e.g. {sorted(retained_ids)[:3]}). GFF and FAA identifiers disagree; the training set would silently "
            "lose all BUSCO genes."
        )
    return written_count


def main():
    busco_gff = snakemake.input.busco_gff
    busco_faa = snakemake.input.busco_faa
    sqanti_gff = snakemake.input.sqanti_gff
    out_clean_gff = snakemake.output.busco_clean_gff
    out_clean_faa = snakemake.output.busco_clean_faa
    flanking_size = int(snakemake.params.flanking_size)

    # Ensure output directories exist
    os.makedirs(os.path.dirname(out_clean_gff), exist_ok=True)
    os.makedirs(os.path.dirname(out_clean_faa), exist_ok=True)

    sys.stderr.write(f"### Loading SQANTI intervals from: [{sqanti_gff}]\n")
    sqanti_intervals = extract_gene_intervals_from_gff(sqanti_gff)

    sys.stderr.write(f"### Loading BUSCO GFF from: [{busco_gff}]\n")
    busco_genes, busco_lines = parse_busco_gff_records(busco_gff)

    retained_ids = filter_busco_genes(busco_genes, sqanti_intervals, flanking_size)
    retained_set = set(retained_ids)

    # Write clean BUSCO GFF
    sys.stderr.write(f"### Writing clean BUSCO GFF to: [{out_clean_gff}]\n")
    with open(out_clean_gff, "w") as out_gff:
        for gid in retained_ids:
            for line in busco_lines[gid]:
                out_gff.write(line)

    # Write clean BUSCO FAA
    sys.stderr.write(f"### Writing clean BUSCO FAA to: [{out_clean_faa}]\n")
    filter_and_write_faa(busco_faa, retained_set, out_clean_faa)


if __name__ == "__main__":
    main()
