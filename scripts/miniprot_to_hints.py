#!/usr/bin/env python3
"""
Convert Miniprot GFF spliced protein alignments into Augustus protein hints (src=P).

Hints generated per alignment:
- CDSpart : one per Miniprot CDS feature, with the GFF phase in column 8.
- intron  : inferred between consecutive CDS features of the same alignment.
- start   : only if the alignment begins at protein residue 1 (Target attribute), at the first codon.
- stop    : only if Miniprot emitted a `stop_codon` feature for the alignment, at that codon.

Design notes (see docs/roadmap_week_2026-09-07.md item 1.6):
- `grp=` is the Miniprot alignment ID (unique per alignment). Using the protein name grouped all
  secondary (paralog) alignments of one protein into a single Augustus hint group, so Augustus
  evaluated distinct loci as one gene. The protein name is kept as `target=` for traceability.
- start/stop hints were previously emitted for every alignment at its extremes with elevated priority,
  which forced wrong gene boundaries wherever a paralog alignment covered only part of the protein.
  Miniprot reports the aligned residue range and marks real stop codons, so those two facts gate them.
- Miniprot CDS features include the stop codon; the `stop_codon` feature coincides with the last codon.

Author: Pablo Atienza & Gemini; revised 2026-09-07.
"""
import os
import sys


def _parse_attrs(attrs_s):
    attrs = {}
    for item in attrs_s.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _parse_target(target_s):
    """'prot_name 12 345' -> ('prot_name', 12, 345); missing numbers -> None."""
    parts = target_s.split()
    name = parts[0] if parts else None
    try:
        return name, int(parts[1]), int(parts[2])
    except (IndexError, ValueError):
        return name, None, None


def parse_miniprot_gff(miniprot_gff_path: str) -> list:
    """
    Parses a Miniprot GFF file and groups CDS segments by alignment (mRNA).

    Returns a list of dicts:
        id        : Miniprot alignment ID (unique per alignment, e.g. MP000123)
        target    : protein name
        target_start, target_end : aligned protein residue range (None if absent)
        rank      : Miniprot Rank (1 = primary), int or None
        identity  : float
        score     : column 6 of the mRNA line
        chrom, strand
        cds       : list of (start, end, phase) sorted by start
        stop_codon: (start, end) of the Miniprot stop_codon feature, or None
    """
    alignments = {}

    if not os.path.exists(miniprot_gff_path) or os.path.getsize(miniprot_gff_path) == 0:
        sys.stderr.write(f"WARNING: Miniprot GFF missing or empty: {miniprot_gff_path}\n")
        return []

    def new_alignment(aid, attrs, chrom, strand, score):
        name, t_start, t_end = _parse_target(attrs.get("Target", aid))
        rank = attrs.get("Rank")
        return {
            "id": aid,
            "target": name or aid,
            "target_start": t_start,
            "target_end": t_end,
            "rank": int(rank) if rank and rank.isdigit() else None,
            "identity": float(attrs.get("Identity", 1.0)),
            "score": score,
            "chrom": chrom,
            "strand": strand,
            "cds": [],
            "stop_codon": None,
        }

    with open(miniprot_gff_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            chrom, _source, ftype, start_s, end_s, score, strand, phase, attrs_s = fields
            start, end = int(start_s), int(end_s)
            attrs = _parse_attrs(attrs_s)

            if ftype in ("mRNA", "transcript"):
                aid = attrs.get("ID", f"MP_{len(alignments) + 1}")
                alignments[aid] = new_alignment(aid, attrs, chrom, strand, score)
                continue

            parent = attrs.get("Parent")
            if not parent:
                continue
            if parent not in alignments:  # CDS before its mRNA line (defensive)
                alignments[parent] = new_alignment(parent, attrs, chrom, strand, score)

            if ftype == "CDS":
                phase_val = phase if phase in ("0", "1", "2") else "0"
                alignments[parent]["cds"].append((start, end, phase_val))
            elif ftype == "stop_codon":
                alignments[parent]["stop_codon"] = (start, end)

    for aln in alignments.values():
        aln["cds"].sort(key=lambda x: x[0])
    return list(alignments.values())


def alignment_hints(aln: dict, src: str = "P", priority: int = 2) -> list:
    """Returns the hint lines (without newline) for one alignment."""
    cds_list = aln["cds"]
    if not cds_list:
        return []
    chrom, strand = aln["chrom"], aln["strand"]
    attr = f"grp={aln['id']};pri={priority};src={src};target={aln['target']}"
    term_attr = f"grp={aln['id']};pri={priority + 1};src={src};target={aln['target']}"
    lines = []

    for start, end, phase in cds_list:
        lines.append(f"{chrom}\tminiprot\tCDSpart\t{start}\t{end}\t0\t{strand}\t{phase}\t{attr}")

    for i in range(len(cds_list) - 1):
        i_start, i_end = cds_list[i][1] + 1, cds_list[i + 1][0] - 1
        if i_start <= i_end:
            lines.append(f"{chrom}\tminiprot\tintron\t{i_start}\t{i_end}\t0\t{strand}\t.\t{attr}")

    # start: only when the protein N-terminus is aligned
    if aln.get("target_start") == 1:
        if strand == "+":
            s = cds_list[0][0]
            lines.append(f"{chrom}\tminiprot\tstart\t{s}\t{s + 2}\t0\t+\t0\t{term_attr}")
        else:
            e = cds_list[-1][1]
            lines.append(f"{chrom}\tminiprot\tstart\t{e - 2}\t{e}\t0\t-\t0\t{term_attr}")

    # stop: only when Miniprot found a stop codon
    if aln.get("stop_codon"):
        s, e = aln["stop_codon"]
        lines.append(f"{chrom}\tminiprot\tstop\t{s}\t{e}\t0\t{strand}\t0\t{term_attr}")

    return lines


def write_augustus_protein_hints(alignments: list, out_hints_path: str, src: str = "P", priority: int = 2):
    """Writes Augustus protein hints for all alignments and reports counts to stderr."""
    os.makedirs(os.path.dirname(os.path.abspath(out_hints_path)), exist_ok=True)
    counts = {"CDSpart": 0, "intron": 0, "start": 0, "stop": 0}
    n_written = 0
    with open(out_hints_path, "w") as out_f:
        for aln in alignments:
            lines = alignment_hints(aln, src=src, priority=priority)
            if not lines:
                continue
            n_written += 1
            for line in lines:
                counts[line.split("\t")[2]] += 1
                out_f.write(line + "\n")
    n_secondary = sum(1 for a in alignments if (a.get("rank") or 1) > 1)
    sys.stderr.write(f"### Converted {n_written} Miniprot alignments ({n_secondary} secondary) into Augustus protein hints ({out_hints_path}):\n")
    for k, v in counts.items():
        sys.stderr.write(f"    - {k:8s}: {v}\n")
    sys.stderr.write(f"    - alignments without start hint (partial N-terminus): {sum(1 for a in alignments if a['cds'] and a.get('target_start') != 1)}\n")
    sys.stderr.write(f"    - alignments without stop hint (no stop codon):      {sum(1 for a in alignments if a['cds'] and not a.get('stop_codon'))}\n")


def main():
    miniprot_gff_file = snakemake.input.miniprot_gff  # noqa: F821
    out_hints_file = snakemake.output.protein_hints  # noqa: F821
    src = getattr(snakemake.params, "src", "P")  # noqa: F821
    priority = int(getattr(snakemake.params, "priority", 2))  # noqa: F821

    sys.stderr.write(f"### Reading Miniprot GFF from: [{miniprot_gff_file}]\n")
    alignments = parse_miniprot_gff(miniprot_gff_file)
    sys.stderr.write(f"### Generating Augustus protein hints to: [{out_hints_file}]\n")
    write_augustus_protein_hints(alignments, out_hints_file, src=src, priority=priority)


if __name__ == "__main__":
    main()
