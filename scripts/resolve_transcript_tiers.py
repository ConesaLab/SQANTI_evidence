#!/usr/bin/env python3
"""
scripts/resolve_transcript_tiers.py
Snakemake-native script for Tier 1 (SQANTI3) Pass-Through + Tier 2 (Augustus) Gap-Filling.

Logic:
  - Tier 1: All SQANTI3-curated models pass directly to output verbatim, preserving full
    transcript structures, UTRs, and empirical splice junctions.
  - Tier 2: Augustus predictions are evaluated against Tier 1 genomic loci.
    - If an Augustus model overlaps an existing Tier 1 gene locus on the same strand -> DISCARD
      (prevents redundant dual-transcripts and precision collapse).
    - If non-overlapping multi-exon -> RETAIN as gap-filler for unexpressed genes.
    - If non-overlapping single-exon -> Retain if length >= min_monoexon_len or hint-supported.
"""

import os
import sys
import re
import bisect
from collections import defaultdict


def extract_gene_id_from_attrs(attrs_str):
    """Extract gene_id from GTF/GFF attributes column."""
    if not attrs_str or attrs_str == ".":
        return None
    # Check standard GTF: gene_id "XYZ"
    m = re.search(r'gene_id\s+"([^"]+)"', attrs_str)
    if m:
        return m.group(1)
    # Check GFF3: ID=XYZ or Parent=XYZ
    m = re.search(r'(?:ID|gene_id)=([^;]+)', attrs_str)
    if m:
        return m.group(1)
    return None


def parse_sqanti_gtf(filepath):
    """
    Parses SQANTI3 GTF file.
    Returns:
        genes (dict): gene_id -> list of raw lines
        intervals (dict): (chrom, strand) -> list of (start, end, gene_id)
    """
    genes = defaultdict(list)
    gene_spans = {}  # gene_id -> [chrom, strand, min_start, max_end]

    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom, feature, start, end, strand, attrs = (
                parts[0],
                parts[2],
                int(parts[3]),
                int(parts[4]),
                parts[6],
                parts[8],
            )

            gid = extract_gene_id_from_attrs(attrs)
            if not gid:
                # Fallback: check transcript_id
                m = re.search(r'transcript_id\s+"([^"]+)"', attrs)
                gid = m.group(1) if m else f"{chrom}_{start}_{end}"

            genes[gid].append(line)

            if gid not in gene_spans:
                gene_spans[gid] = [chrom, strand, start, end]
            else:
                gene_spans[gid][2] = min(gene_spans[gid][2], start)
                gene_spans[gid][3] = max(gene_spans[gid][3], end)

    intervals = defaultdict(list)
    for gid, (chrom, strand, start, end) in gene_spans.items():
        intervals[(chrom, strand)].append((start, end, gid))

    for key in intervals:
        intervals[key].sort(key=lambda x: (x[0], x[1]))

    return genes, intervals


def parse_augustus_gff(filepath):
    """
    Parses Augustus GFF/GTF predictions.
    Returns:
        genes (dict): gene_id -> list of raw lines
        gene_info (dict): gene_id -> {
            'chrom': chrom,
            'strand': strand,
            'start': start,
            'end': end,
            'exons': list of (start, end),
            'has_intron': bool
        }
    """
    genes = defaultdict(list)
    gene_info = {}
    current_gene_id = None

    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                # Augustus comments: # start gene g1
                m = re.search(r"#\s*start\s+gene\s+(\S+)", line)
                if m:
                    current_gene_id = m.group(1)
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom, feature, start, end, strand, attrs = (
                parts[0],
                parts[2],
                int(parts[3]),
                int(parts[4]),
                parts[6],
                parts[8],
            )

            gid = None
            if feature == "gene":
                gid = attrs.strip()
                current_gene_id = gid
            else:
                gid = extract_gene_id_from_attrs(attrs)
                if not gid:
                    gid = current_gene_id
                if not gid:
                    m = re.search(r'transcript_id\s+"([^"]+)"', attrs)
                    if m:
                        gid = m.group(1).split(".t")[0]
                    else:
                        gid = f"{chrom}_{start}_{end}"

            genes[gid].append(line)

            if gid not in gene_info:
                gene_info[gid] = {
                    "chrom": chrom,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "exons": [],
                    "cds": [],
                    "has_intron": False,
                }
            else:
                gene_info[gid]["start"] = min(gene_info[gid]["start"], start)
                gene_info[gid]["end"] = max(gene_info[gid]["end"], end)

            if feature == "exon":
                gene_info[gid]["exons"].append((start, end))
            elif feature == "CDS":
                gene_info[gid]["cds"].append((start, end))
            elif feature == "intron":
                gene_info[gid]["has_intron"] = True

    # Post-process exon counts / introns
    for gid, info in gene_info.items():
        unique_exons = set(info["exons"])
        if not unique_exons:
            unique_exons = set(info["cds"])
        if len(unique_exons) > 1 or info.get("has_intron", False):
            info["has_intron"] = True
        else:
            info["has_intron"] = False

    return genes, gene_info


def load_hints_coords(hints_file):
    """Loads hint start/stop/CDS coordinates for monoexon validation."""
    hint_intervals = defaultdict(list)
    if not hints_file or not os.path.isfile(hints_file):
        return hint_intervals

    with open(hints_file, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, feat, start, end, strand = parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6]
            if feat in ["start", "stop", "start_codon", "stop_codon", "CDS", "CDSpart"]:
                hint_intervals[(chrom, strand)].append((start, end))

    for key in hint_intervals:
        hint_intervals[key].sort()

    return hint_intervals


def check_hint_support(chrom, strand, start, end, hint_intervals):
    """Checks if interval overlaps any hint."""
    hints = hint_intervals.get((chrom, strand), [])
    if not hints:
        return False
    # Binary search for start
    idx = bisect.bisect_right([h[1] for h in hints], start)
    while idx < len(hints):
        h_start, h_end = hints[idx]
        if h_start > end:
            break
        if not (end < h_start or start > h_end):
            return True
        idx += 1
    return False


def check_locus_overlap(chrom, strand, start, end, intervals_dict):
    """
    Checks if genomic interval [start, end] overlaps any interval in intervals_dict[(chrom, strand)].
    """
    key = (chrom, strand)
    if key not in intervals_dict:
        return False

    intervals = intervals_dict[key]
    ends = [iv[1] for iv in intervals]
    idx = bisect.bisect_right(ends, start)

    while idx < len(intervals):
        s_start, s_end, _ = intervals[idx]
        if s_start > end:
            break
        if not (end < s_start or start > s_end):
            return True
        idx += 1

    return False


def resolve_tiers(
    sqanti_gtf,
    augustus_gff,
    output_gtf,
    hints_file=None,
    min_monoexon_len=300,
    filter_mode="medium",
):
    """
    Core resolution function:
      1. Writes 100% of SQANTI3 curated models to output_gtf.
      2. For each Augustus model:
         - Discards if overlapping any SQANTI3 locus on the same strand.
         - For non-overlapping multi-exon: retains.
         - For non-overlapping single-exon: checks length / hint support.
    """
    print(f"Loading Tier 1 SQANTI3 models from: {sqanti_gtf}")
    sqanti_genes, sqanti_intervals = parse_sqanti_gtf(sqanti_gtf)
    print(f"  Loaded {len(sqanti_genes)} SQANTI3 gene loci.")

    print(f"Loading Tier 2 Augustus predictions from: {augustus_gff}")
    aug_genes, aug_info = parse_augustus_gff(augustus_gff)
    print(f"  Loaded {len(aug_info)} Augustus gene predictions.")

    hint_intervals = load_hints_coords(hints_file) if hints_file else {}

    tier1_count = 0
    tier2_gap_count = 0
    dropped_overlap_count = 0
    dropped_monoexon_noise = 0

    os.makedirs(os.path.dirname(os.path.abspath(output_gtf)), exist_ok=True)

    with open(output_gtf, "w") as out:
        # Step 1: Write all Tier 1 SQANTI3 models
        for gid, lines in sqanti_genes.items():
            for line in lines:
                out.write(line)
            tier1_count += 1

        # Step 2: Evaluate Augustus predictions
        for gid, info in aug_info.items():
            chrom = info["chrom"]
            strand = info["strand"]
            start = info["start"]
            end = info["end"]
            has_intron = info["has_intron"]
            span_len = end - start + 1

            # Locus overlap check against Tier 1
            if check_locus_overlap(chrom, strand, start, end, sqanti_intervals):
                dropped_overlap_count += 1
                continue

            # Noise filter for non-overlapping single-exon predictions
            if not has_intron:
                if filter_mode == "strict":
                    # Must have hint support
                    if not check_hint_support(chrom, strand, start, end, hint_intervals):
                        dropped_monoexon_noise += 1
                        continue
                elif filter_mode == "medium":
                    # Must have minimum length OR hint support
                    supported = check_hint_support(chrom, strand, start, end, hint_intervals)
                    if span_len < min_monoexon_len and not supported:
                        dropped_monoexon_noise += 1
                        continue
                elif filter_mode == "none":
                    pass

            # Retain non-overlapping Augustus gap-filler
            for line in aug_genes[gid]:
                out.write(line)
            tier2_gap_count += 1

    print("Resolution Complete:")
    print(f"  Tier 1 (SQANTI3 empirical models retained) : {tier1_count}")
    print(f"  Tier 2 (Augustus novel gap-fillers added)  : {tier2_gap_count}")
    print(f"  Augustus redundant overlaps discarded     : {dropped_overlap_count}")
    print(f"  Augustus noisy monoexons discarded        : {dropped_monoexon_noise}")
    print(f"  Total non-redundant loci in resolved GTF   : {tier1_count + tier2_gap_count}")
    print(f"  Output written to: {output_gtf}")


# Direct Snakemake Execution Entry Point
if "snakemake" in globals():
    resolve_tiers(
        sqanti_gtf=snakemake.input.sqanti_gtf,
        augustus_gff=snakemake.input.augustus_gff,
        output_gtf=snakemake.output.resolved_gtf,
        hints_file=getattr(snakemake.input, "hints", None),
        min_monoexon_len=getattr(snakemake.params, "min_monoexon_len", 300),
        filter_mode=getattr(snakemake.params, "filter_mode", "medium"),
    )
elif __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Resolve Tier 1 SQANTI3 models with Tier 2 Augustus gap-filling predictions."
    )
    parser.add_argument("--sqanti", required=True, help="Input SQANTI3 filtered GTF")
    parser.add_argument("--augustus", required=True, help="Input Augustus prediction GFF/GTF")
    parser.add_argument("--output", required=True, help="Output resolved GTF")
    parser.add_argument("--hints", default=None, help="Optional extrinsic hints GFF")
    parser.add_argument("--min-monoexon-len", type=int, default=300, help="Min length for monoexons")
    parser.add_argument("--filter-mode", default="medium", choices=["strict", "medium", "none"])

    args = parser.parse_args()
    resolve_tiers(
        sqanti_gtf=args.sqanti,
        augustus_gff=args.augustus,
        output_gtf=args.output,
        hints_file=args.hints,
        min_monoexon_len=args.min_monoexon_len,
        filter_mode=args.filter_mode,
    )
