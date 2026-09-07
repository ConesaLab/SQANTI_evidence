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


VALID_FILTER_MODES = ("strict", "medium", "none")


def build_interval_index(intervals_dict):
    """
    Builds a nested-interval-safe search structure per (chrom, strand).

    Intervals are sorted by start; alongside the sorted starts we keep a running maximum of the
    ends (max_ends[i] = largest end among intervals 0..i). Overlap with a query [qs, qe] is then:
        idx = bisect_right(starts, qe)        # intervals starting at or before the query end
        idx > 0 and max_ends[idx - 1] >= qs   # the furthest-reaching one of them reaches the query
    The running maximum carries a long gene's end past every gene nested inside it, which the old
    implementation (bisect over the unsorted list of ends) could miss.
    Accepts values as (start, end, ...) tuples; extra tuple fields are ignored.
    """
    index = {}
    for key, ivs in intervals_dict.items():
        ordered = sorted(ivs, key=lambda x: (x[0], x[1]))
        starts, max_ends, running = [], [], -1
        for iv in ordered:
            running = max(running, iv[1])
            starts.append(iv[0])
            max_ends.append(running)
        index[key] = (starts, max_ends)
    return index


def index_overlaps(index, chrom, strand, start, end):
    """True if [start, end] overlaps any interval of index[(chrom, strand)]."""
    entry = index.get((chrom, strand))
    if not entry:
        return False
    starts, max_ends = entry
    idx = bisect.bisect_right(starts, end)
    return idx > 0 and max_ends[idx - 1] >= start


def check_hint_support(chrom, strand, start, end, hint_intervals):
    """Checks if interval overlaps any hint. Accepts the raw dict from load_hints_coords or a prebuilt index."""
    return index_overlaps(_as_index(hint_intervals), chrom, strand, start, end)


def check_locus_overlap(chrom, strand, start, end, intervals_dict):
    """
    Checks if genomic interval [start, end] overlaps any interval in intervals_dict[(chrom, strand)].
    Accepts the raw dict from parse_sqanti_gtf or a prebuilt index (see build_interval_index).
    """
    return index_overlaps(_as_index(intervals_dict), chrom, strand, start, end)


def _as_index(obj):
    """Returns obj if it already is an index, else builds one (convenience for direct calls/tests)."""
    if obj and all(isinstance(v, tuple) and len(v) == 2 and isinstance(v[0], list) for v in obj.values()):
        return obj
    return build_interval_index(obj or {})


def resolve_tiers(
    sqanti_gtf,
    augustus_gff,
    output_gtf,
    hints_file=None,
    min_monoexon_len=300,
    filter_mode="medium",
    log_path=None,
):
    """
    Core resolution function:
      1. Writes 100% of SQANTI3 curated models to output_gtf.
      2. For each Augustus model:
         - Discards if overlapping any SQANTI3 locus on the same strand.
         - For non-overlapping multi-exon: retains.
         - For non-overlapping single-exon: checks length / hint support.
    """
    if filter_mode not in VALID_FILTER_MODES:
        raise ValueError(f"filter_mode must be one of {VALID_FILTER_MODES}, got {filter_mode!r}")

    messages = []

    def log(msg):
        print(msg)
        messages.append(msg)

    log(f"Loading Tier 1 SQANTI3 models from: {sqanti_gtf}")
    sqanti_genes, sqanti_raw = parse_sqanti_gtf(sqanti_gtf)
    sqanti_intervals = build_interval_index(sqanti_raw)
    log(f"  Loaded {len(sqanti_genes)} SQANTI3 gene loci.")

    log(f"Loading Tier 2 Augustus predictions from: {augustus_gff}")
    aug_genes, aug_info = parse_augustus_gff(augustus_gff)
    log(f"  Loaded {len(aug_info)} Augustus gene predictions.")

    hint_intervals = build_interval_index(load_hints_coords(hints_file) if hints_file else {})
    log(f"Tier 2 single-exon filter: {filter_mode} (min length {min_monoexon_len} bp)")

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

    log("Resolution Complete:")
    log(f"  Tier 1 (SQANTI3 empirical models retained) : {tier1_count}")
    log(f"  Tier 2 (Augustus novel gap-fillers added)  : {tier2_gap_count}")
    log(f"  Augustus redundant overlaps discarded     : {dropped_overlap_count}")
    log(f"  Augustus noisy monoexons discarded        : {dropped_monoexon_noise}")
    log(f"  Total non-redundant loci in resolved GTF   : {tier1_count + tier2_gap_count}")
    log(f"  Output written to: {output_gtf}")

    if log_path:
        os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
        with open(log_path, "w") as log_f:
            log_f.write("\n".join(messages) + "\n")

    return {
        "tier1": tier1_count,
        "tier2": tier2_gap_count,
        "dropped_overlap": dropped_overlap_count,
        "dropped_monoexon": dropped_monoexon_noise,
    }


def main_snakemake(smk):
    """Entry point used by the Snakemake `script:` directive (also callable from tests)."""
    log_files = getattr(smk, "log", None)
    return resolve_tiers(
        sqanti_gtf=smk.input.sqanti_gtf,
        augustus_gff=smk.input.augustus_gff,
        output_gtf=smk.output.resolved_gtf,
        hints_file=getattr(smk.input, "hints", None),
        min_monoexon_len=int(getattr(smk.params, "min_monoexon_len", 300)),
        filter_mode=getattr(smk.params, "filter_mode", "medium"),
        log_path=log_files[0] if log_files else None,
    )


# Direct Snakemake Execution Entry Point
if "snakemake" in globals():
    main_snakemake(snakemake)  # noqa: F821
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
