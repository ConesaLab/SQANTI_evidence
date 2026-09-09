#!/usr/bin/env python3
"""
Convert the SQANTI3-filtered GTF (Tier 1 transcripts) into Augustus hints with src=lrRNA.

Hint types are switched on/off in `envs/hint_config.tsv` (features: exon, exonpart, intron, CDS,
start, stop). Rules:
- exon/exonpart: one hint per exon; terminal exons become `exonpart` when that type is enabled.
- intron: between consecutive exons. When no exon-type hint is enabled, introns are restricted to the
  CDS span of the transcript (deliberate: UTR introns are dropped, see May 2026 commits).
- CDS: one hint per CDS block; the terminal block becomes `CDSpart` when the ORF lacks the start
  (5prime_partial/internal) or the stop (3prime_partial/internal) according to SQANTI's `CDS_type`.
- start/stop: on the first/last codon of the CDS, only for ORFs that have them. SQANTI CDS coordinates
  include the stop codon (verified on Arabidopsis: 3,300/3,302 complete ORFs), so the stop hint lands on it.

Column 6 is always `0` (numeric; Augustus ignores it with numscoreclasses=1). Hints of one transcript
share `grp=<transcript_id>;pri=1;src=lrRNA`.

Transcripts are collected in a dictionary keyed by transcript_id, so the GTF does not need to be
grouped or sorted; output order follows first appearance, which is byte-identical to the previous
streaming implementation on SQANTI's grouped GTFs.
"""
import sys
from collections import OrderedDict

KNOWN_FEATURES = ("exon", "exonpart", "intron", "CDS", "start", "stop")
HINT_SOURCE = "Hints"


def load_hint_config(config_file):
    """
    Reads envs/hint_config.tsv (header `feature<TAB>enabled`, values true/false).
    Blank lines are ignored; unknown features or malformed lines are errors.
    """
    hint_config = {f: False for f in KNOWN_FEATURES}
    with open(config_file, "r") as f:
        for lineno, raw in enumerate(f, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if lineno == 1 and parts[0].lower() == "feature":
                continue  # header
            if len(parts) != 2:
                raise ValueError(f"{config_file}:{lineno}: expected 'feature<TAB>true|false', got {raw.rstrip()!r}")
            feat, enabled = parts[0].strip(), parts[1].strip().lower()
            if feat not in KNOWN_FEATURES:
                raise ValueError(f"{config_file}:{lineno}: unknown hint feature {feat!r}; known: {', '.join(KNOWN_FEATURES)}")
            if enabled not in ("true", "false"):
                raise ValueError(f"{config_file}:{lineno}: value for {feat} must be true or false, got {parts[1]!r}")
            hint_config[feat] = enabled == "true"
    return hint_config


def load_cds_types(table_file, sep="\t"):
    """isoform -> CDS_type (lower-case) from the SQANTI classification table."""
    cds_dict = {}
    with open(table_file, "r") as f:
        header = f.readline().strip("\n").split(sep)
        iso_idx = header.index("isoform")
        cds_idx = header.index("CDS_type")
        for line in f:
            parts = line.strip("\n").split(sep)
            if len(parts) > max(iso_idx, cds_idx):
                cds_dict[parts[iso_idx]] = parts[cds_idx].strip().lower()
    return cds_dict


def transcript_hints(t_id, chrom, strand, exons, cds, cds_dict, config):
    """Returns the hint lines (with newline) for one transcript."""
    if not exons and not cds:
        return []
    lines = []
    hint_attrs = f"grp={t_id};pri=1;src=lrRNA"
    exons = sorted(exons)
    cds = sorted(cds)

    cds_type = cds_dict.get(t_id, "internal")
    has_start = cds_type in ("complete", "3prime_partial")
    has_stop = cds_type in ("complete", "5prime_partial")

    # 1. Exon hints
    if config.get("exon", False) or config.get("exonpart", False):
        for i, (ex_start, ex_end) in enumerate(exons):
            if config.get("exonpart", False) and (i == 0 or i == len(exons) - 1):
                lines.append(f"{chrom}\t{HINT_SOURCE}\texonpart\t{ex_start}\t{ex_end}\t0\t{strand}\t.\t{hint_attrs}\n")
            elif config.get("exon", False):
                lines.append(f"{chrom}\t{HINT_SOURCE}\texon\t{ex_start}\t{ex_end}\t0\t{strand}\t.\t{hint_attrs}\n")

    # 2. Intron hints
    if config.get("intron", False):
        extract_all_introns = config.get("exon", False) or config.get("exonpart", False)
        min_cds = cds[0][0] if cds else None
        max_cds = cds[-1][1] if cds else None
        for i in range(len(exons) - 1):
            intron_start = exons[i][1] + 1
            intron_end = exons[i + 1][0] - 1
            if intron_start > intron_end:
                continue
            if extract_all_introns or (min_cds is not None and intron_start >= min_cds and intron_end <= max_cds):
                lines.append(f"{chrom}\t{HINT_SOURCE}\tintron\t{intron_start}\t{intron_end}\t0\t{strand}\t.\t{hint_attrs}\n")

    # 3. CDS hints
    if cds and config.get("CDS", False):
        for i, (c_start, c_end) in enumerate(cds):
            h_type = "CDS"
            if strand == "+":
                if (i == 0 and not has_start) or (i == len(cds) - 1 and not has_stop):
                    h_type = "CDSpart"
            else:
                if (i == len(cds) - 1 and not has_start) or (i == 0 and not has_stop):
                    h_type = "CDSpart"
            lines.append(f"{chrom}\t{HINT_SOURCE}\t{h_type}\t{c_start}\t{c_end}\t0\t{strand}\t.\t{hint_attrs}\n")

    # 4. Start and stop hints
    if cds:
        min_cds, max_cds = cds[0][0], cds[-1][1]
        if strand == "+":
            if has_start and config.get("start", False):
                lines.append(f"{chrom}\t{HINT_SOURCE}\tstart\t{min_cds}\t{min_cds + 2}\t0\t+\t0\t{hint_attrs}\n")
            if has_stop and config.get("stop", False):
                lines.append(f"{chrom}\t{HINT_SOURCE}\tstop\t{max_cds - 2}\t{max_cds}\t0\t+\t0\t{hint_attrs}\n")
        elif strand == "-":
            if has_start and config.get("start", False):
                lines.append(f"{chrom}\t{HINT_SOURCE}\tstart\t{max_cds - 2}\t{max_cds}\t0\t-\t0\t{hint_attrs}\n")
            if has_stop and config.get("stop", False):
                lines.append(f"{chrom}\t{HINT_SOURCE}\tstop\t{min_cds}\t{min_cds + 2}\t0\t-\t0\t{hint_attrs}\n")
    return lines


def collect_transcripts(gff_file):
    """
    transcript_id -> {chrom, strand, exons: [(s,e)], cds: [(s,e)]}, in order of first appearance.
    Lines without transcript_id are skipped. Works on grouped and ungrouped GTFs alike.
    """
    transcripts = OrderedDict()
    with open(gff_file, "r") as f:
        for line in f:
            if line[0] == "#" or not line.strip():
                continue
            parts = line.strip("\n").split("\t")
            if len(parts) < 9:
                continue
            try:
                t_id = parts[8].split('transcript_id "')[1].split('"')[0]
            except IndexError:
                continue
            rec = transcripts.get(t_id)
            if rec is None:
                rec = transcripts[t_id] = {"chrom": parts[0], "strand": parts[6], "exons": [], "cds": []}
            feature = parts[2]
            if feature == "exon":
                rec["exons"].append((int(parts[3]), int(parts[4])))
            elif feature == "CDS":
                rec["cds"].append((int(parts[3]), int(parts[4])))
    return transcripts


def generate_hints(gff_file, table_file, config_file, out_file_path):
    hint_config = load_hint_config(config_file)
    print("1. Loading CDS metadata into memory...")
    cds_dict = load_cds_types(table_file)
    print("2. Processing GFF and generating hints...")
    transcripts = collect_transcripts(gff_file)
    n_hints = 0
    with open(out_file_path, "w") as out:
        for t_id, rec in transcripts.items():
            lines = transcript_hints(t_id, rec["chrom"], rec["strand"], rec["exons"], rec["cds"], cds_dict, hint_config)
            out.writelines(lines)
            n_hints += len(lines)
    print(f"   {len(transcripts)} transcripts -> {n_hints} hints ({', '.join(k for k, v in hint_config.items() if v)} enabled)")
    return n_hints


def main():
    generate_hints(
        gff_file=snakemake.input.gtf,  # noqa: F821
        table_file=snakemake.input.classification,  # noqa: F821
        config_file=snakemake.input.hint_config,  # noqa: F821
        out_file_path=snakemake.output[0],  # noqa: F821
    )


if __name__ == "__main__":
    main()
