#!/usr/bin/env python3
"""
Script to select dominant complete-ORF isoforms per gene from SQANTI3 filtered results
and extract GFF records and translated protein sequences.

Author: Pablo Atienza & Gemini
"""
import sys
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


def select_dominant_transcripts(classification_file: str) -> dict:
    """
    Filter SQANTI post-filter classification table for complete coding ORFs (filter_result == 'Isoform')
    and select the single dominant isoform per gene locus.

    Args:
        classification_file: Path to SQANTI3 post-filter classification TSV.

    Returns:
        Dictionary mapping selected transcript_id to associated gene_id.
    """
    sys.stderr.write(f"### Loading SQANTI classification: [{classification_file}]\n")
    df = pd.read_csv(classification_file, sep="\t", low_memory=False)

    # Normalize column names
    df.columns = [c.strip() for c in df.columns]

    # Filter for transcripts passing SQANTI3 filter (filter_result == 'Isoform')
    if "filter_result" in df.columns:
        df = df[df["filter_result"].astype(str).str.strip().str.lower() == "isoform"].copy()
        sys.stderr.write(f"### Transcripts passing filter (filter_result == 'Isoform'): {len(df)}\n")

    # Filter for complete coding transcripts
    if "coding" in df.columns:
        df = df[df["coding"] == "coding"]
    if "CDS_type" in df.columns:
        df = df[df["CDS_type"] == "complete"]

    sys.stderr.write(f"### Complete coding ORFs: {len(df)}\n")
    if df.empty:
        sys.stderr.write("WARNING: No complete coding transcripts found after filtering.\n")
        return {}

    # Read count column 'FL'
    if "FL" in df.columns:
        df["FL"] = pd.to_numeric(df["FL"], errors="coerce").fillna(0)
    else:
        df["FL"] = 0

    # Ensure CDS_length is numeric
    if "CDS_length" in df.columns:
        df["CDS_length"] = pd.to_numeric(df["CDS_length"], errors="coerce").fillna(0)
    else:
        df["CDS_length"] = pd.to_numeric(df.get("length", 0), errors="coerce").fillna(0)

    # Determine gene grouping column
    gene_col = "associated_gene" if "associated_gene" in df.columns else "gene_id"
    if gene_col not in df.columns or df[gene_col].isna().all():
        df[gene_col] = df["isoform"].apply(lambda x: x.rsplit(".", 1)[0] if "." in str(x) else str(x))

    # Sort by FL count desc, then CDS_length desc
    df_sorted = df.sort_values(by=["FL", "CDS_length"], ascending=[False, False])

    # Select dominant (first) isoform per gene locus
    dominant_df = df_sorted.drop_duplicates(subset=[gene_col], keep="first")
    sys.stderr.write(f"### Selected dominant loci: {len(dominant_df)}\n")

    return dict(zip(dominant_df["isoform"], dominant_df[gene_col]))


def extract_and_write_data(
    input_gtf: str,
    selected_txs: dict,
    genome_fasta: str,
    output_gff: str,
    output_faa: str,
):
    """
    Extract GTF records for selected transcripts and translate their CDS to protein FASTA.
    """
    sys.stderr.write(f"### Loading Genome FASTA: [{genome_fasta}]\n")
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    cds_coords = {}  # tx_id -> list of (chrom, strand, start, end)
    gff_lines = []

    sys.stderr.write(f"### Filtering GTF records from: [{input_gtf}]\n")
    with open(input_gtf, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            attributes = fields[8]
            # Extract transcript_id
            tx_id = None
            for attr in attributes.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tx_id = attr.split("transcript_id")[-1].strip().strip('"').strip("'").strip()
                    break

            if tx_id and tx_id in selected_txs:
                gff_lines.append(line)
                feature = fields[2]
                if feature == "CDS":
                    chrom = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    if tx_id not in cds_coords:
                        cds_coords[tx_id] = []
                    cds_coords[tx_id].append((chrom, strand, start, end))

    # Write filtered GFF
    sys.stderr.write(f"### Writing dominant GFF: [{output_gff}]\n")
    with open(output_gff, "w") as outfile:
        for line in gff_lines:
            outfile.write(line)

    # Translate CDS and write protein FASTA
    sys.stderr.write(f"### Translating CDS to FAA: [{output_faa}]\n")
    written_proteins = 0
    with open(output_faa, "w") as out_faa:
        for tx_id, coords in cds_coords.items():
            if not coords:
                continue
            chrom = coords[0][0]
            strand = coords[0][1]
            if chrom not in genome:
                sys.stderr.write(f"WARNING: Chromosome {chrom} not found in genome FASTA.\n")
                continue

            # Sort coordinates
            if strand == "+":
                coords.sort(key=lambda x: x[2])
                cds_nuc = "".join(str(genome[chrom].seq[s - 1 : e]) for _, _, s, e in coords)
            else:
                coords.sort(key=lambda x: x[2], reverse=True)
                cds_nuc = "".join(
                    str(genome[chrom].seq[s - 1 : e].reverse_complement())
                    for _, _, s, e in coords
                )

            # Translate CDS
            try:
                prot_seq = str(Seq(cds_nuc).translate(to_stop=True))
            except Exception as e:
                prot_seq = str(Seq(cds_nuc).translate()).rstrip("*")

            gene_id = selected_txs.get(tx_id, tx_id)
            out_faa.write(f">{gene_id} {tx_id}\n{prot_seq}\n")
            written_proteins += 1

    sys.stderr.write(f"### Wrote {written_proteins} protein sequences to FAA\n")


def main():
    gtf_file = snakemake.input.gtf
    class_file = snakemake.input.classification
    genome_file = snakemake.input.genome
    out_gff = snakemake.output.gff
    out_faa = snakemake.output.faa

    # Ensure output directories exist
    os.makedirs(os.path.dirname(out_gff), exist_ok=True)
    os.makedirs(os.path.dirname(out_faa), exist_ok=True)

    selected_txs = select_dominant_transcripts(class_file)
    extract_and_write_data(gtf_file, selected_txs, genome_file, out_gff, out_faa)


if __name__ == "__main__":
    main()
