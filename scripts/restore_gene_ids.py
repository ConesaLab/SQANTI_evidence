#!/usr/bin/env python3
"""Restore the transcript-model gene grouping after the SQANTI3 rules filter.

SQANTI3 assigns a fresh ``novelGene_<n>`` to every isoform it classifies as
``intergenic`` or ``genic_intron`` (``src/helpers.py::rename_novel_genes``: "I don't find it
necessary to cluster these novel genes"). With the placebo reference every isoform is
intergenic, so the filtered GTF carries one gene per transcript and the isoform-per-gene
structure produced by the transcript modeller (IsoQuant) is lost. Downstream this inflates
gene counts and makes every coding isoform a "dominant" isoform of its own gene.

This step maps each surviving transcript back to the ``gene_id`` it had in the transcript
models GTF and rewrites ``gene_id`` in the filtered GTF and ``associated_gene`` in the
filtered classification. Transcripts absent from the models GTF keep their SQANTI id and are
counted in the log.

Runs as a Snakemake ``script:`` (uses the ``snakemake`` object) or standalone::

    restore_gene_ids.py --models isoquant.transcript_models.gtf --gtf sp.filtered.gtf \
        --classification sp_RulesFilter_classification.txt \
        --out-gtf sp.filtered.regrouped.gtf --out-classification sp_RulesFilter_classification.regrouped.txt
"""
import argparse
import re
import sys

TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')
GENE_RE = re.compile(r'gene_id\s+"([^"]+)"')


def load_transcript_genes(models_gtf):
    """Return {transcript_id: gene_id} from a GTF of transcript models (any feature line)."""
    mapping = {}
    with open(models_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            tx = TX_RE.search(fields[8])
            gene = GENE_RE.search(fields[8])
            if tx and gene and tx.group(1) not in mapping:
                mapping[tx.group(1)] = gene.group(1)
    if not mapping:
        raise ValueError(f"no transcript_id/gene_id pairs found in {models_gtf}")
    return mapping


def regroup_gtf(gtf_in, gtf_out, mapping):
    """Rewrite gene_id per transcript; returns (transcripts, genes_before, genes_after, unmapped)."""
    transcripts, unmapped = set(), set()
    genes_before, genes_after = set(), set()
    with open(gtf_in) as fin, open(gtf_out, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                fout.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                fout.write(line)
                continue
            attrs = fields[8]
            tx = TX_RE.search(attrs)
            old = GENE_RE.search(attrs)
            if tx is None or old is None:
                fout.write(line)
                continue
            transcripts.add(tx.group(1))
            genes_before.add(old.group(1))
            new_gene = mapping.get(tx.group(1))
            if new_gene is None:
                unmapped.add(tx.group(1))
                new_gene = old.group(1)
            genes_after.add(new_gene)
            fields[8] = GENE_RE.sub(f'gene_id "{new_gene}"', attrs, count=1)
            fout.write("\t".join(fields) + "\n")
    return len(transcripts), len(genes_before), len(genes_after), sorted(unmapped)


def regroup_classification(class_in, class_out, mapping):
    """Rewrite the associated_gene column; returns (rows, rewritten)."""
    rows = rewritten = 0
    with open(class_in) as fin, open(class_out, "w") as fout:
        header = fin.readline()
        fout.write(header)
        cols = header.rstrip("\n").split("\t")
        try:
            iso_i = cols.index("isoform")
            gene_i = cols.index("associated_gene")
        except ValueError as exc:
            raise ValueError(f"{class_in}: missing column {exc}") from None
        for line in fin:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            rows += 1
            new_gene = mapping.get(fields[iso_i])
            if new_gene is not None:
                fields[gene_i] = new_gene
                rewritten += 1
            fout.write("\t".join(fields) + "\n")
    return rows, rewritten


def restore_gene_ids(models_gtf, gtf_in, class_in, gtf_out, class_out, log=sys.stderr):
    mapping = load_transcript_genes(models_gtf)
    n_tx, g_before, g_after, unmapped = regroup_gtf(gtf_in, gtf_out, mapping)
    rows, rewritten = regroup_classification(class_in, class_out, mapping)
    print(f"Transcript models: {len(mapping)} transcripts in {len(set(mapping.values()))} genes", file=log)
    print(f"Filtered GTF: {n_tx} transcripts; genes {g_before} -> {g_after} "
          f"({n_tx / g_after:.2f} transcripts per gene)" if g_after else "Filtered GTF: empty", file=log)
    print(f"Classification: {rewritten} of {rows} rows rewritten", file=log)
    if unmapped:
        print(f"WARNING: {len(unmapped)} transcripts not found in the transcript models, SQANTI gene_id kept: "
              + ", ".join(unmapped[:10]) + (" ..." if len(unmapped) > 10 else ""), file=log)
    return {"transcripts": n_tx, "genes_before": g_before, "genes_after": g_after, "unmapped": len(unmapped)}


def main_snakemake():
    with open(snakemake.log[0], "w") as log:  # noqa: F821 (injected by Snakemake)
        restore_gene_ids(
            models_gtf=snakemake.input.models,  # noqa: F821
            gtf_in=snakemake.input.gtf,  # noqa: F821
            class_in=snakemake.input.classification,  # noqa: F821
            gtf_out=snakemake.output.gtf,  # noqa: F821
            class_out=snakemake.output.classification,  # noqa: F821
            log=log,
        )


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--models", required=True, help="transcript models GTF (IsoQuant) with the original gene grouping")
    p.add_argument("--gtf", required=True, help="SQANTI3 filtered GTF")
    p.add_argument("--classification", required=True, help="SQANTI3 filtered (RulesFilter) classification")
    p.add_argument("--out-gtf", required=True)
    p.add_argument("--out-classification", required=True)
    a = p.parse_args(argv)
    restore_gene_ids(a.models, a.gtf, a.classification, a.out_gtf, a.out_classification)


if __name__ == "__main__":
    if "snakemake" in globals():
        main_snakemake()
    else:
        main()
