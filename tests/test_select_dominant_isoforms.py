#!/usr/bin/env python3
"""
Unit tests for select_dominant_isoforms.py

To run:
    pytest tests/test_select_dominant_isoforms.py
"""
import os
import sys
import tempfile
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq

# Add scripts directory to PYTHONPATH
SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

import select_dominant_isoforms as sdi


@pytest.fixture
def sample_classification_data():
    """Generates a mock SQANTI classification DataFrame with diverse test cases."""
    return pd.DataFrame([
        # Gene 1: Two passing isoforms -> tx1.1 (FL=10) should beat tx1.2 (FL=2)
        {
            "isoform": "tx1.1",
            "associated_gene": "gene1",
            "filter_result": "Isoform",
            "coding": "coding",
            "CDS_type": "complete",
            "FL": 10,
            "CDS_length": 45,
        },
        {
            "isoform": "tx1.2",
            "associated_gene": "gene1",
            "filter_result": "Isoform",
            "coding": "coding",
            "CDS_type": "complete",
            "FL": 2,
            "CDS_length": 60,
        },
        # Gene 1: Artifact isoform -> should be ignored even though FL is high
        {
            "isoform": "tx1.3",
            "associated_gene": "gene1",
            "filter_result": "Artifact",
            "coding": "coding",
            "CDS_type": "complete",
            "FL": 999,
            "CDS_length": 60,
        },
        # Gene 2: Incomplete CDS -> should be ignored
        {
            "isoform": "tx2.1",
            "associated_gene": "gene2",
            "filter_result": "Isoform",
            "coding": "coding",
            "CDS_type": "5prime_partial",
            "FL": 50,
            "CDS_length": 90,
        },
        # Gene 3: Non-coding transcript -> should be ignored
        {
            "isoform": "tx3.1",
            "associated_gene": "gene3",
            "filter_result": "Isoform",
            "coding": "non_coding",
            "CDS_type": "complete",
            "FL": 20,
            "CDS_length": 0,
        },
        # Gene 4: Equal FL count -> tx4.2 (CDS_length=90) should beat tx4.1 (CDS_length=45)
        {
            "isoform": "tx4.1",
            "associated_gene": "gene4",
            "filter_result": "Isoform",
            "coding": "coding",
            "CDS_type": "complete",
            "FL": 15,
            "CDS_length": 45,
        },
        {
            "isoform": "tx4.2",
            "associated_gene": "gene4",
            "filter_result": "Isoform",
            "coding": "coding",
            "CDS_type": "complete",
            "FL": 15,
            "CDS_length": 90,
        },
    ])


def test_select_dominant_transcripts(sample_classification_data):
    """Test filtering for complete coding ORFs and selecting dominant isoforms."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp_class:
        sample_classification_data.to_csv(tmp_class.name, sep="\t", index=False)
        tmp_class_path = tmp_class.name

    try:
        selected = sdi.select_dominant_transcripts(tmp_class_path)
        # Expected winners:
        # gene1 -> tx1.1 (FL=10 beats FL=2)
        # gene4 -> tx4.2 (CDS_length=90 beats CDS_length=45 on FL tie)
        assert "tx1.1" in selected
        assert selected["tx1.1"] == "gene1"
        assert "tx1.2" not in selected
        assert "tx1.3" not in selected
        assert "tx2.1" not in selected
        assert "tx3.1" not in selected
        assert "tx4.2" in selected
        assert selected["tx4.2"] == "gene4"
        assert "tx4.1" not in selected
        assert len(selected) == 2
    finally:
        if os.path.exists(tmp_class_path):
            os.remove(tmp_class_path)


def test_extract_and_write_data_forward_and_reverse_strand():
    """Test GFF writing and CDS protein translation on both + and - strands."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Mock reference genome (chr1: 200 bp)
        # Forward CDS: ATGAAACCC (9bp) + TTTAAATAA (9bp) = 18 bp -> MKP FKN* -> MKPFK
        # Coordinates: 10..18 and 31..39
        # Reverse CDS: complement of ATGAAACCCTTTAAATAA -> TTATTTAAAGGGTTTCAT
        # Coordinates: 100..108 and 121..129
        fwd_exon1 = "ATGAAACCC"
        fwd_exon2 = "TTTAAATAA"
        rev_exon1_rc = "ATGAAACCC"  # transcribed from revcomp
        rev_exon2_rc = "TTTAAATAA"

        seq_list = ["N"] * 200
        # Forward exons
        for i, ch in enumerate(fwd_exon1):
            seq_list[9 + i] = ch
        for i, ch in enumerate(fwd_exon2):
            seq_list[30 + i] = ch

        # Reverse exons (store rev-comp in genome)
        rev_exon1_dna = str(Seq(rev_exon1_rc).reverse_complement())
        rev_exon2_dna = str(Seq(rev_exon2_rc).reverse_complement())
        for i, ch in enumerate(rev_exon2_dna):
            seq_list[99 + i] = ch
        for i, ch in enumerate(rev_exon1_dna):
            seq_list[120 + i] = ch

        genome_fasta = os.path.join(tmpdir, "genome.fa")
        with open(genome_fasta, "w") as f:
            f.write(f">chr1\n{''.join(seq_list)}\n")

        # Mock GTF
        gtf_file = os.path.join(tmpdir, "sample.gtf")
        with open(gtf_file, "w") as f:
            # tx1 (+)
            f.write('chr1\tSQANTI3\ttranscript\t10\t39\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";\n')
            f.write('chr1\tSQANTI3\texon\t10\t18\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";\n')
            f.write('chr1\tSQANTI3\tCDS\t10\t18\t.\t+\t0\tgene_id "gene1"; transcript_id "tx1";\n')
            f.write('chr1\tSQANTI3\texon\t31\t39\t.\t+\t.\tgene_id "gene1"; transcript_id "tx1";\n')
            f.write('chr1\tSQANTI3\tCDS\t31\t39\t.\t+\t0\tgene_id "gene1"; transcript_id "tx1";\n')
            # tx2 (-)
            f.write('chr1\tSQANTI3\ttranscript\t100\t129\t.\t-\t.\tgene_id "gene2"; transcript_id "tx2";\n')
            f.write('chr1\tSQANTI3\texon\t100\t108\t.\t-\t.\tgene_id "gene2"; transcript_id "tx2";\n')
            f.write('chr1\tSQANTI3\tCDS\t100\t108\t.\t-\t0\tgene_id "gene2"; transcript_id "tx2";\n')
            f.write('chr1\tSQANTI3\texon\t121\t129\t.\t-\t.\tgene_id "gene2"; transcript_id "tx2";\n')
            f.write('chr1\tSQANTI3\tCDS\t121\t129\t.\t-\t0\tgene_id "gene2"; transcript_id "tx2";\n')

        out_gff = os.path.join(tmpdir, "output.gff")
        out_faa = os.path.join(tmpdir, "output.faa")

        selected_txs = {"tx1": "gene1", "tx2": "gene2"}
        sdi.extract_and_write_data(gtf_file, selected_txs, genome_fasta, out_gff, out_faa)

        # Check GFF
        with open(out_gff, "r") as f:
            gff_lines = f.readlines()
        assert len(gff_lines) == 10

        # Check FAA
        records = list(SeqIO.parse(out_faa, "fasta"))
        assert len(records) == 2
        assert records[0].id == "gene1"
        assert str(records[0].seq) == "MKPFK"
        assert records[1].id == "gene2"
        assert str(records[1].seq) == "MKPFK"
