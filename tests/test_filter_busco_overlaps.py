#!/usr/bin/env python3
"""
Unit tests for filter_busco_overlaps.py

To run:
    pytest tests/test_filter_busco_overlaps.py
"""
import os
import sys
import tempfile
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Add scripts directory to PYTHONPATH
SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

import filter_busco_overlaps as fbo


def test_filter_busco_overlaps_direct_and_flanking():
    """
    Test direct overlap dropping, flanking proximity dropping (< 1000 bp),
    and retention of clean, distant BUSCO models.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # 1. SQANTI GFF (chr1: 1000..3000)
        sqanti_gff = os.path.join(tmpdir, "sqanti.gff")
        with open(sqanti_gff, "w") as f:
            f.write('chr1\tSQANTI3\ttranscript\t1000\t3000\t.\t+\t.\tgene_id "sqanti_gene1"; transcript_id "tx1";\n')
            f.write('chr1\tSQANTI3\tCDS\t1000\t3000\t.\t+\t0\tgene_id "sqanti_gene1"; transcript_id "tx1";\n')

        # 2. BUSCO GFF:
        # - busco_direct_inside (chr1: 1500..2500) -> should drop (direct overlap)
        # - busco_direct_partial (chr1: 500..1200) -> should drop (direct overlap)
        # - busco_too_close (chr1: 3500..4500) -> distance = 500 bp (< 1000 bp) -> should drop (flanking)
        # - busco_clean (chr1: 6000..7000) -> distance = 3000 bp (>= 1000 bp) -> should KEEP
        # - busco_diff_chr (chr2: 1000..2000) -> different chromosome -> should KEEP
        busco_gff = os.path.join(tmpdir, "busco.gff")
        with open(busco_gff, "w") as f:
            f.write('chr1\tminiprot\tCDS\t1500\t2500\t.\t+\t0\tID=busco_direct_inside.t1;Parent=busco_direct_inside\n')
            f.write('chr1\tminiprot\tCDS\t500\t1200\t.\t+\t0\tID=busco_direct_partial.t1;Parent=busco_direct_partial\n')
            f.write('chr1\tminiprot\tCDS\t3500\t4500\t.\t+\t0\tID=busco_too_close.t1;Parent=busco_too_close\n')
            f.write('chr1\tminiprot\tCDS\t6000\t7000\t.\t+\t0\tID=busco_clean.t1;Parent=busco_clean\n')
            f.write('chr2\tminiprot\tCDS\t1000\t2000\t.\t+\t0\tID=busco_diff_chr.t1;Parent=busco_diff_chr\n')

        # 3. BUSCO FAA
        busco_faa = os.path.join(tmpdir, "busco.faa")
        records = [
            SeqRecord(Seq("MAAA"), id="busco_direct_inside", description=""),
            SeqRecord(Seq("MBBB"), id="busco_direct_partial", description=""),
            SeqRecord(Seq("MCCC"), id="busco_too_close", description=""),
            SeqRecord(Seq("MDDD"), id="busco_clean", description=""),
            SeqRecord(Seq("MEEE"), id="busco_diff_chr", description=""),
        ]
        SeqIO.write(records, busco_faa, "fasta")

        out_clean_gff = os.path.join(tmpdir, "busco_clean.gff")
        out_clean_faa = os.path.join(tmpdir, "busco_clean.faa")

        # Run logic
        sqanti_intervals = fbo.extract_gene_intervals_from_gff(sqanti_gff)
        busco_genes, busco_lines = fbo.parse_busco_gff_records(busco_gff)
        retained_ids = fbo.filter_busco_genes(busco_genes, sqanti_intervals, flanking_size=1000)

        assert set(retained_ids) == {"busco_clean", "busco_diff_chr"}

        # Write outputs
        with open(out_clean_gff, "w") as out_gff:
            for gid in retained_ids:
                for line in busco_lines[gid]:
                    out_gff.write(line)

        fbo.filter_and_write_faa(busco_faa, set(retained_ids), out_clean_faa)

        # Verify clean GFF
        with open(out_clean_gff) as f:
            gff_text = f.read()
        assert "busco_clean" in gff_text
        assert "busco_diff_chr" in gff_text
        assert "busco_direct_inside" not in gff_text
        assert "busco_direct_partial" not in gff_text
        assert "busco_too_close" not in gff_text

        # Verify clean FAA
        clean_faa_records = list(SeqIO.parse(out_clean_faa, "fasta"))
        assert len(clean_faa_records) == 2
        assert {r.id for r in clean_faa_records} == {"busco_clean", "busco_diff_chr"}
