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


# ---------------------------------------------------------------------------
# Real-shaped ids (roadmap 1.7b): BUSCO's miniprot GFF uses ID=MP######/Parent=MP######
# and identifies the BUSCO through Target=<busco_id>_<taxid>_<...>; the FAA headers are the
# bare BUSCO id. The old key (MP id) never matched a FAA header -> empty FAA.
# ---------------------------------------------------------------------------
REAL_BUSCO_GFF = (
    "5\tminiprot\tmRNA\t20152898\t20155574\t2551\t+\t.\tID=MP062483;Rank=1;Identity=0.8260;Target=10052at3699_29727_0:004f4a 1 586\n"
    "5\tminiprot\tCDS\t20152898\t20152956\t101\t+\t0\tParent=MP062483;Rank=1;Identity=1.0000;Target=10052at3699_29727_0:004f4a 1 19\n"
    "5\tminiprot\tCDS\t20153075\t20155574\t943\t+\t1\tParent=MP062483;Rank=1;Identity=0.8211;Target=10052at3699_29727_0:004f4a 20 586\n"
    "5\tminiprot\tstop_codon\t20155572\t20155574\t.\t+\t0\tParent=MP062483;Rank=1\n"
    "2\tminiprot\tmRNA\t15449839\t15451307\t607\t-\t.\tID=MP249840;Rank=1;Identity=0.9600;Target=10055at3699_90675_0:00349a 1 125\n"
    "2\tminiprot\tCDS\t15450518\t15451307\t212\t-\t0\tParent=MP249840;Rank=1;Identity=0.9184;Target=10055at3699_90675_0:00349a 1 125\n"
    "2\tminiprot\tstop_codon\t15450518\t15450520\t.\t-\t0\tParent=MP249840;Rank=1\n"
)


def test_busco_gene_key_from_target():
    assert fbo.busco_gene_key("ID=MP062483;Rank=1;Target=10052at3699_29727_0:004f4a 1 586") == "10052at3699"
    assert fbo.busco_gene_key("Parent=MP062483;Rank=1") is None
    assert fbo.generic_gene_key("Parent=MP062483;Rank=1") == "MP062483"
    assert fbo.generic_gene_key('gene_id "novelGene_1"; transcript_id "t1";') == "novelGene_1"


def test_real_busco_ids_reach_the_faa(tmp_path):
    sqanti = tmp_path / "sqanti.gff"
    # SQANTI gene overlapping the chr2 BUSCO only
    sqanti.write_text('2\tSQANTI3\tCDS\t15450000\t15452000\t.\t-\t0\tgene_id "novelGene_9"; transcript_id "t9";\n')
    busco_gff = tmp_path / "busco.gff"; busco_gff.write_text(REAL_BUSCO_GFF)
    busco_faa = tmp_path / "busco.faa"
    SeqIO.write([SeqRecord(Seq("MAAA"), id="10052at3699", description=""),
                 SeqRecord(Seq("MBBB"), id="10055at3699", description="")], str(busco_faa), "fasta")
    genes, lines = fbo.parse_busco_gff_records(str(busco_gff))
    assert set(genes) == {"10052at3699", "10055at3699"}          # keyed by BUSCO id, not MP id
    assert len(lines["10052at3699"]) == 4                        # mRNA + 2 CDS + stop_codon kept together
    retained = fbo.filter_busco_genes(genes, fbo.extract_gene_intervals_from_gff(str(sqanti)), 1000)
    assert retained == ["10052at3699"]
    out_faa = tmp_path / "clean.faa"
    n = fbo.filter_and_write_faa(str(busco_faa), set(retained), str(out_faa))
    assert n == 1 and [r.id for r in SeqIO.parse(str(out_faa), "fasta")] == ["10052at3699"]


def test_id_mismatch_between_gff_and_faa_raises(tmp_path):
    busco_faa = tmp_path / "busco.faa"
    SeqIO.write([SeqRecord(Seq("MAAA"), id="something_else", description="")], str(busco_faa), "fasta")
    with pytest.raises(ValueError, match="identifiers disagree"):
        fbo.filter_and_write_faa(str(busco_faa), {"10052at3699"}, str(tmp_path / "clean.faa"))
