#!/usr/bin/env python3
"""
Unit tests for resolve_transcript_tiers.py

To run:
    pytest tests/test_resolve_transcript_tiers.py -v
"""
import os
import sys
import tempfile
import pytest
from types import SimpleNamespace

# Add scripts directory to PYTHONPATH
SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

import resolve_transcript_tiers as rtt


@pytest.fixture
def mock_sqanti_gtf():
    content = """chr1\tSQANTI3\ttranscript\t1000\t2500\t.\t+\t.\ttranscript_id "sq_tx1"; gene_id "sq_gene1";
chr1\tSQANTI3\texon\t1000\t1400\t.\t+\t.\ttranscript_id "sq_tx1"; gene_id "sq_gene1";
chr1\tSQANTI3\texon\t1800\t2500\t.\t+\t.\ttranscript_id "sq_tx1"; gene_id "sq_gene1";
chr1\tSQANTI3\tCDS\t1100\t1400\t.\t+\t0\ttranscript_id "sq_tx1"; gene_id "sq_gene1";
chr1\tSQANTI3\tCDS\t1800\t2300\t.\t+\t1\ttranscript_id "sq_tx1"; gene_id "sq_gene1";
chr1\tSQANTI3\ttranscript\t5000\t7000\t.\t-\t.\ttranscript_id "sq_tx2"; gene_id "sq_gene2";
chr1\tSQANTI3\texon\t5000\t7000\t.\t-\t.\ttranscript_id "sq_tx2"; gene_id "sq_gene2";
chr1\tSQANTI3\tCDS\t5200\t6800\t.\t-\t0\ttranscript_id "sq_tx2"; gene_id "sq_gene2";
"""
    return content


@pytest.fixture
def mock_augustus_gff():
    content = """# start gene g1
chr1\tAUGUSTUS\tgene\t1050\t2400\t1\t+\t.\tg1
chr1\tAUGUSTUS\ttranscript\t1100\t2300\t1\t+\t.\tg1.t1
chr1\tAUGUSTUS\texon\t1100\t1400\t1\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";
chr1\tAUGUSTUS\texon\t1800\t2300\t1\t+\t.\ttranscript_id "g1.t1"; gene_id "g1";
chr1\tAUGUSTUS\tCDS\t1100\t1400\t1\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";
chr1\tAUGUSTUS\tCDS\t1800\t2300\t1\t+\t1\ttranscript_id "g1.t1"; gene_id "g1";
# end gene g1
# start gene g2_gap
chr1\tAUGUSTUS\tgene\t3000\t4500\t1\t+\t.\tg2_gap
chr1\tAUGUSTUS\ttranscript\t3000\t4500\t1\t+\t.\tg2_gap.t1
chr1\tAUGUSTUS\texon\t3000\t3500\t1\t+\t.\ttranscript_id "g2_gap.t1"; gene_id "g2_gap";
chr1\tAUGUSTUS\texon\t4000\t4500\t1\t+\t.\ttranscript_id "g2_gap.t1"; gene_id "g2_gap";
chr1\tAUGUSTUS\tCDS\t3000\t3500\t1\t+\t0\ttranscript_id "g2_gap.t1"; gene_id "g2_gap";
chr1\tAUGUSTUS\tCDS\t4000\t4500\t1\t+\t2\ttranscript_id "g2_gap.t1"; gene_id "g2_gap";
# end gene g2_gap
# start gene g3_mono_noise
chr1\tAUGUSTUS\tgene\t8000\t8200\t0.5\t+\t.\tg3_mono_noise
chr1\tAUGUSTUS\ttranscript\t8000\t8200\t0.5\t+\t.\tg3_mono_noise.t1
chr1\tAUGUSTUS\texon\t8000\t8200\t0.5\t+\t.\ttranscript_id "g3_mono_noise.t1"; gene_id "g3_mono_noise";
chr1\tAUGUSTUS\tCDS\t8000\t8200\t0.5\t+\t0\ttranscript_id "g3_mono_noise.t1"; gene_id "g3_mono_noise";
# end gene g3_mono_noise
# start gene g4_mono_long
chr1\tAUGUSTUS\tgene\t9000\t9600\t0.9\t+\t.\tg4_mono_long
chr1\tAUGUSTUS\ttranscript\t9000\t9600\t0.9\t+\t.\tg4_mono_long.t1
chr1\tAUGUSTUS\texon\t9000\t9600\t0.9\t+\t.\ttranscript_id "g4_mono_long.t1"; gene_id "g4_mono_long";
chr1\tAUGUSTUS\tCDS\t9000\t9600\t0.9\t+\t0\ttranscript_id "g4_mono_long.t1"; gene_id "g4_mono_long";
# end gene g4_mono_long
"""
    return content


def test_resolve_transcript_tiers_logic(mock_sqanti_gtf, mock_augustus_gff):
    """Test core Tier 1 pass-through, Tier 2 overlap discard, and gap-filling."""
    with tempfile.TemporaryDirectory() as tmpdir:
        sqanti_path = os.path.join(tmpdir, "sqanti.filtered.gtf")
        augustus_path = os.path.join(tmpdir, "augustus.gff")
        out_path = os.path.join(tmpdir, "resolved.gtf")

        with open(sqanti_path, "w") as f:
            f.write(mock_sqanti_gtf)
        with open(augustus_path, "w") as f:
            f.write(mock_augustus_gff)

        rtt.resolve_tiers(
            sqanti_gtf=sqanti_path,
            augustus_gff=augustus_path,
            output_gtf=out_path,
            min_monoexon_len=300,
            filter_mode="medium",
        )

        with open(out_path) as f:
            lines = [l.strip() for l in f if l.strip()]

        # Check that SQANTI models (sq_gene1 and sq_gene2) are present
        sq1_found = any('gene_id "sq_gene1"' in l for l in lines)
        sq2_found = any('gene_id "sq_gene2"' in l for l in lines)
        assert sq1_found, "Tier 1 sq_gene1 must be preserved"
        assert sq2_found, "Tier 1 sq_gene2 must be preserved"

        # Check that overlapping Augustus g1 is discarded
        g1_found = any('gene_id "g1"' in l for l in lines)
        assert not g1_found, "Overlapping Augustus gene g1 must be discarded"

        # Check that non-overlapping multi-exon Augustus g2_gap is retained
        g2_found = any('gene_id "g2_gap"' in l for l in lines)
        assert g2_found, "Non-overlapping gap-filler g2_gap must be retained"

        # Check that short monoexon noise g3 (<300 bp) is discarded
        g3_found = any('gene_id "g3_mono_noise"' in l for l in lines)
        assert not g3_found, "Short monoexon noise g3 must be discarded"

        # Check that long monoexon g4 (601 bp >= 300 bp) is retained
        g4_found = any('gene_id "g4_mono_long"' in l for l in lines)
        assert g4_found, "Long monoexon g4 must be retained under medium mode"


def test_resolve_transcript_tiers_hints_support(mock_sqanti_gtf, mock_augustus_gff):
    """Test that monoexons with hint support are retained even in strict mode."""
    with tempfile.TemporaryDirectory() as tmpdir:
        sqanti_path = os.path.join(tmpdir, "sqanti.filtered.gtf")
        augustus_path = os.path.join(tmpdir, "augustus.gff")
        hints_path = os.path.join(tmpdir, "hints.gff")
        out_path = os.path.join(tmpdir, "resolved.gtf")

        with open(sqanti_path, "w") as f:
            f.write(mock_sqanti_gtf)
        with open(augustus_path, "w") as f:
            f.write(mock_augustus_gff)
        with open(hints_path, "w") as f:
            # Add hint overlapping g3_mono_noise (8000..8200)
            f.write("chr1\tHints\tstart\t8000\t8002\t0\t+\t0\tgrp=h1;pri=1;src=P\n")

        rtt.resolve_tiers(
            sqanti_gtf=sqanti_path,
            augustus_gff=augustus_path,
            output_gtf=out_path,
            hints_file=hints_path,
            min_monoexon_len=300,
            filter_mode="strict",
        )

        with open(out_path) as f:
            lines = [l.strip() for l in f if l.strip()]

        # In strict mode, g3_mono_noise has a hint so it should be retained
        g3_found = any('gene_id "g3_mono_noise"' in l for l in lines)
        assert g3_found, "Hint-supported monoexon must be retained in strict mode"

        # g4_mono_long has NO hints, so in strict mode it should be dropped
        g4_found = any('gene_id "g4_mono_long"' in l for l in lines)
        assert not g4_found, "Unsupported monoexon must be dropped in strict mode"


def test_snakemake_object_integration(mock_sqanti_gtf, mock_augustus_gff):
    """Test execution via a mock Snakemake context."""
    with tempfile.TemporaryDirectory() as tmpdir:
        sqanti_path = os.path.join(tmpdir, "sqanti.filtered.gtf")
        augustus_path = os.path.join(tmpdir, "augustus.gff")
        out_path = os.path.join(tmpdir, "resolved.gtf")

        with open(sqanti_path, "w") as f:
            f.write(mock_sqanti_gtf)
        with open(augustus_path, "w") as f:
            f.write(mock_augustus_gff)

        # Create mock snakemake object
        mock_snakemake = SimpleNamespace(
            input=SimpleNamespace(
                sqanti_gtf=sqanti_path,
                augustus_gff=augustus_path,
                hints=None,
            ),
            output=SimpleNamespace(
                resolved_gtf=out_path,
            ),
            params=SimpleNamespace(
                min_monoexon_len=300,
                filter_mode="medium",
            ),
        )

        # Inject into module and execute core function
        rtt.resolve_tiers(
            sqanti_gtf=mock_snakemake.input.sqanti_gtf,
            augustus_gff=mock_snakemake.input.augustus_gff,
            output_gtf=mock_snakemake.output.resolved_gtf,
            hints_file=getattr(mock_snakemake.input, "hints", None),
            min_monoexon_len=getattr(mock_snakemake.params, "min_monoexon_len", 300),
            filter_mode=getattr(mock_snakemake.params, "filter_mode", "medium"),
        )

        assert os.path.isfile(out_path)
        with open(out_path) as f:
            content = f.read()
            assert 'gene_id "sq_gene1"' in content
            assert 'gene_id "g2_gap"' in content
