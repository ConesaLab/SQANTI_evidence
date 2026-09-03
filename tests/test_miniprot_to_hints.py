#!/usr/bin/env python3
"""
Unit tests for miniprot_to_hints.py

To run:
    pytest tests/test_miniprot_to_hints.py
"""
import os
import sys
import tempfile
import pytest

# Add scripts directory to PYTHONPATH
SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

import miniprot_to_hints as mth


def test_miniprot_to_hints_positive_strand():
    """Test Miniprot parsing and Augustus hint generation on positive strand."""
    with tempfile.TemporaryDirectory() as tmpdir:
        miniprot_gff = os.path.join(tmpdir, "miniprot.gff")
        with open(miniprot_gff, "w") as f:
            f.write("chr1\tminiprot\tmRNA\t1000\t3000\t150\t+\t.\tID=MP001;Target=prot_gene1 1 500;Identity=0.98\n")
            f.write("chr1\tminiprot\tCDS\t1000\t1200\t.\t+\t0\tParent=MP001;Target=prot_gene1 1 67\n")
            f.write("chr1\tminiprot\tCDS\t1800\t2200\t.\t+\t1\tParent=MP001;Target=prot_gene1 68 200\n")
            f.write("chr1\tminiprot\tCDS\t2800\t3000\t.\t+\t2\tParent=MP001;Target=prot_gene1 201 267\n")

        out_hints = os.path.join(tmpdir, "protein.hints.gff")

        alignments = mth.parse_miniprot_gff(miniprot_gff)
        assert len(alignments) == 1
        assert alignments[0]["target"] == "prot_gene1"
        assert len(alignments[0]["cds"]) == 3

        mth.write_augustus_protein_hints(alignments, out_hints, src="P", priority=2)

        with open(out_hints) as f:
            lines = [l.strip() for l in f if l.strip()]

        # Check hints: 3 CDSpart + 2 introns + 1 start + 1 stop = 7 hints
        assert len(lines) == 7

        # Check CDSpart phases
        assert "chr1\tminiprot\tCDSpart\t1000\t1200\t0\t+\t0\tgrp=prot_gene1;pri=2;src=P" in lines
        assert "chr1\tminiprot\tCDSpart\t1800\t2200\t0\t+\t1\tgrp=prot_gene1;pri=2;src=P" in lines
        assert "chr1\tminiprot\tCDSpart\t2800\t3000\t0\t+\t2\tgrp=prot_gene1;pri=2;src=P" in lines

        # Check introns
        assert "chr1\tminiprot\tintron\t1201\t1799\t0\t+\t.\tgrp=prot_gene1;pri=2;src=P" in lines
        assert "chr1\tminiprot\tintron\t2201\t2799\t0\t+\t.\tgrp=prot_gene1;pri=2;src=P" in lines

        # Check start & stop on + strand
        assert "chr1\tminiprot\tstart\t1000\t1002\t0\t+\t0\tgrp=prot_gene1;pri=3;src=P" in lines
        assert "chr1\tminiprot\tstop\t2998\t3000\t0\t+\t0\tgrp=prot_gene1;pri=3;src=P" in lines


def test_miniprot_to_hints_negative_strand():
    """Test Miniprot parsing and Augustus hint generation on negative strand."""
    with tempfile.TemporaryDirectory() as tmpdir:
        miniprot_gff = os.path.join(tmpdir, "miniprot.gff")
        with open(miniprot_gff, "w") as f:
            f.write("chr2\tminiprot\tmRNA\t5000\t7000\t200\t-\t.\tID=MP002;Target=prot_gene2 1 400;Identity=0.95\n")
            f.write("chr2\tminiprot\tCDS\t5000\t5500\t.\t-\t0\tParent=MP002;Target=prot_gene2 201 367\n")
            f.write("chr2\tminiprot\tCDS\t6500\t7000\t.\t-\t0\tParent=MP002;Target=prot_gene2 1 167\n")

        out_hints = os.path.join(tmpdir, "protein.hints.gff")
        alignments = mth.parse_miniprot_gff(miniprot_gff)
        mth.write_augustus_protein_hints(alignments, out_hints, src="P", priority=2)

        with open(out_hints) as f:
            lines = [l.strip() for l in f if l.strip()]

        # 2 CDSpart + 1 intron + 1 start + 1 stop = 5 hints
        assert len(lines) == 5

        # Check intron
        assert "chr2\tminiprot\tintron\t5501\t6499\t0\t-\t.\tgrp=prot_gene2;pri=2;src=P" in lines

        # On '-' strand: start is at highest coordinate (6998..7000), stop is at lowest coordinate (5000..5002)
        assert "chr2\tminiprot\tstart\t6998\t7000\t0\t-\t0\tgrp=prot_gene2;pri=3;src=P" in lines
        assert "chr2\tminiprot\tstop\t5000\t5002\t0\t-\t0\tgrp=prot_gene2;pri=3;src=P" in lines
