#!/usr/bin/env python3
"""
Unit tests for assemble_training_gff.py

To run:
    pytest tests/test_assemble_training_gff.py
"""
import os
import sys
import tempfile
import pytest

# Add scripts directory to PYTHONPATH
SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

import assemble_training_gff as atg


def test_assemble_training_genes_busco_core():
    """Test Option 1: BUSCO core guaranteed + SQANTI majority."""
    cdhit_ids = {"busco1", "busco2", "sq1", "sq2", "sq3"}
    busco_order = ["busco1", "busco2"]
    busco_dict = {"busco1": ["b1_line"], "busco2": ["b2_line"]}
    sqanti_order = ["sq1", "sq2", "sq3"]
    sqanti_dict = {"sq1": ["s1_line"], "sq2": ["s2_line"], "sq3": ["s3_line"]}

    # Case 1: max_genes = 4 (takes 2 BUSCO + 2 SQANTI)
    b_sel, s_sel = atg.assemble_training_genes(
        cdhit_ids, busco_order, busco_dict, sqanti_order, sqanti_dict, max_genes=4, strategy="busco_core"
    )
    assert b_sel == ["busco1", "busco2"]
    assert s_sel == ["sq1", "sq2"]

    # Case 2: max_genes = 10 (takes all 2 BUSCO + all 3 SQANTI = 5)
    b_sel, s_sel = atg.assemble_training_genes(
        cdhit_ids, busco_order, busco_dict, sqanti_order, sqanti_dict, max_genes=10, strategy="busco_core"
    )
    assert b_sel == ["busco1", "busco2"]
    assert s_sel == ["sq1", "sq2", "sq3"]

    # Case 3: max_genes = 1 (takes 1 BUSCO)
    b_sel, s_sel = atg.assemble_training_genes(
        cdhit_ids, busco_order, busco_dict, sqanti_order, sqanti_dict, max_genes=1, strategy="busco_core"
    )
    assert b_sel == ["busco1"]
    assert s_sel == []


def test_assemble_training_genes_busco_and_sqanti_only():
    """Test dedicated single-mode strategies (busco_only and sqanti_only)."""
    cdhit_ids = {"busco1", "busco2", "sq1", "sq2", "sq3"}
    busco_order = ["busco1", "busco2"]
    busco_dict = {"busco1": ["b1_line"], "busco2": ["b2_line"]}
    sqanti_order = ["sq1", "sq2", "sq3"]
    sqanti_dict = {"sq1": ["s1_line"], "sq2": ["s2_line"], "sq3": ["s3_line"]}

    # busco_only
    b_sel, s_sel = atg.assemble_training_genes(
        cdhit_ids, busco_order, busco_dict, sqanti_order, sqanti_dict, max_genes=5000, strategy="busco_only"
    )
    assert b_sel == ["busco1", "busco2"]
    assert s_sel == []

    # sqanti_only
    b_sel, s_sel = atg.assemble_training_genes(
        cdhit_ids, busco_order, busco_dict, sqanti_order, sqanti_dict, max_genes=5000, strategy="sqanti_only"
    )
    assert b_sel == []
    assert s_sel == ["sq1", "sq2", "sq3"]


def test_assemble_training_genes_alternative_strategies():
    """Test future fallback strategies (sqanti_priority and proportional)."""
    cdhit_ids = {"busco1", "busco2", "sq1", "sq2", "sq3"}
    busco_order = ["busco1", "busco2"]
    busco_dict = {"busco1": ["b1_line"], "busco2": ["b2_line"]}
    sqanti_order = ["sq1", "sq2", "sq3"]
    sqanti_dict = {"sq1": ["s1_line"], "sq2": ["s2_line"], "sq3": ["s3_line"]}

    # SQANTI priority with max_genes = 4 (takes 3 SQANTI + 1 BUSCO)
    b_sel, s_sel = atg.assemble_training_genes(
        cdhit_ids, busco_order, busco_dict, sqanti_order, sqanti_dict, max_genes=4, strategy="sqanti_priority"
    )
    assert s_sel == ["sq1", "sq2", "sq3"]
    assert b_sel == ["busco1"]


def test_end_to_end_gff_assembly():
    """Test full file parsing and assembled GFF writing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        cdhit_file = os.path.join(tmpdir, "cdhit.lst")
        with open(cdhit_file, "w") as f:
            f.write("busco_1\nsqanti_1\nsqanti_2\n")

        busco_gff = os.path.join(tmpdir, "busco.gff")
        with open(busco_gff, "w") as f:
            f.write('chr1\tminiprot\tCDS\t100\t200\t.\t+\t0\tID=busco_1.t1;Parent=busco_1\n')

        sqanti_gff = os.path.join(tmpdir, "sqanti.gff")
        with open(sqanti_gff, "w") as f:
            f.write('chr1\tSQANTI3\tCDS\t1000\t2000\t.\t+\t0\tgene_id "sqanti_1"; transcript_id "tx1";\n')
            f.write('chr1\tSQANTI3\tCDS\t3000\t4000\t.\t+\t0\tgene_id "sqanti_2"; transcript_id "tx2";\n')

        out_gff = os.path.join(tmpdir, "training.gff")

        cdhit_ids = atg.parse_cdhit_list(cdhit_file)
        b_order, b_dict = atg.parse_gff_by_gene(busco_gff)
        s_order, s_dict = atg.parse_gff_by_gene(sqanti_gff)

        b_sel, s_sel = atg.assemble_training_genes(
            cdhit_ids, b_order, b_dict, s_order, s_dict, max_genes=5000, strategy="busco_core"
        )

        with open(out_gff, "w") as out_f:
            for gid in b_sel:
                for line in b_dict[gid]:
                    out_f.write(line)
            for gid in s_sel:
                for line in s_dict[gid]:
                    out_f.write(line)

        with open(out_gff) as f:
            gff_content = f.read()

        assert "busco_1" in gff_content
        assert "sqanti_1" in gff_content
        assert "sqanti_2" in gff_content


def test_real_busco_ids_match_cdhit_representatives(tmp_path):
    """BUSCO GFF keyed by Target= BUSCO id must match CD-HIT ids taken from the FAA headers (roadmap 1.7b)."""
    busco_gff = tmp_path / "busco.gff"
    busco_gff.write_text(
        "5\tminiprot\tmRNA\t100\t900\t.\t+\t.\tID=MP062483;Rank=1;Target=10052at3699_29727_0:004f4a 1 200\n"
        "5\tminiprot\tCDS\t100\t900\t.\t+\t0\tParent=MP062483;Rank=1;Target=10052at3699_29727_0:004f4a 1 200\n"
        "2\tminiprot\tmRNA\t100\t900\t.\t-\t.\tID=MP249840;Rank=1;Target=10055at3699_90675_0:00349a 1 125\n"
        "2\tminiprot\tCDS\t100\t900\t.\t-\t0\tParent=MP249840;Rank=1;Target=10055at3699_90675_0:00349a 1 125\n"
    )
    sqanti_gff = tmp_path / "sqanti.gff"
    sqanti_gff.write_text('1\tSQANTI3\tCDS\t1000\t2000\t.\t+\t0\ttranscript_id "t1"; gene_id "novelGene_1";\n')
    cdhit = tmp_path / "cdhit.lst"
    cdhit.write_text("novelGene_1 t1\n10052at3699\n")   # 10055at3699 clustered away
    b_order, b_dict = atg.parse_gff_by_gene(str(busco_gff))
    s_order, s_dict = atg.parse_gff_by_gene(str(sqanti_gff))
    assert b_order == ["10052at3699", "10055at3699"]
    sel_b, sel_s = atg.assemble_training_genes(atg.parse_cdhit_list(str(cdhit)), b_order, b_dict, s_order, s_dict,
                                               max_genes=5000, strategy="busco_core")
    assert sel_b == ["10052at3699"] and sel_s == ["novelGene_1"]


def test_unknown_strategy_raises():
    with pytest.raises(ValueError, match="Unknown training-set assembly strategy"):
        atg.assemble_training_genes({"b1"}, ["b1"], {"b1": []}, [], {}, max_genes=10, strategy="random_pick")
