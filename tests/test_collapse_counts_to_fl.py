import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import collapse_counts_to_fl as cc  # noqa: E402

ABUNDANCE = (
    "#count_fl: Number of input FL reads associated with isoform\n"
    "#norm_fl : fl_assoc / total number of FL reads supporting a mapped isoform\n"
    "#fl_assoc: Total number of FL reads associated with isosform\n"
    "pbid\tcount_fl\tfl_assoc\tcell_barcodes\n"
    "PB.1.1\t1\t3\tNA\n"
    "PB.2.1\t1\t4\tNA\n"
    "PB.3.2\t2\t13\tNA\n"
)


@pytest.fixture
def files(tmp_path):
    ab = tmp_path / "sample.collapsed.abundance.txt"
    ab.write_text(ABUNDANCE)
    return ab, tmp_path / "sample.fl_counts.tsv"


def test_writes_the_two_column_matrix_sqanti_expects(files):
    ab, out = files
    n, total = cc.convert(str(ab), str(out), log=open(os.devnull, "w"))
    assert (n, total) == (3, 20)
    assert out.read_text().splitlines() == ["feature_id\tcount", "PB.1.1\t3", "PB.2.1\t4", "PB.3.2\t13"]


def test_uses_fl_assoc_not_count_fl(files):
    """count_fl counts clusters; fl_assoc counts reads, which is what IsoQuant's column means."""
    ab, out = files
    cc.convert(str(ab), str(out), log=open(os.devnull, "w"))
    counts = [l.split("\t")[1] for l in out.read_text().splitlines()[1:]]
    assert counts == ["3", "4", "13"]          # fl_assoc
    assert counts != ["1", "1", "2"]           # count_fl


def test_comment_lines_are_not_mistaken_for_the_header(files):
    ab, out = files
    rows = list(cc.read_abundance(str(ab)))
    assert rows[0] == ("PB.1.1", "3") and len(rows) == 3


@pytest.mark.parametrize("content,msg", [
    ("pbid\tcount_fl\tcell_barcodes\nPB.1.1\t1\tNA\n", "missing column 'fl_assoc'"),
    ("count_fl\tfl_assoc\n1\t3\n", "missing column 'pbid'"),
    ("#only a comment\n", "no header found"),
])
def test_malformed_abundance_names_the_problem(tmp_path, content, msg):
    ab = tmp_path / "bad.txt"
    ab.write_text(content)
    with pytest.raises(ValueError, match=msg):
        list(cc.read_abundance(str(ab)))


def test_empty_transcript_list_raises_instead_of_writing_an_empty_matrix(tmp_path):
    ab = tmp_path / "empty.txt"
    ab.write_text("pbid\tcount_fl\tfl_assoc\tcell_barcodes\n")
    with pytest.raises(ValueError, match="no transcripts listed"):
        cc.convert(str(ab), str(tmp_path / "out.tsv"), log=open(os.devnull, "w"))


def test_non_numeric_count_is_reported_with_its_transcript(tmp_path):
    ab = tmp_path / "bad.txt"
    ab.write_text("pbid\tcount_fl\tfl_assoc\tcell_barcodes\nPB.1.1\t1\tmany\tNA\n")
    with pytest.raises(ValueError, match="non-numeric fl_assoc 'many' for PB.1.1"):
        cc.convert(str(ab), str(tmp_path / "out.tsv"), log=open(os.devnull, "w"))
