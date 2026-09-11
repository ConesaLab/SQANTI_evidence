import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import fastq2bam as f2b  # noqa: E402

pysam = pytest.importorskip("pysam", reason="pysam lives in envs/isoseq.yaml, not the test environment")
pytest.importorskip("Bio", reason="biopython lives in envs/isoseq.yaml, not the test environment")

TEMPLATE = os.path.join(os.path.dirname(__file__), "..", "envs", "pacbio_mock.bam")

FASTQ = "@read_one\nACGTACGTAC\n+\nIIIIIIIIII\n@read_two\nTTTTGGGGCCCCAA\n+\nIIIIIIIIIIIIII\n"
FASTA = ">read_one\nACGTACGTAC\n>read_two\nTTTTGGGGCCCCAA\n"


def test_detect_format(tmp_path):
    fq, fa = tmp_path / "r.fastq", tmp_path / "r.fasta"
    fq.write_text(FASTQ)
    fa.write_text(FASTA)
    assert f2b.detect_format(str(fq)) == "fastq"
    assert f2b.detect_format(str(fa)) == "fasta"


@pytest.mark.parametrize("content,msg", [
    ("", "file is empty"),
    ("not a sequence file\n", "neither"),
])
def test_detect_format_errors_name_the_problem(tmp_path, content, msg):
    bad = tmp_path / "bad.txt"
    bad.write_text(content)
    with pytest.raises(ValueError, match=msg):
        f2b.detect_format(str(bad))


def test_read_identity_and_sequences_survive_the_conversion(tmp_path):
    """Names, sequences and qualities are the reads' own; everything else comes from the template."""
    reads, out = tmp_path / "r.fastq", tmp_path / "out.bam"
    reads.write_text(FASTQ)
    n = f2b.convert(str(reads), TEMPLATE, str(out), log=open(os.devnull, "w"))
    assert n == 2
    with pysam.AlignmentFile(str(out), "rb", check_sq=False) as bam:
        records = list(bam)
        assert [r.query_name for r in records] == ["read_one", "read_two"]
        assert [r.query_sequence for r in records] == ["ACGTACGTAC", "TTTTGGGGCCCCAA"]
        assert records[1].query_qualities[0] == 40                      # 'I' == Q40
        assert dict(records[1].tags).get("qe") == 14                    # qe tracks the new length
        assert "PACBIO" in str(bam.header)                              # template header reused


def test_fasta_input_gets_placeholder_qualities(tmp_path):
    reads, out = tmp_path / "r.fasta", tmp_path / "out.bam"
    reads.write_text(FASTA)
    f2b.convert(str(reads), TEMPLATE, str(out), log=open(os.devnull, "w"))
    with pysam.AlignmentFile(str(out), "rb", check_sq=False) as bam:
        first = next(iter(bam))
        assert set(first.query_qualities) == {93}                       # '~', i.e. synthetic


def test_truncated_record_fails_loudly(tmp_path):
    """A malformed FASTQ must stop the rule, not yield a silently short BAM.

    Biopython rejects the record before the 'no reads found' guard in convert() is reached; the guard
    stays as a backstop for parsers that would return an empty iterator instead.
    """
    reads, out = tmp_path / "r.fastq", tmp_path / "out.bam"
    reads.write_text("@only_a_name\n")
    with pytest.raises(ValueError):
        f2b.convert(str(reads), TEMPLATE, str(out), log=open(os.devnull, "w"))


def test_missing_template_is_reported(tmp_path):
    reads = tmp_path / "r.fastq"
    reads.write_text(FASTQ)
    with pytest.raises(Exception):
        f2b.convert(str(reads), str(tmp_path / "absent.bam"), str(tmp_path / "out.bam"),
                    log=open(os.devnull, "w"))
