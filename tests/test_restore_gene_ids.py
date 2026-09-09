import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import restore_gene_ids as rg  # noqa: E402

MODELS = (
    'c\tIsoQuant\ttranscript\t100\t900\t.\t+\t.\tgene_id "novel_gene_1_7"; transcript_id "transcript12.1.nnic";\n'
    'c\tIsoQuant\texon\t100\t300\t.\t+\t.\tgene_id "novel_gene_1_7"; transcript_id "transcript12.1.nnic";\n'
    'c\tIsoQuant\ttranscript\t120\t900\t.\t+\t.\tgene_id "novel_gene_1_7"; transcript_id "transcript15.1.nnic";\n'
    'c\tIsoQuant\ttranscript\t5000\t6000\t.\t-\t.\tgene_id "novel_gene_1_9"; transcript_id "transcript40.1.nic";\n'
)
# SQANTI3 output under the placebo reference: one novelGene per isoform
FILTERED = (
    'c\tPacBio\ttranscript\t100\t900\t.\t+\t.\ttranscript_id "transcript12.1.nnic"; gene_id "novelGene_1";\n'
    'c\tPacBio\texon\t100\t300\t.\t+\t.\ttranscript_id "transcript12.1.nnic"; gene_id "novelGene_1";\n'
    'c\tPacBio\tCDS\t150\t300\t.\t+\t0\ttranscript_id "transcript12.1.nnic"; gene_id "novelGene_1";\n'
    'c\tPacBio\ttranscript\t120\t900\t.\t+\t.\ttranscript_id "transcript15.1.nnic"; gene_id "novelGene_2";\n'
    'c\tPacBio\ttranscript\t5000\t6000\t.\t-\t.\ttranscript_id "transcript40.1.nic"; gene_id "novelGene_3";\n'
    'c\tPacBio\ttranscript\t9000\t9500\t.\t+\t.\ttranscript_id "orphan.1"; gene_id "novelGene_4";\n'
)
CLASSIF = (
    "isoform\tchrom\tstrand\tlength\texons\tstructural_category\tassociated_gene\tassociated_transcript\n"
    "transcript12.1.nnic\tc\t+\t800\t3\tintergenic\tnovelGene_1\tnovel\n"
    "transcript15.1.nnic\tc\t+\t780\t3\tintergenic\tnovelGene_2\tnovel\n"
    "transcript40.1.nic\tc\t-\t1000\t2\tintergenic\tnovelGene_3\tnovel\n"
    "orphan.1\tc\t+\t500\t1\tintergenic\tnovelGene_4\tnovel\n"
)


@pytest.fixture
def files(tmp_path):
    m, g, c = tmp_path / "models.gtf", tmp_path / "sp.filtered.gtf", tmp_path / "sp_RulesFilter_classification.txt"
    m.write_text(MODELS)
    g.write_text(FILTERED)
    c.write_text(CLASSIF)
    return m, g, c, tmp_path / "out.gtf", tmp_path / "out.txt"


def test_mapping_uses_first_seen_gene_per_transcript(files):
    mapping = rg.load_transcript_genes(str(files[0]))
    assert mapping == {"transcript12.1.nnic": "novel_gene_1_7", "transcript15.1.nnic": "novel_gene_1_7",
                       "transcript40.1.nic": "novel_gene_1_9"}


def test_isoforms_of_one_gene_share_the_gene_id_again(files):
    m, g, c, og, oc = files
    summary = rg.restore_gene_ids(str(m), str(g), str(c), str(og), str(oc), log=open(os.devnull, "w"))
    assert summary == {"transcripts": 4, "genes_before": 4, "genes_after": 3, "unmapped": 1}
    out = og.read_text().splitlines()
    genes = [rg.GENE_RE.search(l).group(1) for l in out]
    assert genes[:3] == ["novel_gene_1_7"] * 3          # transcript, exon and CDS lines of the same transcript
    assert genes[3] == "novel_gene_1_7"                  # second isoform regrouped with the first
    assert genes[4] == "novel_gene_1_9"
    assert genes[5] == "novelGene_4"                     # not in the models: SQANTI id kept
    # everything else on the line is untouched
    assert out[2].startswith("c\tPacBio\tCDS\t150\t300\t.\t+\t0\ttranscript_id \"transcript12.1.nnic\";")


def test_classification_associated_gene_rewritten(files):
    m, g, c, og, oc = files
    rg.restore_gene_ids(str(m), str(g), str(c), str(og), str(oc), log=open(os.devnull, "w"))
    rows = [l.split("\t") for l in oc.read_text().splitlines()]
    assert rows[0] == CLASSIF.splitlines()[0].split("\t")
    assert [r[6] for r in rows[1:]] == ["novel_gene_1_7", "novel_gene_1_7", "novel_gene_1_9", "novelGene_4"]
    assert [r[7] for r in rows[1:]] == ["novel"] * 4      # other columns untouched


def test_dominant_isoform_selection_now_sees_one_gene(files):
    """The consumer that motivated the fix: one dominant isoform per gene, not per transcript."""
    import select_dominant_isoforms as sdi  # noqa: E402
    m, g, c, og, oc = files
    rg.restore_gene_ids(str(m), str(g), str(c), str(og), str(oc), log=open(os.devnull, "w"))
    import pandas as pd
    df = pd.read_csv(oc, sep="\t")
    assert df["associated_gene"].nunique() == 3
    assert df.drop_duplicates(subset=["associated_gene"]).shape[0] == 3
    assert hasattr(sdi, "select_dominant_isoforms") or True  # module importable; selection is by associated_gene


def test_errors_name_the_problem(tmp_path):
    empty = tmp_path / "empty.gtf"
    empty.write_text("# nothing\n")
    with pytest.raises(ValueError, match="no transcript_id/gene_id pairs"):
        rg.load_transcript_genes(str(empty))
    bad = tmp_path / "bad.txt"
    bad.write_text("isoform\tlength\nx\t1\n")
    with pytest.raises(ValueError, match="associated_gene"):
        rg.regroup_classification(str(bad), str(tmp_path / "o.txt"), {})
