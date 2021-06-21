import pytest
from lingrex.align import (
    gap_free_pairwise,
    align_to_template,
    shrink_alignments,
    template_alignment,
    shrink_template,
)
from lingpy import Wordlist


def test_gap_free_pairwise():

    seqA, seqB = list("andra"), list("an-ra")

    almA, almB = gap_free_pairwise(seqA, seqB)
    assert almA[1] == "n<d"

    almA, almB = gap_free_pairwise(seqA, seqB, syllables=[0, 2])
    assert almA[2] == "d>r"

    seqA, seqB = list("este"), list("-ste")
    almA, almB = gap_free_pairwise(seqA, seqB)
    assert almA[0] == "e>s"

    seqA, seqB = list("euste"), list("--ste")
    almA, almB = gap_free_pairwise(seqA, seqB)
    assert almA[0] == "e>u>s"


def test_align_to_template():
    out = align_to_template("ka", "Cv", "Cvc")
    assert out[-1] == "-"

    with pytest.raises(ValueError):
        align_to_template("ka", "c", "Cvc")
    with pytest.raises(ValueError):
        align_to_template("ka", "cv", "Cv")


def test_shrink_alignments():
    out = shrink_alignments([["a", "b", "-"], ["a", "b", "-"]])
    assert len(out[0]) == 2


@pytest.fixture
def wldata():
    return {
        0: ["doculect", "concept", "tokens", "structure", "cogid"],
        1: ["a", "a", "b au".split(), "i n".split(), 1],
        2: ["b", "a", "b o k".split(), "i n c".split(), 1],
        3: ["c", "a", "b w a k".split(), "i m n c".split(), 1],
    }


@pytest.fixture
def wldata_listvalued_cogid(wldata):
    return {k: v if k == 0 else v[:-1] + [[v[-1]]] for k, v in wldata.items()}


def test_template_alignment(wldata, wldata_listvalued_cogid):
    wl = Wordlist(wldata)
    template_alignment(wl, fuzzy=False, template="imnc")
    assert "alignment" in wl.columns
    wl = Wordlist(wldata_listvalued_cogid)
    template_alignment(wl, fuzzy=True, template="imnc")
    assert "alignment" in wl.columns


def test_shrink_template(wldata_listvalued_cogid):
    wl = Wordlist(wldata_listvalued_cogid)
    shrink_template(wl)
    assert wl[2, "tokens2"][-1] == "ok"
