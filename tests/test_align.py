from lingrex.align import (
        gap_free_pairwise, align_to_template,
        shrink_alignments,
        template_alignment,
        shrink_template
        )
from lingpy import Wordlist
import pytest

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
    out = shrink_alignments(
            [["a", "b", "-"],
            ["a", "b", "-"]])
    assert len(out[0]) == 2


def test_template_alignment():
    wl = Wordlist(
            {
                0: ["doculect", "concept", "tokens", "structure", "cogid"],
                1: ["a", "a", "b au".split(), "i n".split(), 1],
                2: ["b", "a", "b o k".split(), "i n c".split(), 1],
                3: ["c", "a", "b w a k".split(), "i m n c".split(), 1]}
            )
    template_alignment(wl, fuzzy=False, template="imnc")
    assert "alignment" in wl.columns
    wl = Wordlist(
            {
                0: ["doculect", "concept", "tokens", "structure", "cogid"],
                1: ["a", "a", "b au".split(), "i n".split(), [1]],
                2: ["b", "a", "b o k".split(), "i n c".split(), [1]],
                3: ["c", "a", "b w a k".split(), "i m n c".split(), [1]]}
            )
    template_alignment(wl, fuzzy=True, template="imnc")
    assert "alignment" in wl.columns

def test_shrink_template():
    wl = Wordlist(
            {
                0: ["doculect", "concept", "tokens", "structure", "cogid"],
                1: ["a", "a", "b au".split(), "i n".split(), [1]],
                2: ["b", "a", "b o k".split(), "i n c".split(), [1]],
                3: ["c", "a", "b w a k".split(), "i m n c".split(), [1]]}
            )
    shrink_template(wl)
    assert wl[2, "tokens2"][-1] == "ok"







