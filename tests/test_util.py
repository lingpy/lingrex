import pytest
from lingrex.util import lingrex_path, add_structure
from lingpy import Wordlist, Alignments
from lingrex.util import ungap, clean_sound, unjoin, alm2tok, bleu_score


def test_bleu_score():
    candidate = "this is a test".split()
    reference = "this is a small test".split()

    assert round(
            bleu_score(
                candidate, 
                reference, 
                weights=[0.5, 0.5],
                n=2,
                trim=True
                ),
            2) == 0.64

    assert round(
            bleu_score(
                candidate,
                reference,
                weights=[0.5, 0.5],
                n=2,
                trim=False),
            2) == 0.70

    

def test_ungap():

    matrix = ungap([['a', 'b'], ['x', '-'], ['y', '-']], ['proto', 'l1', 'l2'],
            'proto')
    assert round(
        bleu_score(
            candidate,
            reference,
            n=2,
            trim=False),
        2) == 0.70


def test_ungap():
    matrix = ungap([['a', 'b'], ['x', '-'], ['y', '-']], ['proto', 'l1', 'l2'], 'proto')
    assert matrix[0][0] == 'a.b'
    assert matrix[1][0] == 'x'
    assert matrix[2][0] == "y"
    matrix2 = ungap([['a', 'b'], ['x', '-'], ['y', 'h']], ['proto', 'l1', 'l2'], 'proto')
    assert matrix2[0][1] == ["a", "b"][1]
    assert matrix2[1][1] == ["x", "-"][1]
    assert matrix2[2][1] == ["y", "h"][1]

    out = ungap([["p", "-", "a"], ["p", "j", "a"]], ["German", "E"], "E")
    assert out[1][0] == "p.j"

def test_clean_sound():
    
    alm = [['a', 'b'], ['-', '-'], ['-', '-']]
    assert ungap(alm, ['p', 'l1', 'l2'], 'p') == alm

def test_clean_sound():
    assert clean_sound("a/b") == "b"
    assert clean_sound("a") == "a"
    assert clean_sound("a/b.c/d") == "b.d"


def test_unjoin():
    assert unjoin("k.p a p u k.a/b".split())[0] == "k"


def test_lingrex_path():
    lingrex_path("test")


def test_add_structure():

    with pytest.raises(ValueError):
        add_structure(
            Wordlist(
                {
                    0: ["doculect", "concept", "tokens", "cogid"],
                    1: ["a", "b", "b l a".split(), 1],
                    2: ["b", "b", "b l a x".split(), 1],
                    3: ["c", "b", "b l i k u s".split(), 1],
                }
            ),
            model="bla",
        )

    for m in ["cv", "c", "CcV", "nogap", "ps"]:
        D = {
            0: ["doculect", "concept", "tokens", "cogid"],
            1: ["a", "b", "b l a".split(), 1],
            2: ["b", "b", "b l a x".split(), 1],
            3: ["c", "b", "b l i k u s".split(), 1],
            4: ["d", "b", "b l u k", 2],
        }
        wl = Alignments(D, transcription="tokens")
        add_structure(wl, m)

    for m in ["cv", "c", "CcV", "nogap", "ps"]:
        D = {
            0: ["doculect", "concept", "tokens", "cogids"],
            1: ["a", "b", "b l a".split(), [1]],
            2: ["b", "b", "b l a x".split(), [1]],
            3: ["c", "b", "b l i k u s".split(), [1]],
        }
        wl = Alignments(D, ref="cogids", transcription="tokens")
        add_structure(wl, m, ref="cogids")
