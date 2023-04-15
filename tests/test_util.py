import pytest
from lingpy import Wordlist, Alignments
from lingrex.util import lingrex_path, add_structure
from lingrex.util import ungap, clean_sound, unjoin, alm2tok, bleu_score
from lingrex.util import prep_wordlist, subsequence_of


@pytest.mark.parametrize(
    'source,target,result',
    [
        ('cvc', 'cvcvc', True),
        ('bla', 'bla', True),
        ('bla', 'bxlyaz', True),
        ('bla', 'abxlyaz', True),
        ('abc', 'ab', False),
    ]
)
def test_subsequence_of(source, target, result):
    assert subsequence_of(source, target) == result


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


def test_prep_wordlist(wl_with_alignments):
    test_wl = prep_wordlist(Wordlist(wl_with_alignments))

    assert len(test_wl) == 4
    assert "+" not in test_wl[1, "tokens"]
    assert "_" not in test_wl[2, "tokens"]
