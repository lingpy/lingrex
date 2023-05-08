import pytest

from lingpy import Alignments

from lingrex.trimming import *


def test_Site():
    site = Site([GAP, 'a', GAP, 't'])
    assert site.gap_ratio() == pytest.approx(0.5)
    assert site.gap_ratio(gap='#') == pytest.approx(0.0)
    assert site.soundclass() == 'V'
    assert site.soundclass(gap='a') == '0'


@pytest.mark.parametrize(
    'alms,gap,ratios',
    [
        (["aaa", "aa-"], '-', [0.0, 0.0, 0.5]),
        (["aa", "a#"], '#', [0.0, 0.5]),
    ]
)
def test_gap_ratio(alms, gap, ratios):
    assert Sites([list(w) for w in alms], gap=gap).gap_ratios == ratios


def test_trimmed():
    alm = [list("toxta-"), list("to-tir"), list("to-t-r"), list("do--ar")]
    assert " ".join(Sites(alm)._trimmed([2, 5]).to_alignment()[0]) == "t o t a"


def test_soundclasses():
    assert Sites([list("-bc"), list("ab-")], gap="-").soundclasses == ["V", "C", "C"]


@pytest.mark.parametrize(
    'alms,kw,result',
    [
        (["abc", "a-c", "--c"], {}, 'ac'),
        (["abc", "a-c", "--c"], dict(skeletons=['VCC']), 'abc'),
        (["a+bco", "-+cco", "-+cco"], {}, 'bco'),
        (["a+b", "-+c", "-+c"], dict(exclude=""), 'a+b'),
    ]
)
def test_trim_by_gap(alms, kw, result):
    assert Sites([list(w) for w in alms]).trimmed(**kw).to_alignment()[0] == list(result)


@pytest.mark.parametrize(
    'alms,kw,result',
    [
        (["--mat", "-xmut", "--mit", "m-xit"], {}, 'mat'),
        (["--mat--", "-xmut--", "--mitx-", "m-xit-x"], {}, 'mat'),
    ]
)
def test_trim_by_core(alms, kw, result):
    sites = Sites([list(w) for w in alms])
    assert sites.trimmed(strategy='core', **kw).to_alignment()[0] == list(result)
    assert str(sites)


def test_trim_random(mocker):
    mocker.patch('lingrex.trimming.random', mocker.Mock(sample=lambda pop, k: list(pop)[:k]))
    alms = [list(w) for w in ["--mat", "-xmut", "m-xut", "--xit"]]
    assert len(Sites(alms).trimmed()) == len(Sites(alms).trimmed_random())
    assert set(Sites(alms).trimmed().soundclasses) == \
        set(Sites(alms).trimmed_random().soundclasses)
    assert Sites(alms).trimmed_random(strategy='core')


def test_prep_alignments(wl_with_alignments):
    test_wl = prep_alignments(Alignments(wl_with_alignments, transcription="form"))
    assert test_wl[4, "structure"] == "C V C V"
