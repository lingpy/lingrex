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
        (["abc", "a-c", "--c"], {}, list('ac')),
        (["abc", "a-c", "--c"], dict(skeletons=['VCC']), list('abc')),
        (["a+bco", "-+cco", "-+cco"], {}, list('bco')),
        (["a+b", "-+c", "-+c"], dict(exclude=""), list('a+b')),
        ([
             #"- - n u - - 'b/b a".split(),
             '- - - - d ù/u - -'.split(),
             '- - - - d ú/u - -'.split(),
             '- - - - d ù/u - -'.split(),
             "ɾ u 'w/w a s i ɾ a".split(),
             '- - - - s u - e'.split(),
             "- - n u - - 'b/b a".split(),
             '- - - - d u l -'.split(),
             '- - n u k - w ɔ'.split(),
         ], {}, ['d', 'ù/u']),
        ([
             "- - n u - - 'b/b a".split(),
             '- - - - d ù/u - -'.split(),
             '- - - - d ú/u - -'.split(),
             '- - - - d ù/u - -'.split(),
             "ɾ u 'w/w a s i ɾ a".split(),
             '- - - - s u - e'.split(),
             #"- - n u - - 'b/b a".split(),
             '- - - - d u l -'.split(),
             '- - n u k - w ɔ'.split(),
         ], {}, ['-', '-']),
        # Non-overlapping alignments:
        (['- - a b'.split(), 'a b - -'.split(), 'a b - -'.split()], {}, ['-', '-']),
        #
        (['- a b'.split(), 'b a -'.split(), 'b a -'.split()], {}, ['-', 'a']),
        (['- a b'.split(), 'b - -'.split(), 'b - -'.split()], {}, ['-', 'a']),
        ([
             '- a b c'.split(),
             'b - - -'.split(),
             'b - - -'.split(),
             'b - - d'.split()
         ], {}, ['-', 'a', 'c']),
        ([list('bbabb'), list('bb-bb'), list('-b-b-'), list('-b-b-')], {}, list('bbabb')),
        ([list('bbabb'), list('bb-bb'), list('-b-b-'), list('-b-b-')], {'strict_ratio': False}, list('bab')),
    ]
)
def test_trim_by_gap(alms, kw, result):
    assert Sites([list(w) for w in alms]).trimmed(**kw).to_alignment()[0] == result


@pytest.mark.parametrize(
    'alms,kw,result',
    [
        (["--mat", "-xmut", "--mit", "m-xit"], {}, list('mat')),
        (["--mat--", "-xmut--", "--mitx-", "m-xit-x"], {}, list('mat')),
        ([
            "- - n u - - 'b/b a".split(),
            '- - - - d ù/u - -'.split(),
            '- - - - d ú/u - -'.split(),
            '- - - - d ù/u - -'.split(),
            "ɾ u 'w/w a s i ɾ a".split(),
            '- - - - s u - e'.split(),
            "- - n u - - 'b/b a".split(),
            '- - - - d u l -'.split(),
            '- - n u k - w ɔ'.split(),
         ], {}, ['-', '-', "'b/b", 'a']),
        ([list('bbabb'), list('bb-bb'), list('-b-b-'), list('-b-b-')], {}, list('bab')),
    ]
)
def test_trim_by_core(alms, kw, result):
    sites = Sites([list(w) for w in alms])
    assert sites.trimmed(strategy='core', **kw).to_alignment()[0] == result
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
