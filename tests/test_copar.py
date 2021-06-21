import pytest
from lingrex.copar import (
    CoPaR,
    consensus_pattern,
    incompatible_columns,
    score_patterns,
    density,
)
from lingpy import Wordlist, Alignments
from lingrex.util import add_structure


def test_consensus_pattern():
    missing = "?"
    out = consensus_pattern(
        [["a", "b", "c"], ["a", "b", missing], [missing, missing, "c"]], missing=missing
    )
    assert out == ("a", "b", "c")
    with pytest.raises(ValueError):
        consensus_pattern([["a", "b"], ["a", "c"]])


def test_incompatible_columns():
    missing = "?"
    out = incompatible_columns(
        [
            ["a", "b", "c"],
            ["I", "b", "c"],
            ["a", "b", missing],
            [missing, missing, "c"],
        ],
        missing=missing,
    )
    assert out[0] == "*"


@pytest.mark.parametrize(
    'patterns,mode,result',
    [
        ([["a", "b", "c"], ["a", "b", "d"], ["a", "b", "?"], ["?", "?", "c"]], 'coverage', -1),
        (["a", "b", "c"], 'coverage', -1),
        ([["a", "b", "c"], ["a", "b", "c"], ["a", "b", "?"], ["?", "?", "c"]], 'ranked', 0.75),
        ([["a", "b", "c"], ["a", "b", "c"], ["a", "b", "?"], ["?", "?", "c"]], 'squared', 0.64),
        ([["a", "b", "c"], ["a", "b", "c"], ["a", "b", "?"], ["?", "?", "c"]], 'pairs', 0.44),
        ([["a", "b", "c"], ["a", "b", "c"], ["a", "b", "?"], ["?", "?", "c"]], 'coverage', 0.75),
    ]
)
def test_score_patterns(patterns, mode, result):
    assert result == pytest.approx(score_patterns(patterns, missing='?', mode=mode), abs=1e-2)


def test_score_patterns_error():
    with pytest.raises(ValueError):
        score_patterns([["a", "b"], ["a", "b"]], mode="bla")


def test_density():
    D = {
        0: ["doculect", "concept", "tokens", "ipa", "cogid"],
        1: ["a", "b", "t o x t ə".split(), "tochter", 1],
        2: ["b", "b", "t o x t ə".split(), "tochter", 1],
        3: ["c", "b", "t o x t ə".split(), "tochter", 1],
        4: ["a", "c", "t o x t ə".split(), "tochter", 2],
        5: ["b", "c", "t o x t ə".split(), "tochter", 2],
        6: ["c", "c", "t o x t ə".split(), "tochter", 2],
    }
    assert 0.67 == pytest.approx(density(Wordlist(D), ref="cogid"), abs=1e-2)


def test_CoPaR_fuzzy():
    D = {
        0: ["doculect", "concept", "ipa", "tokens", "cogids", "alignment"],
        1: ["a", "a", "pla", "p l a", [1], "p l a -".split()],
        2: ["b", "a", "pla", "p l a t", [1], "p l a t".split()],
        3: ["c", "a", "pla", "p l u p", [1], "p l u p".split()],
        4: ["d", "a", "pla", "p l a k", [1], "p l a k".split()],
        5: ["a", "b", "pla", "t r a", [2], "t r a -".split()],
        6: ["b", "b", "pla", "t a t", [2], "t - a t".split()],
        7: ["c", "b", "pla", "d r ə p", [2], "d r ə p".split()],
        # 8: ["a", "b", "pla", "p l a k", [1], "d x a k".split()],
        9: ["a", "c", "pla", "k l a", [3], "k r a -".split()],
        # 10: ["a", "c", "pla", "p l a t", [1], "k  a t".split()],
        11: ["c", "c", "pla", "k l ə p", [3], "k l ə p".split()],
        12: ["d", "c", "pla", "g l a k", [3], "g l a k".split()],
        13: ["d", "f", "buk", "b u k", [4], "b u k".split()],
    }
    alms = Alignments(D, ref="cogids", transcription="ipa")
    with pytest.raises(ValueError):
        CoPaR(alms, ref="cogids", structure="structure", minrefs=2)
    add_structure(alms, model="cv", structure="structure")
    cop = CoPaR(alms, ref="cogids", structure="structure", minrefs=1)
    cop.get_sites()
    assert len(cop.sites) == 12
    cop.cluster_sites()
    assert len(cop.clusters) == 9
    cop.sites_to_pattern()
    cop.add_patterns()
    cop.irregular_patterns()
    cop.fuzziness()
    # get the cluster graph
    G = cop.get_cluster_graph()
    assert len(G.nodes) == len(cop.sites)

    # compute the purity of the cluster graph
    assert round(cop.purity(), 2) == 0.42
    cop.load_patterns()


def test_CoPaR_plain(tmp_path):
    D = {
        0: ["doculect", "concept", "ipa", "tokens", "cogid", "alignment"],
        1: ["a", "a", "pla", "p l a", 1, "p l a -".split()],
        2: ["b", "a", "pla", "p l a t", 1, "p l a t".split()],
        3: ["c", "a", "pla", "p l u p", 1, "p l u p".split()],
        4: ["d", "a", "pla", "p l a k", 1, "p l a k".split()],
        5: ["a", "b", "pla", "t r a", 2, "t r a -".split()],
        6: ["b", "b", "pla", "t a t", 2, "t - a t".split()],
        7: ["c", "b", "pla", "d r ə p", 2, "d r ə p".split()],
        # 8: ["a", "b", "pla", "p l a k", [1], "d x a k".split()],
        9: ["a", "c", "pla", "k l a", 3, "k r a -".split()],
        # 10: ["a", "c", "pla", "p l a t", [1], "k  a t".split()],
        11: ["c", "c", "pla", "k l ə p", 3, "k l ə p".split()],
        12: ["d", "c", "pla", "g l a k", 3, "g l a k".split()],
        13: ["d", "f", "buk", "b u k", 4, "b u k".split()],
    }
    alms = Alignments(D, ref="cogid", transcription="ipa")
    add_structure(alms, model="cv", structure="structure")
    cop = CoPaR(alms, ref="cogid", structure="structure", minrefs=1)

    with pytest.raises(ValueError):
        cop.write_patterns("f")
    with pytest.raises(ValueError):
        cop.predict_words()

    cop.get_sites()
    assert len(cop.sites) == 12
    cop.cluster_sites()
    assert len(cop.clusters) == 9
    cop.sites_to_pattern()
    cop.irregular_patterns()
    cop.add_patterns(proto="a", irregular_patterns=True)
    cop.fuzziness()
    # get the cluster graph
    G = cop.get_cluster_graph()
    assert len(G.nodes) == len(cop.sites)

    # compute the purity of the cluster graph
    assert round(cop.purity(), 2) == 0.42

    assert cop.upper_bound() > 1

    cop.write_patterns(tmp_path / 'test')
    cop.write_patterns(tmp_path / 'test', proto="a", irregular_patterns=True)
    cop.predict_words()
    cop.load_patterns()


def test_polynesian(data):
    cop = CoPaR(str(data / "east-polynesian.tsv"), ref="cogid", segments="segments")
    cop.align()
    cop.get_sites()
    cop.cluster_sites()
    cop.add_patterns()


def test_warnings():
    D = {
        0: ["doculect", "concept", "ipa", "tokens", "cogids", "alignment", "structure"],
        1: ["a", "a", "pla", "p l a", [1], "p l a -".split(), "i m n c".split()],
        2: ["b", "a", "pla", "p l a t", [1], "p l a t".split(), "i n c".split()],
        3: ["c", "a", "pla", "p l u p", [1], "p l u p".split(), "i m n c".split()],
        4: ["d", "a", "pla", "p l a k", [1], "p l a k".split(), "i m n c".split()],
    }
    alms = Alignments(D, ref="cogids")
    cop = CoPaR(alms, structure="structure")
    with pytest.raises(ValueError):
        cop.get_sites()
    D = {
        0: ["doculect", "concept", "ipa", "tokens", "cogids", "alignment", "structure"],
        1: ["a", "a", "pla", "p l a", [1], "p !l a -".split(), "i m n".split()],
        2: ["b", "a", "pla", "p l a t", [1], "p f/l a t".split(), "i m n c".split()],
        3: ["c", "a", "pla", "p l u p", [1], "p l u p".split(), "i m n c".split()],
        4: ["d", "a", "pla", "p l a k", [1], "p l a k".split(), "i m n c".split()],
    }
    cop = CoPaR(D)
    cop.get_sites()
