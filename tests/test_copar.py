import pytest
from lingrex.copar import (
        CoPaR, consensus_pattern,
        incompatible_columns,
        score_patterns,
        density
        )
from lingpy import Wordlist, Alignments
from lingrex.util import add_structure


def test_consensus_pattern():
    
    missing = "?"
    out = consensus_pattern(
            [
                ["a", "b", "c"],
                ["a", "b", missing],
                [missing, missing, "c"]
                ], missing=missing)
    assert out == ("a", "b", "c")
    with pytest.raises(ValueError):
        consensus_pattern([
            ["a", "b"],
            ["a", "c"]])

def test_incompatible_columns():

    missing = "?"
    out = incompatible_columns(
            [
                ["a", "b", "c"],
                ["I", "b", "c"],
                ["a", "b", missing],
                [missing, missing, "c"]
                ], missing=missing)
    assert out[0] == "*"


def test_score_patterns():
    
    missing = "?"
    out = score_patterns(
            [
                ["a", "b", "c"],
                ["a", "b", "d"],
                ["a", "b", missing],
                [missing, missing, "c"]
                ], 
            missing=missing, 
            mode="coverage"
            )
    assert out == -1

    out = score_patterns(
            [
                ["a", "b", "c"],
                ["a", "b", "d"],
                ["a", "b", missing],
                [missing, missing, "c"]
                ], 
            missing=missing, 
            mode="ranked"
            )
    assert score_patterns(["a", "b", "c"]) == -1
    assert out == -1


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
    assert round(density(Wordlist(D), ref="cogid"), 2) == 0.67 


def test_CoPaR():
    D = {
            0: ["doculect", "concept", "ipa", "tokens", "cogids", "alignment"],
            1: ["a", "a", "pla", "p l a", [1], "p l a -".split()],
            2: ["b", "a", "pla", "p l a t", [1], "p l a t".split()],
            3: ["c", "a", "pla", "p l u p", [1], "p l u p".split()],
            4: ["d", "a", "pla", "p l a k", [1], "p l a k".split()],
            5: ["a", "b", "pla", "t r a", [2], "t r a -".split()],
            6: ["b", "b", "pla", "t a t", [2], "t - a t".split()],
            7: ["c", "b", "pla", "d r ə p", [2], "d r ə p".split()],
            #8: ["a", "b", "pla", "p l a k", [1], "d x a k".split()],
            9: ["a", "c", "pla", "k l a", [3], "k r a -".split()],
            #10: ["a", "c", "pla", "p l a t", [1], "k  a t".split()],
            11: ["c", "c", "pla", "k l ə p", [3], "k l ə p".split()],
            12: ["d", "c", "pla", "g l a k", [3], "g l a k".split()],
            }
    alms = Alignments(D, ref="cogids", transcription="ipa")
    add_structure(alms, model="cv", structure="structure")
    cop = CoPaR(alms, ref="cogids", structure="structure", minrefs=1)
    cop.get_sites()
    assert len(cop.sites) == 12
    cop.cluster_sites()
    assert len(cop.clusters) == 9
    cop.sites_to_pattern()
    cop.add_patterns()
    cop.refine_patterns()
    cop.irregular_patterns()

    cop.fuzziness()



if __name__ == "__main__":
    missing = "?"
    out = score_patterns(
            [
                ["a", "b", "c"],
                ["a", "b", "d"],
                ["a", "b", missing],
                [missing, missing, "c"]
                ], 
            missing=missing, 
            mode="ranked"
            )
    print(out)

    D = {
            0: ["doculect", "concept", "ipa", "tokens", "cogids", "alignment"],
            1: ["a", "a", "pla", "p l a", [1], "p l a -".split()],
            2: ["b", "a", "pla", "p l a t", [1], "p l a t".split()],
            3: ["c", "a", "pla", "p l u p", [1], "p l u p".split()],
            4: ["d", "a", "pla", "p l a k", [1], "F/p l a k".split()],
            5: ["a", "b", "pla", "t r a", [2], "t r a -".split()],
            6: ["b", "b", "pla", "t a t", [2], "t - a t".split()],
            7: ["c", "b", "pla", "d r ə p", [2], "d r ə p".split()],
            #8: ["a", "b", "pla", "p l a k", [1], "d x a k".split()],
            9: ["a", "c", "pla", "k l a", [3], "k l a -".split()],
            #10: ["a", "c", "pla", "p l a t", [1], "k  a t".split()],
            11: ["c", "c", "pla", "k l ə p", [3], "k l ə p".split()],
            12: ["d", "c", "pla", "g l a k", [3], "g l a k".split()],
            }
    alms = Alignments(D, ref="cogids", transcription="ipa")
    add_structure(alms, model="cv", structure="structure")
    cop = CoPaR(alms, ref="cogids", structure="structure", minrefs=1)
    cop._check()
    print(cop._str_type)
    cop.get_sites()
    print(len(cop.sites))
    cop.cluster_sites()
    print(cop.clusters)
    print(len(cop.clusters))
    cop.sites_to_pattern()
    cop.refine_patterns()
    cop.irregular_patterns()
    print(len(cop.clusters))


