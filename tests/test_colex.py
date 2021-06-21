import pytest
from lingrex.colex import (
    expand_alignment,
    find_bad_internal_alignments,
    compatible,
    merge_alignments,
    find_colexified_alignments,
)
from lingpy import Alignments


def test_find_bad_internal_alignments():

    wl = Alignments(
        {
            0: ["doculect", "concept", "ipa", "tokens", "alignment", "cogids"],
            1: ["a", "a", "bla", "b l a", "b l a -".split(), [1]],
            2: ["b", "a", "bla", "b l a k", "b l a k".split(), [1]],
            3: ["c", "a", "bla", "b a k", "b - a k".split(), [1]],
            4: ["a", "b", "bla", "b l a k", "b l a k".split(), [1]],
            5: ["b", "b", "bla", "b l a k", "b l a k".split(), [1]],
            6: ["a", "c", "bla", "b l a", "b l a -".split(), [1]],
        },
        ref="cogids",
    )
    find_bad_internal_alignments(wl)
    assert wl[4, "cogids"][0] != 1


def test_expand_alignment():

    missing = "?"
    out = expand_alignment(
        {"taxa": ["a", "b", "c"], "alignment": [["t", "a"], ["t/p", "u"], ["t", "-"]]},
        ["a", "d", "b", "c"],
        missing=missing,
    )
    assert out[1][1] == missing
    assert out[2][0] == "p"


def test_compatible():
    missing = "?"
    matches = compatible(
        [["a", "b"], [missing, missing], ["a", "c"], ["a", "d"]],
        [
            ["a", "-", "b"],
            ["a", "x", "b"],
            ["a", "-", "c"],
            [missing, missing, missing],
        ],
        missing=missing,
    )
    assert matches == 2

    matches = compatible(
        [["a", "b"], [missing, missing], ["a", "c"], ["a", "d"]],
        [
            ["a", "-", "c"],
            ["a", "x", "b"],
            ["a", "-", "c"],
            [missing, missing, missing],
        ],
        missing=missing,
    )
    assert not matches


def test_merge_alignments():
    missing = "?"
    matches = merge_alignments(
        [
            ["-", "a", "b"],
            [missing, missing, missing],
            ["-", "a", "c"],
            ["x", "a", "d"],
        ],
        [
            ["a", "-", "b"],
            ["a", "x", "b"],
            ["a", "-", "c"],
            [missing, missing, missing],
        ],
        missing=missing,
    )
    assert len(matches[0]) == 4

    missing = "?"
    matches = merge_alignments(
        [["a", "b"], ["a", "c"], ["a", "d"]],
        [
            ["a", "-", "b"],
            ["a", "x", "b"],
            ["a", "-", "c"],
        ],
        missing=missing,
    )
    assert len(matches[0]) == 3

    missing = "?"
    matches = merge_alignments(
        [["a", "b"], ["a", "c"], [missing, missing], ["a", "d"]],
        [
            ["a", "-", "b", "-"],
            ["a", "x", "b", "-"],
            ["a", "-", "c", "e"],
            [missing, missing, missing, missing],
        ],
        missing=missing,
    )
    assert len(matches[0]) == 4


def test_find_colexified_alignments():
    wl = Alignments(
        {
            0: ["doculect", "concept", "ipa", "tokens", "alignment", "cogids"],
            1: ["a", "a", "bla", "b l a", "b l a -".split(), [1]],
            2: ["b", "a", "bla", "b l a k", "b l a k".split(), [1]],
            3: ["c", "a", "bla", "b a k", "b - a k".split(), [1]],
            4: ["a", "b", "bla", "b l a k", "b l a -".split(), [2]],
            5: ["b", "b", "bla", "b l a k", "b l a k".split(), [2]],
            6: ["a", "c", "bla", "b l a", "- b l a".split(), [3]],
            7: ["d", "c", "bla", "a b l", "a b l -".split(), [3]],
        },
        ref="cogids",
    )

    find_colexified_alignments(wl)
    assert wl[1, "crossids"][0] == 1

    wl = Alignments(
        {
            0: ["doculect", "concept", "ipa", "tokens", "alignment", "cogid"],
            1: ["a", "a", "bla", "b l a", "b l a -".split(), 1],
            2: ["b", "a", "bla", "b l a k", "b l a k".split(), 1],
            3: ["c", "a", "bla", "b a k", "b - a k".split(), 1],
            4: ["a", "b", "bla", "b l a k", "b l a -".split(), 2],
            5: ["b", "b", "bla", "b l a k", "b l a k".split(), 2],
            6: ["a", "c", "bla", "b l a", "- b l a".split(), 3],
            7: ["d", "c", "bla", "a b l", "a b l -".split(), 3],
        },
        ref="cogid",
    )

    find_colexified_alignments(wl, cognates="cogid", ref="crossid")
    assert wl[1, "crossid"] == 1


#if __name__ == "__main__":
#
#    wl = Alignments(
#        {
#            0: ["doculect", "concept", "ipa", "tokens", "alignment", "cogids"],
#            1: ["a", "a", "bla", "b l a", "b l a -".split(), [1]],
#            2: ["b", "a", "bla", "b l a k", "b l a k".split(), [1]],
#            3: ["c", "a", "bla", "b a k", "b - a k".split(), [1]],
#            4: ["a", "b", "bla", "b l a k", "b l a k".split(), [1]],
#            5: ["b", "b", "bla", "b l a k", "b l a k".split(), [1]],
#            6: ["a", "c", "bla", "b l a", "b l a -".split(), [1]],
#        },
#        ref="cogids",
#    )
#    find_bad_internal_alignments(wl)
#    for idx in wl:
#        print(idx, wl[idx, "doculect"], wl[idx, "cogids"][0], wl[idx, "alignment"])
#
#    missing = "Ø"
#    matches = merge_alignments(
#        [
#            ["-", "a", "b"],
#            [missing, missing, missing],
#            ["-", "a", "c"],
#            ["x", "a", "d"],
#        ],
#        [
#            ["a", "-", "b"],
#            ["a", "x", "b"],
#            ["a", "-", "c"],
#            [missing, missing, missing],
#        ],
#        missing=missing,
#    )
#    for row in matches:
#        print(" ".join(row))
#
#    wl = Alignments(
#        {
#            0: ["doculect", "concept", "ipa", "tokens", "alignment", "cogids"],
#            1: ["a", "a", "bla", "b l a", "b l a -".split(), [1]],
#            2: ["b", "a", "bla", "b l a k", "b l a k".split(), [1]],
#            3: ["c", "a", "bla", "b a k", "b - a k".split(), [1]],
#            4: ["a", "b", "bla", "b l a k", "b l a -".split(), [2]],
#            5: ["b", "b", "bla", "b l a k", "b l a k".split(), [2]],
#            6: ["a", "c", "bla", "b l a", "- b l a".split(), [3]],
#            7: ["d", "c", "bla", "a b l", "a b l -".split(), [3]],
#        },
#        ref="cogids",
#    )
#
#    find_colexified_alignments(wl)
#    for idx in wl:
#        print(
#            idx,
#            "\t",
#            wl[idx, "doculect"],
#            "\t",
#            wl[idx, "alignment"],
#            "\t",
#            wl[idx, "cogids"][0],
#            "\t",
#            wl[idx, "crossids"][0],
#        )
