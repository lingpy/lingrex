"""
Test evaluate module of lingrex.
"""
from lingrex.evaluate import (
        compare_cognate_sets,
        cross_semantic_cognate_statistics
        )
from lingpy import Wordlist


def test_compare_cognate_sets():

    wordlist = Wordlist({
            0: ["doculect", "concept", "form", "looseid", "strictid"],
            1: ["a", "a", "b", "1", "2"],
            2: ["b", "a", "c", "1", "3"],
            3: ["c", "a", "c", "1", "2"],
            4: ["d", "a", "d", "1", "4"]
            })
    ranks = compare_cognate_sets(
            wordlist, "strictid", "looseid")
    assert len(ranks) == 1
    assert ranks[0][0] == "a"
    assert ranks[0][1] == 1
    assert ranks[0][2] == 0.375 


def test_cross_semantic_cognate_statistics():

    pass
