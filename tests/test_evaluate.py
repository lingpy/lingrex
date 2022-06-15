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

    wordlist = Wordlist({
            0: ["doculect", "concept", "form", "cogids", "morphemes"],
            1: ["a", "A", "a + b", "1 2", "a _suf"],
            2: ["b", "A", "a + b", "1 2", "a _suf"],
            3: ["c", "A", "c + d + a", "3 4 1", "_suf d a"],
            4: ["d", "A", "d + e", "4 5", "d e"],
            5: ["a", "B", "a + f", "1 6", "a f"],
            6: ["b", "B", "a + f", "1 6", "a f"],
            7: ["c", "C", "g + h + a", "7 8 1", "g h a"],
            8: ["d", "C", "h + i", "8 9", "h i"],
            })
    ranks = cross_semantic_cognate_statistics(
            wordlist,
            concept="concept",
            morpheme_glosses="morphemes",
            ignore_affixes=True
            )
    assert len(ranks) == 3
    assert ranks[0][0] == "C"
    assert ranks[2][1] == 0.625
    wordlist = Wordlist({
            0: ["doculect", "concept", "form", "cogids", "morphemes"],
            1: ["a", "A", "a + b", "1 2", "a _suf"],
            2: ["b", "A", "a + b", "1 2", "a _suf"],
            3: ["c", "A", "c + d + a", "3 4 1", "_suf d a"],
            4: ["d", "A", "d + e", "4 5", "d e"],
            5: ["a", "B", "a + f", "1 6", "a f"],
            6: ["b", "B", "a + f", "1 6", "a f"],
            7: ["c", "C", "g + h + a", "7 8 1", "g h a"],
            8: ["d", "C", "h + i", "8 9", "h i"],
            })
    ranks2 = cross_semantic_cognate_statistics(
            wordlist,
            concept="concept",
            morpheme_glosses="morphemes",
            ignore_affixes=False
            )
    assert ranks2[2][1] != ranks[2][1]
    assert ranks2[2][1] == 0.5



