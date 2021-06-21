import pytest
from lingrex.cognates import common_morpheme_cognates, salient_cognates
from lingpy import Wordlist

def test_common_morpheme_cognates():

    D = {
        0: ["doculect", "concept", "ipa", "tokens", "cogids"],
        1: ["a", "a", "pla", "p l a + p u", [1, 2]],
        2: ["b", "a", "pla", "p l a t + k i", [1, 3]],
        3: ["c", "a", "pla", "k i + p l u p", [4, 1]],
        4: ["d", "a", "pla", "p l a k", [1]],
        5: ["a", "b", "pla", "t r a", [2]],
        6: ["b", "b", "pla", "t a t", [2]],
        7: ["c", "b", "pla", "d r ə p", [2]],
    }

    wl = Wordlist(D)
    common_morpheme_cognates(wl)
    assert wl[1, "autocogid"] == wl[2, "autocogid"]

def test_salient_cognates():

    D = {
        0: ["doculect", "concept", "morphemes", "tokens", "cogids",],
        1: ["a", "a", "pla _pi".split(), "p l a + p u".split(), [1, 2]],
        2: ["b", "a", "pla _po".split(), "p l a t + k i".split(), [1, 3]],
        3: ["c", "a", "_po pla".split(), "k i + p l u p".split(), [4, 1]],
        4: ["d", "a", "pla".split(), "p l a k".split(), [1]],
        5: ["a", "b", "pla".split(), "t r a".split(), [2]],
        6: ["b", "b", "pla".split(), "t a t".split(), [2]],
        7: ["c", "b", "pla".split(), "d r ə p".split(), [2]],
    }

    wl = Wordlist(D)
    salient_cognates(wl)
    assert wl[1, "newcogid"] == wl[2, "newcogid"]
