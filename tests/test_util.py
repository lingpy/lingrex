import pytest
from lingrex.util import lingrex_path, add_structure
from lingpy import Wordlist, Alignments


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
