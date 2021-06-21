import pytest
from lingrex.borrowing import internal_cognates, external_cognates
from lingpy import Wordlist
from pathlib import Path

def test_internal_cognates():

    wl = Wordlist(
            Path(__file__).parent.joinpath("data", "wordlist.tsv").as_posix(),
            )
    internal_cognates(wl, ref="autocogids", partial=True, method="lexstat",
            runs=10)
    assert "autocogids" in wl.columns
    wl = Wordlist(
            Path(__file__).parent.joinpath("data", "wordlist.tsv").as_posix(),
            )
    internal_cognates(wl, ref="autocogid", partial=False, method="lexstat",
            runs=10)
    assert "autocogid" in wl.columns
    wl = Wordlist(
            Path(__file__).parent.joinpath("data", "wordlist.tsv").as_posix(),
            )

    internal_cognates(wl, ref="autocogids", partial=True, method="sca",
            runs=10)
    assert "autocogids" in wl.columns


def test_external_cognates():
    wl = Wordlist(
            Path(__file__).parent.joinpath("data", "wordlist.tsv").as_posix(),
            )
    external_cognates(wl, cognates="cogid", ref="borrids")
    assert "borrids" in wl.columns


