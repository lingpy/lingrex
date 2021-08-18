import pytest
from lingrex.borrowing import internal_cognates, external_cognates
from lingpy import Wordlist


@pytest.fixture
def wl(data):
    return Wordlist(str(data / 'wordlist.tsv'))


@pytest.mark.parametrize(
    'kw,success',
    [
        (
                dict(ref="autocogids", partial=True, method="lexstat"),
                lambda wl: "autocogids" in wl.columns),
        (
                dict(ref="autocogid", partial=False, method="lexstat"),
                lambda wl: "autocogid" in wl.columns),
        (
                dict(ref="autocogids", partial=True, method="sca"),
                lambda wl: "autocogids" in wl.columns),
    ]
)
def test_internal_cognates(kw, success, wl):
    internal_cognates(wl, runs=10, **kw)
    assert success(wl)
    
    etd = wl.get_etymdict(ref=kw["ref"])
    
    for cogid, vals in etd.items():
        concepts = []
        for idx_ in vals:
            if idx_:
                for idx in idx_:
                    concepts += [wl[idx, "concept"]]
        assert len(set(concepts)) == 1



def test_external_cognates(wl):
    external_cognates(wl, cognates="cogid", ref="borrids")
    assert "borrids" in wl.columns
