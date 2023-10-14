"""
Test fuzzy reconstruction.
"""
import pytest
from lingrex.fuzzy import FuzzyReconstructor, ntile
from lingrex.reconstruct import CorPaRClassifier
import random
import lingpy


def test_ntile():
    assert ntile(["kap", "kap", "kup", "kup"], 3) == 'k|k a|u p|p'
    # counting is not the same for missing data!
    assert ntile(["kap", "kØp", "kØp"], n=2) == 'k|k a|a p|p'
    

def test_FuzzyReconstructor(data):
    random.seed(1234)
    
    pytest.raises(ValueError, FuzzyReconstructor, 1, "ProtoBurmish")
    pt = FuzzyReconstructor(str(data / "hillburmish.tsv"), "ProtoBurmish", ref="cogids",
            fuzzy=False)
    alms = lingpy.align.sca.Alignments(
            str(data / "hillburmish.tsv"),
            transcription="form", ref="cogids")
    pt = FuzzyReconstructor(alms, "ProtoBurmish", ref="cogids",
            fuzzy=False)
    pt.random_splits()
    assert hasattr(pt, "wordlists")
    
    clf = lambda: CorPaRClassifier()
    pt.fit_samples(clf)
    predis = pt.predict(
            pt.wordlist.msa["cogids"][665]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True
            )
    assert predis[0] == "ŋ:100"
    predis = pt.predict(
            pt.wordlist.msa["cogids"][666]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True
            )
    assert predis[-1] == "?:90¦⁴:10"

    predis = pt.predict(
            pt.wordlist.msa["cogids"][665]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True,
            output="percentiles"
            )
    assert predis[0] == "ŋ:100"

    words, predis = pt.predict(
            pt.wordlist.msa["cogids"][665]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True,
            output="wp"
            )
    assert predis[0] == "ŋ:100"

    words = pt.predict(
            pt.wordlist.msa["cogids"][665]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True,
            output="words"
            )
    assert words[0][0] == "ŋ"




    
