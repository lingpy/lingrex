"""
Test fuzzy reconstruction.
"""
import pytest
from lingrex.fuzzy import FuzzyReconstructor
from lingrex.reconstruct import (
        CorPaRClassifier,
        #OneHot,
        #ReconstructionBase,
        #PatternReconstructor,
        #transform_alignment,
        #eval_by_dist,
        #eval_by_bcubes
        )
import random
import lingpy

random.seed(1234)

def test_FuzzyReconstructor(data):

    pt = FuzzyReconstructor(str(data / "hillburmish.tsv"), target="ProtoBurmish", ref="cogids",
            fuzzy=False)
    alms = lingpy.align.sca.Alignments(
            str(data / "hillburmish.tsv"),
            transcription="form", ref="cogids")
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

    
