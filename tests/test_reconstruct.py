"""
Test the reconstruction module of lingrex.
"""
import pytest
from lingrex.reconstruct import (
        CorPaRClassifier,
        OneHot,
        ReconstructionBase,
        PatternReconstructor,
        transform_alignment,
        eval_by_dist,
        eval_by_bcubes
        )
from functools import partial



def test_transform_alignment():

    out = transform_alignment(
            [["b", "a", "k"], ["b", "a"]],
            ["a", "b"],
            ["a", "b", "u"],
            training=False
            )
    assert len(out) == 3

    out = transform_alignment(
            [["b", "k"], ["b", "a", "k"]],
            ["a", "b"],
            ["a", "b", "u"],
            training=True,

            )
    assert len(out) == 2

    out = transform_alignment(
            [["b", "k"], ["b", "a", "k"]],
            ["a", "b"],
            ["a", "b", "u"],
            training=True,
            firstlast=True

            )
    assert out[0][-1] == "k"

    out = transform_alignment(
            [["b", "k"], ["b", "a", "k"]],
            ["a", "b"],
            ["a", "b", "u"],
            training=True,
            startend=True
            )
    assert out[0][-1] == 0


def test_PatternReconstructor(data):

    pt = PatternReconstructor(str(data / "hillburmish.tsv"), "ProtoBurmish", ref="cogids",
            )
    t1 = partial(transform_alignment, align=True, position=False,
            prosody=False, startend=False, firstlast=False)
    t2 = partial(transform_alignment, align=True, position=True,
            prosody=True, startend=True, firstlast=True)
    pt.fit(func=t1)
    assert pt.predict(
            pt.msa["cogids"][665]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True
            ) == ['ŋ', 'a', '¹']
    pt.fit(func=t2)
    assert pt.predict(
            pt.msa["cogids"][665]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True
            ) == ['ŋ', 'a', '¹']

    pt.fit(func=t1, onehot=True)
    assert pt.predict(
            pt.msa["cogids"][665]["seqs"][:3],
            ["Atsi", "Lashi", "OldBurmese"],
            desegment=True
            ) == ['ŋ', 'a', '¹']

def test_eval_by_dist():
    assert eval_by_dist([[["t", "a"], ["t", "o"]]]) == 1
    assert eval_by_dist([[["t", "a"], []]]) == 2

    assert eval_by_dist([[["t", "a"], ["t", "o"]]], normalized=True) == 0.5

def test_eval_by_bcubes():
    assert eval_by_bcubes([[["t", "a"], ["t", "a"]]]) == 1
    assert eval_by_bcubes([
        [["t", "a"], ["t", "o"]]
        ]) == 1.0
    assert eval_by_bcubes([
        [["t", "a"], []]
        ]) == 1


