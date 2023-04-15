from lingrex.trimming import (
        apply_trim, revert, get_skeleton, consecutive_gaps, gap_profile,
        trim_by_gap,trim_by_core, trim_random, subsequence_of
        )
import random
random.seed(1234)


def test_gap_profile():
    assert gap_profile(["aaa", "aaa"], gap="-") == [0.0, 0.0, 0.0]
    assert gap_profile(["aa", "a-"], gap="-") == [0.0, 0.5]

def test_apply_trim():

    alm = [
            list("toxta-"), list("to-tir"), 
            list("to-t-r"), list("do--ar")]
    trimmed = apply_trim(alm, [2, 5])
    assert " ".join(trimmed[0]) == "t o t a" 


def test_get_skeleton():
    assert get_skeleton(["-bc", "ab-"], gap="-") == ["V", "C", "C"]


def test_revert():
    assert revert(["ab", "ab"]) == [["a", "a"], ["b", "b"]]


def test_trim_by_gap():
    assert trim_by_gap(["abc", "a-c", "--c"]) == [1]
    assert trim_by_gap(["a+bco", "-+cco", "-+cco"], exclude="_+") == [0, 1]
    assert trim_by_gap(["a+b", "-+c", "-+c"], exclude="") == []


def test_subsequence_of():
    assert subsequence_of("cvc", "cvcvc")
    assert subsequence_of("bla", "bla")
    assert not subsequence_of("abc", "ab")


def test_consecutive_gaps():
    left, right = consecutive_gaps("--mat-gu-go--", gap="-")
    assert " ".join([c for i, c in enumerate("--mat-gu-go--") if i not in
                     left+right]
                    ) == "m a t - g u - g o"


def test_trim_by_core():
    assert trim_by_core(["--mat", "-xmut", "--mit", "m-xit"]) == [0, 1]
    assert trim_by_core(["--mat--", "-xmut--", "--mitx-", "m-xit-x"]) == [0, 1, 5, 6]



def test_trim_random():
    
    alms = ["--mat", "-xmut", "m-xut", "--xit"]
    assert len(trim_by_gap(alms)) == len(trim_random(alms))
    assert set(get_skeleton(apply_trim(alms, trim_by_gap(alms)))) == \
            set(get_skeleton(apply_trim(alms, trim_random(alms))))

