from pytest import raises
from lingpy import Wordlist, Alignments
from lingrex.copar import CoPaR
from lingrex.util import add_structure
from lingrex.regularity import regularity


dummy_wl = {
    0: ["doculect", "concept", "form", "ipa", "alignment", "cogid"],
    1: ["A", "one", "atawu", "atawu", "a t a w u", 1],
    2: ["B", "one", "atwu", "atwu", "a t - w u", 1],
    3: ["C", "one", "tawu", "tawu", "- t a w u", 1],
    4: ["D", "one", "tefu", "tefu", "- t e f u", 1],
    5: ["A", "two", "satu", "satu", "s a t u", 2],
    6: ["B", "two", "setu", "setu", "s e t u", 2],
    7: ["C", "two", "situ", "situ", "s i t u", 2]
}


def test_regularity():
    test_wl = Wordlist(dummy_wl)
    with raises(ValueError):
        regularity(test_wl)

    test_alg = Alignments(test_wl)
    add_structure(test_alg, model="cv", structure="structure")
    test_alg = CoPaR(test_alg, ref="cogid")
    test_alg.get_sites()
    test_alg.cluster_sites()
    test_alg.sites_to_pattern()
    output = regularity(test_alg, threshold=2, word_threshold=0.5)

    assert output == (2, 5, 7, 0.29, 4, 5, 9, 0.44, 3, 4, 7, 0.43)
