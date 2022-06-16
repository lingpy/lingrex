"""
Miscellaneous evaluation functions.
"""
import statistics
from lingpy.evaluate.acd import _get_bcubed_score as bcs
import lingpy


def compare_cognate_sets(wordlist, refA, refB):
    """
    Compute cognate set comparison statistics by computing B-Cubed Scores.

    Note
    ----
    This check was first described in Wu and List (forthcoming).

    > Wu, M.-S. and J.-M. List (2022): Annotating cognates in phylogenetic
    > studies of South-East Asian languages. Preprint: https://doi.org/10.17613/rabq-7z45
    """
    ranks = []
    for concept in wordlist.rows:
        cogsA = wordlist.get_list(row=concept, flat=True, entry=refA)
        cogsB = wordlist.get_list(row=concept, flat=True, entry=refB)
        p, r = bcs(cogsA, cogsB), bcs(cogsB, cogsA)
        f = 2 * (p * r) / (p + r)
        ranks += [[concept, p, r, f]]
    return ranks


def cross_semantic_cognate_statistics(
    wordlist,
    ref="cogids",
    concept="concept",
    morpheme_glosses="morphemes",
    ignore_affixes=True,
    affixes=("suf", "suffix", "SUF", "SUFFIX"),
):
    """
    Calculate colexification statistics for partial colexifications.

    :param wordlist: A LingPy wordlist.
    :param ref: Reference to the column with cognate identifiers.
    :param concept: Reference to the concept column.
    :param morpheme_glosses: Reference to the morpheme glosses.
    :param ignore_affixes: If set to True, will ignore morphemes flagged as affixes.
    :param affixes: List of strings that trigger that a morpheme gloss is
        ignored if it contains one of them as a substring.

    Note
    ----
    This check was first described in Wu and List (forthcoming).

    > Wu, M.-S. and J.-M. List (2022): Annotating cognates in phylogenetic
    > studies of South-East Asian languages. Preprint: https://doi.org/10.17613/rabq-7z45
    """

    # type check for basic types if they are not there
    for idx, cogids, morphemes in wordlist.iter_rows(ref, morpheme_glosses):
        wordlist[idx, ref] = lingpy.basictypes.ints(cogids)
        wordlist[idx, morpheme_glosses] = lingpy.basictypes.strings(morphemes)

    if ignore_affixes:
        D = {}
        for idx, cogids, morphemes in wordlist.iter_rows(ref, morpheme_glosses):
            new_cogids = []
            for cogid, morpheme in zip(cogids, morphemes):
                if not sum([1 if s in morpheme else 0 for s in affixes]):
                    new_cogids += [cogid]
            D[idx] = lingpy.basictypes.ints(new_cogids)
        wordlist.add_entries(ref + "_derived", D, lambda x: x)
        new_ref = ref + "_derived"
    else:
        new_ref = ref

    etd = wordlist.get_etymdict(ref=new_ref)
    indices = {ln: {} for ln in wordlist.cols}
    for i, ln in enumerate(wordlist.cols):
        for cogid, reflexes in etd.items():
            if reflexes[i]:
                concepts = [wordlist[idx, concept] for idx in reflexes[i]]
                indices[ln][cogid] = len(set(concepts)) - 1

    all_scores = []
    for cnc in wordlist.rows:
        # Loop through all the concepts in the data
        reflexes = wordlist.get_list(
            row=cnc, flat=True
        )  # The lexical entries of the concept.
        scores = []
        for idx in reflexes:
            doculect, cogids = wordlist[idx, "doculect"], wordlist[idx, new_ref]
            scores += [statistics.mean([indices[doculect][cogid] for cogid in cogids])]
        all_scores += [[cnc, statistics.mean(scores)]]
    return sorted(all_scores, key=lambda x: (x[1], x[0]))
