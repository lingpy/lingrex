"""
Calculate regularity metrics on dataset.
"""
import statistics

from lingpy import log


def regularity(wordlist, threshold=3, ref="cogid", min_refs=3,
               word_threshold=0.75, sound_classes="cv"):
    """
    Check regularity in three flavors.

    - regularity based on the number of correspondence patterns that have more
      or the same number of sites as threshold
    - the proportion of correspondence patterns identified as regular via
      threshold counting all alignment sites
    - the proportion of words that we judge regular, judging words to be
      regular when more than the proportion word_threshold of sites are judged
      to be regular since they can be assigned to patterns that are covered by
      more than threshol sites

    :param wordlist: A lingpy Wordlist.
    :type wordlist: :class:lingpy.Wordlist
    :param threshold: The minimum number of alignment sites for a cognate set
        to be considered in the computation of regular words. Defaults to '3'.
    :type threshold: int
    :param ref: The column which stores the cognate sets, defaults to 'cogid'
    :type ref: str
    :param min_refs: The minimum number of occurrences a correspondence pattern
        to be considered recurring. Defaults to '3'.
    :type min_refs: int
    :param word_threshold: The relative threshold of patterns that need to be regular
        in order for a word to be considered regular as well. Defaults to '0.75'.
    :type word_threshold: float
    :param sound_classes: A string of characters or a list or a set of strings
        that contain the sound classes that the regularity should concentrate on.
    :type sound_clasess: str, list, set, tuple
    :return: Different scores of regularity.
    :rtype: tuple


    Note
    ----
    These regularity checks were first introduced in a study by Blum and List (2023):

    > Blum, F. and J.-M. List (2023): Trimming phonetic alignments improves the inference of
    > sound correspondence patterns from multilingual wordlists.
    > In: Proceedings of the 5th Workshop on Computational Typology and Multilingual NLP.
    > Association for Computational Linguistics 52-64. https://aclanthology.org/2023.sigtyp-1.6
    """
    if not hasattr(wordlist, "clusters"):
        raise ValueError("need a CoPaR object with clusters")
    patterns = len({p: len(vals) for p, vals in wordlist.clusters.items() \
            if p[0] in sound_classes})
    regular_patterns = len(
        [p for p, vals in wordlist.clusters.items() \
                if len(vals) >= threshold and p[0] in sound_classes])
    regular_proportion = sum(
        [len(vals) for p, vals in wordlist.clusters.items() \
                if len(vals) >= threshold and p[0] in sound_classes]
    )
    full_proportion = sum([len(vals) for p, vals in wordlist.clusters.items() \
            if p[0] in sound_classes])

    # get the proportion of words
    regular_words, irregular_words = 0, 0
    for cogid, msa in filter(
        lambda x: len(set(x[1]["taxa"])) >= min_refs, wordlist.msa[ref].items()
    ):
        scores = []
        for idx in range(len(msa["alignment"][0])):
            if (cogid, idx) not in wordlist.patterns:  # pragma: no cover
                log.warning("duplicate cognate in {0} / {1}".format(cogid, idx))
            else:
                if wordlist.patterns[cogid, idx][0][1] in sound_classes:
                    if (
                        max(
                            [
                                len(wordlist.clusters[b, c])
                                for a, b, c in wordlist.patterns[cogid, idx]
                            ]
                        )
                        >= threshold
                    ):
                        scores.append(1)
                    else:
                        scores.append(0)
        if scores:
            if statistics.mean(scores) >= word_threshold:
                regular_words += len(set(msa["taxa"]))
            else:
                irregular_words += len(set(msa["taxa"]))

    return (
        regular_patterns,
        patterns - regular_patterns,
        patterns,
        round((regular_patterns / patterns), 2),
        regular_proportion,
        full_proportion - regular_proportion,
        full_proportion,
        round((regular_proportion / full_proportion), 2),
        regular_words,
        irregular_words,
        regular_words + irregular_words,
        round((regular_words / (regular_words + irregular_words)), 2),
    )
