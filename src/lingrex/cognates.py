"""
Operations with cognate sets.
"""
import collections

from clldutils.text import strip_brackets, split_text
import lingpy


def common_morpheme_cognates(
    wordlist, cognates="cogids", ref="autoid", morphemes="automorphemes", override=True
):
    """
    Convert partial cognates to full cognates.
    """

    C, M = {}, {}
    current = 1
    for concept in wordlist.rows:
        base = split_text(strip_brackets(concept))[0].upper().replace(" ", "_")
        idxs = wordlist.get_list(row=concept, flat=True)
        cogids = collections.defaultdict(list)
        for idx in idxs:
            M[idx] = [c for c in wordlist[idx, cognates]]
            for cogid in lingpy.basictypes.ints(wordlist[idx, cognates]):
                cogids[cogid] += [idx]
        for i, (cogid, idxs) in enumerate(
            sorted(cogids.items(), key=lambda x: len(x[1]), reverse=True)
        ):
            for idx in idxs:
                if idx not in C:
                    C[idx] = current
                    M[idx][M[idx].index(cogid)] = base
                else:
                    M[idx][M[idx].index(cogid)] = "_" + base.lower()
            current += 1
    wordlist.add_entries(ref, C, lambda x: x)
    if morphemes:
        wordlist.add_entries(morphemes, M, lambda x: x, override=override)


def salient_cognates(
    wordlist, cognates="cogids", ref="newcogid", morphemes="morphemes", override=True
):
    """
    Convert partial cognates to full cognates ignoring non-salient cognate sets.
    """

    lookup, D = {}, {}
    for idx, cogids, morphemes in wordlist.iter_rows(cognates, morphemes):
        selected_cogids = []
        for cogid, morpheme in zip(cogids, morphemes):
            if not morpheme.startswith("_"):
                selected_cogids += [cogid]
        salient = tuple(selected_cogids)
        if salient in lookup:
            D[idx] = lookup[salient]
        elif D.values():
            next_cogid = max(D.values()) + 1
            lookup[salient] = next_cogid
            D[idx] = next_cogid
        else:
            lookup[salient] = 1
            D[idx] = 1

    wordlist.add_entries(ref, D, lambda x: x, override=override)
