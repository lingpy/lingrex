"""
Various phonetic alignment functions.
"""
from lingpy import basictypes as bt


def gap_free_pairwise(seqA, seqB, syllables=None, gap="-"):
    """
    Carry out a gap-free alignment in which segments are merged instead of gapped.
    """
    syllables = [] if syllables is None else syllables
    start = True
    merge = False
    outA, outB = [], []
    for i, (charA, charB) in enumerate(zip(seqA, seqB)):
        if i in syllables:
            start = True
        if start and charB == gap:
            outA.append(charA + ">")
            merge = True
        elif not merge and charB == gap:
            outA[-1] += "<" + charA
        elif merge:
            if charB == gap:
                outA[-1] += charA + ">"
            else:
                outA[-1] += charA
                outB.append(charB)
                merge = False
        else:
            outA.append(charA)
            outB.append(charB)
        start = False
    return outA, outB


def align_to_template(sequence, structures, template, gap="-"):
    """
    Align a sequence to a template.
    """
    if (len(sequence) != len(structures)) or (len(template) < len(sequence)):
        raise ValueError(
            "sequence {0} and structure {1} have different length".format(
                repr(sequence), repr(structures)
            )
        )
    if len([x for x in structures if x not in template]) != 0:
        raise ValueError(
            "{0} items in the structure {1} is not in the template".format(
                len([x for x in structures if x not in template]), repr(structures)
            )
        )

    out = []
    idxA, idxB = 0, 0
    while idxB < len(template):
        if idxA < len(sequence):
            segment, structure = sequence[idxA], structures[idxA]
        else:
            segment, structure = gap, ""
        current_structure = template[idxB]
        if current_structure == structure:
            out.append(segment)
            idxA += 1
        else:
            out.append(gap)
        idxB += 1

    return out


def shrink_alignments(alignments):
    """
    Remove columns from alignment which all consist of gaps.
    """
    excludes = []
    for i in range(len(alignments[0])):
        col = set([line[i] for line in alignments])
        if "-" in col and len(col) == 1:
            excludes.append(i)
    return [
        [site for i, site in enumerate(alignment) if i not in excludes]
        for alignment in alignments
    ]


def shrink(tokens, structures, converter):
    """
    Shrink tokens according to the converter.

    .. note:: Works only for shrinking two structure elements so far.
    """
    outt, outs = [], []
    sm, merge = None, False
    for i in range(len(tokens)):
        if i > 0:
            sm = " ".join([structures[i - 1], structures[i]])
            if sm in converter:
                outt += [tokens[i - 1] + tokens[i]]
                outs += [converter[sm]]
                merge = True
            elif not merge:
                outt += [tokens[i - 1]]
                outs += [converter.get(structures[i - 1], structures[i - 1])]
            else:
                merge = False
    if sm not in converter:
        outt += [tokens[i]]
        outs += [converter.get(structures[i], structures[i])]
    return outt, outs


def shrink_template(
    wordlist,
    structure="structure",
    segments="tokens",
    converter={"i m": "I", "i": "I", "n c": "R", "n": "R", "c": "R"},
    new_structure="structure2",
    new_tokens="tokens2",
    override=False,
):
    """
    Reduce a template by merging certain parts of the structure.
    """
    D = {}
    for idx, strucs, tokens in wordlist.iter_rows(structure, segments):
        D[idx] = shrink(tokens, strucs, converter)
    wordlist.add_entries(new_structure, D, lambda x: bt.lists(x[1]), override=override)
    wordlist.add_entries(new_tokens, D, lambda x: bt.lists(x[0]), override=override)


def template_alignment(
    wordlist,
    ref="cogid",
    template="CCCCVVccccT_CCCCVVccccT_CCCCVVccccT_CCCCVVccccT_CCCCvvT",
    structure="structure",
    fuzzy=False,
    segments="tokens",
    gap="-",
    alignment="alignment",
    override=True,
):
    """
    Function aligns the cognate sets in a wordlist to a template.
    """

    for idx, tokens, structures in wordlist.iter_rows(segments, structure):
        wordlist[idx, segments], wordlist[idx, structure] = bt.lists(tokens), bt.lists(
            structures
        )

    etd = wordlist.get_etymdict(ref)
    A = {}
    if not fuzzy:
        for cogid, vals in etd.items():
            idxs = []
            for val in vals:
                if val:
                    idxs += val
            alignments = shrink_alignments(
                [
                    align_to_template(
                        wordlist[idx, segments],
                        wordlist[idx, structure],
                        template,
                        gap=gap,
                    )
                    for idx in idxs
                ]
            )
            for idx, alm in zip(idxs, alignments):
                A[idx] = alm
    if fuzzy:
        cogid2alm = {}
        # only align the first item
        for cogid, vals in etd.items():
            idxs, alms, strucs = [], [], []
            for val in vals:
                if val:
                    idxs += val
                    alms += [
                        wordlist[idx, segments].n[wordlist[idx, ref].index(cogid)]
                        for idx in val
                    ]
                    strucs += [
                        wordlist[idx, structure].n[wordlist[idx, ref].index(cogid)]
                        for idx in val
                    ]
            alignments = shrink_alignments(
                [
                    align_to_template(alm, struc, template, gap=gap)
                    for alm, struc in zip(alms, strucs)
                ]
            )
            for idx, alm in zip(idxs, alignments):
                cogid2alm[cogid, idx] = " ".join(alm)
        # second iteration, add the alignments per cogid
        for idx, cogids in wordlist.iter_rows(ref):
            A[idx] = bt.lists(
                " + ".join([cogid2alm.get((cogid, idx), "?") for cogid in cogids])
            )
    wordlist.add_entries(alignment, A, lambda x: x, override=override)
