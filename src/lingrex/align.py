from lingpy import *
from lingpy import basictypes as bt

def align_to_template(sequence, structures, template, gap="-"):
    
    assert len(sequence) == len(structures) and len(template) >= len(sequence)
    assert len([x for x in structures if x not in template]) == 0

    out = []
    idxA, idxB = 0, 0
    while idxB < len(template):
        if idxA < len(sequence):
            segment, structure = sequence[idxA], structures[idxA]
        else:
            segment, structure = gap, ''
        current_structure = template[idxB]
        if current_structure == structure:
            out += [segment]
            idxA += 1
        else:
            out += [gap]
        idxB += 1
        
    return out

def shrink_alignments(alignments, gap="-"):
    excludes = []
    for i in range(len(alignments[0])):
        col = set([line[i] for line in alignments])
        if '-' in col and len(col) == 1:
            excludes += [i]
    out = []
    for alignment in alignments:
        out += [[site for i, site in enumerate(alignment) if i not in
            excludes]]
    return out

def template_alignment(
        wordlist,
        ref='cogid',
        template='CCCCVVccccT_CCCCVVccccT_CCCCVVccccT_CCCCVVccccT_CCCCvvT',
        structure='structure',
        fuzzy=False,
        segments='tokens',
        gap="-",
        alignment='alignment'
        ):
    
    for idx, tokens, structures in wordlist.iter_rows(segments, structure):
        wordlist[idx, segments], wordlist[idx, structure] = bt.lists(
                tokens), bt.lists(structures)

    etd = wordlist.get_etymdict(ref)
    A = {}
    if not fuzzy:
        for cogid, vals in etd.items():
            idxs = []
            for val in vals:
                if val:
                    idxs += val
            alignments = shrink_alignments([align_to_template(
                        wordlist[idx, segments],
                        wordlist[idx, structure],
                        template,
                        gap=gap) for idx in idxs])
            for idx, alm in zip(idxs, alignments):
                A[idx] = alm
    wordlist.add_entries(alignment, A, lambda x: x)
                

