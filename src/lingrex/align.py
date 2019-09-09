from lingpy import *
from lingpy import basictypes as bt
from tqdm import tqdm
from sinopy.segmentize import segmentize

def align_to_template(sequence, structures, template, gap="-"):
    try: 
        assert len(sequence) == len(structures) and len(template) >= len(sequence)
        assert len([x for x in structures if x not in template]) == 0
    except:
        print(sequence)
        print(structures)
        raise ValueError

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


def _shrink(segments, structure, converter={
            'i m': {'new': 'I'}, 
            'i':   {'new': 'I'}, 
            'n c': {'new': 'R'}, 
            'n':   {'new': 'R'}, 
            'c':   {'new': 'R'},
            ' ':   {'new': ''}}
            ):
    converted = segmentize(str(structure), converter)
    pos, out = 0, ['']
    for elms in converted:
        if elms.strip():
            before, before_ = '', False
            for elm in elms.split():
                if '/' in segments[pos]:
                    before_, next_char = segments[pos].split('/')
                    before +=  before_
                else:
                    before += segments[pos]
                    next_char = segments[pos]
                out[-1] += next_char
                pos += 1
            if before_:
                out[-1] = before +'/'+out[-1]
        else:
            out += ['']

    new_structure = segmentize(str(structure), converter, column='new')
    return ' '.join(out), ' '.join([n for n in new_structure if n])


def shrink_template(wordlist, structure='structure', segments='tokens',
        converter={
            'i m': {'new': 'I'}, 
            'i':   {'new': 'I'}, 
            'n c': {'new': 'R'}, 
            'n':   {'new': 'R'}, 
            'c':   {'new': 'R'},
            ' ':   {'new': ''},
            '+':   {'new': ' + '}},
        new_structure='structure2', 
        new_tokens='tokens2',
        new_alignment=False,
        override=False
            ):
    D = {}
    for idx, strucs, tokens in wordlist.iter_rows(structure, segments):
        D[idx] = _shrink(
            tokens, strucs, converter)
    wordlist.add_entries(new_structure, D, lambda x: bt.lists(x[1]),
            override=override)
    wordlist.add_entries(new_tokens, D, lambda x: bt.lists(x[0]),
            override=override)

    if new_alignment:
        wordlist.add_entries(new_alignment, D, lambda x: bt.lists(x[0]),
                override=override)



def template_alignment(
        wordlist,
        ref='cogid',
        template='CCCCVVccccT_CCCCVVccccT_CCCCVVccccT_CCCCVVccccT_CCCCvvT',
        structure='structure',
        fuzzy=False,
        segments='tokens',
        gap="-",
        alignment='alignment',
        override=True
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
    if fuzzy:
        cogid2alm = {}
        # only align the first item
        for cogid, vals in tqdm(etd.items(), desc='aligning cognates'):
            idxs, alms, strucs = [], [], []
            for val in vals:
                if val:
                    idxs += val
                    alms += [wordlist[idx, segments].n[wordlist[idx,
                        ref].index(cogid)] for idx in val]
                    strucs += [wordlist[idx, structure].n[wordlist[idx,
                        ref].index(cogid)] for idx in val]
            alignments = shrink_alignments([align_to_template(alm, struc,
                template, gap=gap) for alm, struc in zip(alms, strucs)])
            for idx, alm in zip(idxs, alignments):
                cogid2alm[cogid, idx] = ' '.join(alm)
        # second iteration, add the alignments per cogid
        for idx, cogids in wordlist.iter_rows(ref):
            A[idx] = bt.lists(' + '.join([cogid2alm.get((cogid, idx), '?') for cogid in
                cogids]))
    wordlist.add_entries(alignment, A, lambda x: x, override=override)

