from collections import OrderedDict
from lingpy.read.csv import csv2list
from lingpy.sequence.sound_classes import tokens2morphemes
from lingrex.util import *
from lingrex.copar import CoPaR, compatible_columns

def read_patterns(filename, taxa, structure='structure', cognates='cognates',
        proto=False, note=False):

    data = csv2list(filename, strip_lines=False)
    header = [h.lower() for h in data[0]]

    # get the taxa, to make sure to have the right order
    taxa_ = [t for t in taxa if t in data[0]]
    idxs_ = [i for i in range(len(header)) if data[0][i] in taxa_]

    # get the index of the structure
    sidx = header.index(structure.lower())
    cidx = header.index(cognates.lower())
    pidx = header.index(proto.lower()) if proto else -1
    nidx = header.index(note.lower()) if note else -1

    # create an index to change the taxonomic order
    tax2idx = {taxon: taxa.index(taxon) for taxon in taxa}

    # write patterns to ordered dict
    patterns = OrderedDict()
    patterns['patterns'] = []
    for line in data[1:]:
        reflexes = [line[i] for i in idxs_]

        if cognates:
            cogids = [x.split(':') for x in line[cidx].split(', ')]
        else:
            cogids = []

        pform = line[pidx] if proto else '?'
        struc = line[sidx] if structure else ''
        note_txt = line[nidx] if note else ''

        pattern = ['' for t in taxa]
        for i, (taxon, value) in enumerate(zip(taxa_, reflexes)):
            pattern[tax2idx[taxon]] = value

        patterns[tuple(pattern)] = {
                "cognates": [(int(x), int(y)) for x, y in cogids],
                "proto": pform,
                "structure": struc,
                "note": note_txt,
                "size": len(cogids)
                }
        patterns['patterns'] += [tuple(pattern)]

    patterns[tuple(taxa)] = {}
    patterns['taxa'] = tuple(taxa)

    return patterns

def to_super(number):
    string = str(number)
    st = list(zip(
        '¹²³⁴⁵⁶⁷⁸⁹⁰',
        '1234567890'))
    for s, t in st:
        string = string.replace(t, s)
    return string


def consensus_from_structures(strucs):
    """
    Helper function for creating a quick consensus.
    """
    out = []
    for i, s in enumerate('imnct'):
        if [x for x in strucs if s in x]:
            out += [s]
    return out

def reconstruct(wordlist, proto, patterns, ref='cogid', segments='tokens',
        alignment='alignment', minrefs=3, missing="Ø", threshold=3,
        structure='structure', uncertainty=2, frequency=2):

    # check for some things
    # TODO

    protos = {cogid: {} for cogid in wordlist.msa[ref]}
    freqs = defaultdict(list)
    for cogid, msa in wordlist.msa[ref].items():
        tax2idx = {taxon: patterns['taxa'].index(taxon) for taxon in patterns['taxa']}
        proto_segments = []

        # get the structures
        strucs = []
        for idx in msa['ID']:
            if wordlist._mode == 'fuzzy':
                cogidx = wordlist[idx, ref].index(cogid)
                strucs += [
                        tokens2morphemes(
                            wordlist[idx, structure].split())[cogidx]
                        ]
        scons = consensus_from_structures(strucs)

        for i in range(len(msa['alignment'][0])):
            protos[cogid][i] = []
            pattern = [missing for t in wordlist.cols]
            visited = set()
            for taxon, alm in zip(msa['taxa'], msa['alignment']):
                if taxon not in visited:
                    pattern[tax2idx[taxon]] = alm[i].split('/')[1] if '/' in \
                            alm[i] else alm[i]
                    visited.add(taxon)

            pattern = tuple(pattern)
            
            proto_form = defaultdict(int)
            # search for pattern
            for pt in sorted(patterns['patterns'], key=lambda x:
                    patterns[x]['size'], reverse=True):
                try:
                    scons[i]
                except:
                    print(scons, msa['ID'], cogid)
                if patterns[pt]['structure'] == scons[i] and \
                        patterns[pt]['size'] >= frequency and \
                        len([p for p in pt if p not in missing]) >= minrefs:
                    m1, m2 = compatible_columns(pattern, pt)
                    if m1 >= 1 and m2 == 0:
                        proto_form['{0}'.format(patterns[pt]['proto'] or '!')] += \
                                patterns[pt]['size']
                        freqs[pt] += [(cogid, i)]
                        protos[cogid][i] += [pt]

    
    out = {}
    for cogid, indices in protos.items():
        pform = []
        for i in indices:
            variants = defaultdict(int)
            count = 1
            for pt in sorted(protos[cogid][i], key=lambda x: len(freqs[x]),
                    reverse=True):
                variants[patterns[pt]['proto'].strip() or
                        '!'+':'.join(pt)] += len(freqs[pt])
                if not patterns[pt]['proto'].strip():
                    count += 1
            if not variants:
                pform += ['⁰?']
            else:
                start = 0
                variants_ = []
                for s, f in sorted(variants.items(), key=lambda x: x[1],
                        reverse=True):
                    if uncertainty > start:
                        variants_ += ["{0}{1}".format(to_super(f), s)]
                        start += 1
                pform += ['|'.join(variants_)]
        out[cogid] = pform
    return out


                

    
    
