# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
from collections import defaultdict
from itertools import combinations
import os
# TODO export this function to lingpy or something else
from sinopy.segments import get_structure
from lingpy.sequence.sound_classes import tokens2morphemes
from lingpy import *
from lingpy.align.sca import get_consensus
from lingpy import basictypes as bt
import networkx as nx
import html
import codecs


def complexity(wordlist, pos='A', segments='tokens', structure='structure'):
    """
    Return the worst-case scenario for correspondence patterns.
    
    Note
    ----
    Worst-case here means: for n languages with n₁, n₂, ..., sounds each, we
    would calculate the number of different calculations between them all.
    """
    # retrieve all sounds in the structure of the language
    languages = defaultdict(lambda : defaultdict(int))
    for idx, doc, tokens, strucs in wordlist.iter_rows(
            'doculect', segments, structure):
        if isinstance(strucs, str):
            strucs = strucs.split()
        for t, s in zip(tokens, strucs):
            if pos == s:
                languages[doc][t] += 1
    score = 1
    for doc in wordlist.cols:
        print('{0:20} | {1:20} '.format(doc, len(languages[doc])))
        score *= len(languages[doc])
    print('{0:20} | {1:20}'.format('TOTAL', score))
    return languages
    


def lingrex_path(*comps):
    """
    Our data-path in CLICS.
    """
    return os.path.join(os.path.dirname(__file__), os.pardir, *comps)

def data_path(*comps):

    return lingrex_path('data', *comps)

def get_c_structure(seg, cldf=True):
    """Bad-ass function to get the major structures for alignments in chinese
    data"""
    cls = ''.join(tokens2class(seg, 'cv', cldf=cldf))
    mapper = {
            'C': "i",
            "CVCC": "inNc",
            "CVV": "inN",
            "CVCV": "inIN",
            'CCV': 'imn',
            'CCVC': 'imnc',
            'CCVCT': 'imnct',
            'CCVT': 'imnt',
            'CCVVT' : '',
            'CVC': 'inc',
            'CT' : 'nt',
            'CVCT': 'inct',
            'VVCT': 'mnct',
            'CCCVT': '',
            'CCCVCT': '',
            'VC': 'nc',
            'VCT': 'nct',
            'CVVCT': 'inNct',
            'VVT': 'nNt',
            'CVT': 'int',
            'CVVT': '',
            'CCT': '',
            'V': 'n',
            'CV': 'in',
            'CCVT': 'imnt',
            'VT': 'nt',
            'CCVVCT': 'imnNct',
            'CVCCT': 'incCt',
            }
    if not cls in mapper:
        return '?' * len(seg)
    if not mapper[cls]:
        # our problem are VV instances, so we need to extract these
        dlg = ''.join(tokens2class(seg, 'sca', cldf=True))
        ncls, tcls = '', ''
        for c, d in zip(cls, dlg):
            if c == 'V':
                ncls += d
                tcls += c
            elif c == 'C':
                ncls += c
                tcls += d
            else:
                ncls += c
                tcls += c
        

        nmapper = {
            "UIT": 'nNt',
            "CCT": "int",
            "CUYT": 'inNt',
            "CYACT": 'imnct',
            "CYET": "imnt",
            "CYAT": 'imnt',
            "AYT": 'nNt',
            "CCUAT": 'imMnt',
            "CCYAT": 'imMnt',
            "CCAUT": 'imnNt',
            "CCAYT": 'imnNt',
            "CCIAT": 'imMnt',
            "CCAIT": 'imnNt',
            "CCEAT": 'imMnt',
            "CCYIT": 'imnNt',
            "CCEIT": 'imnNt',
            "CCUYT": 'imnNt',
            "CIYT": 'imnt',
            "CYIT": 'inNt',
            "CAIT": 'inNt',
            "CEIT": 'inNt',
            "CIUCT": 'imnct',
            "CAYT": 'inNt',
            "CIACT": 'imnct',
            "CAYCT": 'inNct',
            "CAICT": 'inNct',
            "CUIT": 'inNt',
            "YIT": 'nNt',
            "CYUT": "imnt",
            "CIET": "imnt",
            "CIUT": "imnt",
            "CIAT": "imnt",
            "CAUT": "inNt",
            "CIIT": "imnt",
            "EET": 'mnt' 
            }

        cmapper = {
            "CNT": "int",
            "TNT": "int",
            "KNT": "int",
            "MNT": "int",
            "SNT": "int",
            "MMT": 'int',
            "NNT": "int",
            "MSWVT": "iMmnt",
            "MSWVTT": "iMmnct",
            "HMT": "int",
            "KSWVT": "iMmnt",
            "HNT": "int",
            "GNT": "int",
            "PLWVT": "iMmnt",
            "PSWVT": "iMmnt",
            "MSWVT": "iMmnt",
            "TLWVT": "iMmnt",
            "GSWVT": "iMmnt",
            "KSWVHT": "iMmnct",
            "PSWVHT": "iMmnct",
            "MSWVHT": "iMmnct",
            "KSWVTT": "iMmnct",
            "MSWVNT": "iMmnct",
            "LLT": "int",
                }
        if ncls in nmapper:
            return nmapper[ncls]
        
        if tcls in cmapper:
            return cmapper[tcls]
        
        print(ncls,seg,tcls)
        #input()
        return '?' * len(seg)

    return mapper[cls]

def get_segments_and_structure(wordlist, etd, cogid, ref='cogids', segments='segments',
        structure='structure'):
    """shortcut function to retrieve a couple of things from a wordlist"""
    idxs = []
    for v in etd[cogid]:
        if v: idxs += v
    out = {}
    for idx in idxs:
        segs = wordlist[idx, segments]
        cogids = wordlist[idx, ref]
        cogidx = cogids.index(cogid)
        morpheme = tokens2morphemes(segs)[cogidx]
        if isinstance(wordlist[idx, structure], (list, tuple)):
            struct = ' '.join(wordlist[idx, structure]).split(' + ')[cogidx]
        else:
            struct = wordlist[idx, structure].split(' + ')[cogidx]
        out[idx] = [wordlist[idx, 'doculect'], wordlist[idx, 'concept'], 
            morpheme, struct.split(' ')]

    return out


def renumber_partials(wordlist, ref):
    """Renumber function for partial cognates"""
    newidx = 1
    partials = {}
    for idx, concept, cogids in iter_rows(wordlist, 'concept', ref):
        new_cogs = []
        for cogid in cogids:
            _cogid = str(cogid)+'-'+concept
            if _cogid not in partials:
                partials[_cogid] = newidx
                newidx += 1
            new_cogs += [partials[_cogid]]
        wordlist[idx, ref] = new_cogs
   

def align_by_structure(wordlist, template='imMnNct', segments='segments',
        ref='cogids', structure='structure', alignment='alignment',
        override=True):
    """Align patterns simply by following the template"""
    etd = wordlist.get_etymdict(ref=ref)
    alms = {}
    for key in etd:
        # get the essential data
        alm, idxs, cidxs = [], [], []
        for idx, (d, c, m, s) in get_segments_and_structure(wordlist, etd, key,
                ref, segments, structure).items():
            alm += [[]]
            mapper = dict(zip(s, m))
            for t in template:
                alm[-1] += [mapper.get(t, '-')]
            idxs += [idx]
        ignore = []

        for i in range(len(alm[0])):
            col = [alm[j][i] for j in range(len(idxs))]
            if not [x for x in col if x != '-']:
                ignore += [i]
        out = []
        for idx, row in zip(idxs, alm):
            out = [row[i] for i in range(len(alm[0])) if i not in ignore]
            alms[key, idx] = out

    alignments = {}
    for idx, cogids in iter_rows(wordlist, ref):
        alignments[idx] = ' + '.join([' '.join(alms[cogid, idx]) for cogid in
            cogids]).split(' ')
    wordlist.add_entries(alignment, alignments, lambda x: x, override=override)
    if hasattr(wordlist, 'add_alignments'):
        wordlist.add_alignments(override=True)


def add_structure(wordlist, model='cv', segments='tokens',
        structure='structure', ref='cogid', gap='-'):
    """Add structure to a wordlist to make sure correspondence patterns can be
    inferred"""
    if model not in ['cv', 'c', 'CcV', 'ps', 'nogap']:
        raise ValueError('[i] you need to select a valid model')
    D = {}
    if model == 'cv':
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = ' '.join(tokens2class(tks, 'cv')).lower()
    
    if model == 'c':
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = ' '.join(tokens2class(tks, 'cv')).lower().replace('v',
                    'c').replace('t', 'c')
    if model == 'nogap':
        assert hasattr(wordlist, 'msa')
        for cogid, msa in wordlist.msa[ref].items():
            cons = ['c' if c != gap else gap for c in get_consensus(msa['alignment'], 
                gaps=True)]
            for idx, alm in zip(msa['ID'], msa['alignment']):
                struc = []
                for a, b in zip(cons, alm):
                    if b != '-':
                        struc += [a]
                D[idx] = ' '.join(struc)
        for idx, tks in wordlist.iter_rows(segments):
            if idx not in D:
                D[idx] = ' '.join(['c' if c != '+' else c for c in tks])
    if model == 'CcV':
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = ' '.join(list(prosodic_string(tks,
                _output='CcV').replace('_', '+')))
    if model == 'ps':
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = ' '.join(list(prosodic_string(tks)))

    struc_ = bt.lists if wordlist._mode == 'fuzzy' else bt.strings
    wordlist.add_entries(structure, D, lambda x: struc_(x))


def add_c_structure(wordlist, segments='tokens', structure='structure', sep=' +'):
    """Add the structure data on sounds we need for our analysis"""
    structures = {}
    ssep = ' '+sep+' '
    for idx, segs in iter_rows(wordlist, segments):
        strucs = []
        for mrp in tokens2morphemes(segs):
            strucs += [' '.join(get_c_structure(mrp))]
        structures[idx] = ssep.join(strucs)
    wordlist.add_entries(structure, structures, lambda x: bt.strings(x))
    

def save_network(filename, graph):
    with codecs.open(filename, 'w', 'utf-8') as f:
        for line in nx.generate_gml(graph):
            f.write(html.unescape(line)+'\n')


def load_network(filename):
    with codecs.open(filename, 'r', 'utf-8') as f:
        lines = [l.encode('ascii', 'xmlcharrefreplace').decode('utf-8') for l in f]
        return nx.parse_gml('\n'.join(lines))
