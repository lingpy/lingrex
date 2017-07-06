# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
from collections import defaultdict
from itertools import combinations
from lingpy.sequence.sound_classes import tokens2morphemes

import networkx as nx
from networkx.algorithms.approximation.clique import clique_removal
from lingpy import *
from lingpy.util import pb
from lingpy.align.sca import get_consensus, SCA
from lingpy.compare.partial import Partial, _get_slices
from .util import save_network
from six import text_type
from lingpy.read.qlc import normalize_alignment
from lingrex.util import get_simple_structure

def diversity(wordlist, ref='cogid', concept='concept'):
    if not ref in wordlist.header or not concept in wordlist.header:
        raise ValueError("Wrong arguments for function diversity")
    etd = wordlist.get_etymdict(ref=ref, entry=concept)
    counter = defaultdict(list)
    for key, vals in etd.items():
        for val in filter(None, vals):
            for v in val:
                counter[v] += [key]
    return sum([(len(set(vals))-1)/(len(vals)-1) for vals in
        counter.values()])/len(counter)
    
def regularity(wordlist, graph, ref='cogid', clique='clique', debug=False,
        threshold=2, gap="Ø"):
    """Central measure for neogrammarian correspondence regularity.

    Note
    ----
    The regularity measure measures how many words in a dataset can be assigned
    to regular cognate sets. Regularity of cognate sets is defined on the
    cliques in a compatibility graph, by applying a numeric threshold that
    defines the average number of examples for a given pattern, excluding
    languages which do not take part in the given cognate sets. The measure is
    applied to all sounds in a given sound template, following the rules of
    construction for the compatibility graph.
    """
    etd = wordlist.get_etymdict(ref=ref)
    counter = defaultdict(list)
    wc = rc = 0 
    for k, vals in etd.items():
        if str(k) in graph:
            counter[graph.node[str(k)]['clique']] += [k]
        wc += len(list(filter(None, vals)))

    regular_words = {}
    for k, v in counter.items():
        nodes = [n for n in graph.nodes() if graph.node[n]['clique'] == k]
        if score_pattern(graph, nodes) >= threshold:
            for cogid in v:
                rc += len(list(filter(None, etd[cogid])))
                for wset in [list(filter(None, etd[cogid]))]:
                    for word in wset:
                        regular_words[word[0]] = k
    

    return regular_words, rc / wc

def apply_regularity(wordlist, regular_words, wwc, ref='cogid', refb='cogid'):
    tmp = {}
    idx = 1
    for k, cogid in iter_rows(wordlist, ref):
        if k not in regular_words and wwc[k] == 1:
            tmp[k] = 'n-{0}'.format(idx)
            idx += 1
        else:
            tmp[k] = cogid
    wordlist.add_entries('_tmp', tmp, lambda x: x, override=True)
    wordlist.renumber('_tmp', refb, override=True)

def slice_segments(segments):
    return list(map(lambda x: segments[x[0]:x[1]], _get_slices(segments)))

def align_by_structure(wl, ref='cogids', segments='segments',
        structure='structure', template='i m n c t', override=True):
    """
    Compute alignments based on manually assigned structures for each segment.

    Parameters
    ----------
    ref : str (default="cogids")
        Column storing the cognate identifiers.
    segments : str (default='segments')
        Column storing the segment identifiers.
    structure : str (default='structure')
        Column storing the segment structure (the "context" in a broader
        sense).
    template : str (default='i m n c t')
        The template which is going to be searched for the alignments, that is,
        the maximal syllable structure.

    Notes
    -----
    This function requires accurately annotated data. So be careful when using
    it and make sure your wordlist file contains all the required data. The
    function adds to more columns to the wordlist: one containing alignments
    (default name: "alignments"), one containing the structural consensus,
    (default name "consensus").
    """
    if wl._class.get(ref) == text_type:
        for k in wl:
            if not isinstance(wl[k, ref], (tuple, list)):
                wl[k][wl.header[ref]] = [int(x) for x in wl[k, ref].split(' ')]
    # check for same number of cogids and slices
    errors = []
    for k in wl:
        if len(wl[k, ref]) != len(_get_slices(wl[k, segments])):
            errors += [k]
            print(k, wl[k, segments], wl[k, ref])
    if errors:
        raise ValueError('There were {0} errors in your wl cognates'.format(
                    len(errors)))
    etd = wl.get_etymdict(ref)

    # assemble alignments
    alms = {k: ['0' for x in wl[k, ref]] for k in wl}
    strucs = {k: ['0' for x in wl[k, ref]] for k in wl}
    for key, vals in etd.items():
        keys = list(map(lambda x: x[0], filter(None, vals)))
        taxa = list(map(lambda x: wl[x, 'doculect'], keys))
        idxs = list(map(lambda x: wl[x, ref].index(key), keys))
        morphemes = list(map(lambda x: slice_segments(wl[x[0],
            segments])[x[1]],
            zip(keys, idxs)))
        structures = list(map(lambda x: wl[x[0], structure].split(' + ')[x[1]],
            zip(keys, idxs)))

        _alms, _struc = [], []
        for m, s in zip(morphemes, structures):
            alm, struc = [], []
            _s = s.split(' ')
            for t in template.split(' '):
                if t in _s:
                    alm += [m[_s.index(t)]]
                    struc += [t]
                else:
                    alm += ['-']
                    struc += ['-']
            _alms += [alm]
            _struc += [struc]
        _alms = normalize_alignment(_alms)
        _struc = normalize_alignment(_struc)
        
        # get the structure consensus
        cons = []
        for i in range(len(_struc[0])):
            c = set([line[i] for line in _struc if line[i] != '-'])
            if c:
                cons += [c.pop()]
        cons = ' '.join(cons)
        for k, alm, i in zip(keys, _alms, idxs):
            alms[k][i] = ' '.join(alm)
            strucs[k][i] = cons
    wl.add_entries('alignment', alms, lambda x: ' + '.join(x),
            override=override)
    wl.add_entries('prostring', strucs, lambda x: ' + '.join(x),
            override=override)

def _get_matrix(graph, clique, gap="Ø"):
    matrix = []
    for n in clique:
        matrix += [[1 if x != gap else 0 for x in graph.node[n]['column'].split(' ')]]
    return matrix

def weight_clique(graph, clique):
    matrix = _get_matrix(graph, clique)
    return sum([sum(line) for line in matrix]) / (len(matrix) * len(matrix[0]))

def score_pattern(graph, clique, gap="Ø"):
    """
    Weight the density of a given pattern.

    Note
    ----
    Essentially: check how well each column in a given pattern is filled by
    enough values. Empty columns are excluded in this analysis.
    """
    matrix = _get_matrix(graph, clique, gap=gap)
    m = len(matrix[0])
    zeros = len([
        0 for i in range(m) if not sum(
            [row[i] for row in matrix])]) 
    return sum([sum(line) for line in matrix]) / (m-zeros)

def _has_null(matrix):
    return any([not sum([line[i] for line in matrix]) for i in
        range(len(matrix[0]))])

def consensus_pattern(patterns, gap="Ø"):
    """Return consensus pattern of multiple patterns.

    Parameters
    ----------
    patterns : list
        List of patterns (each pattern should be a sequence, preferably a
        list).
    gap : str (default="Ø")
        A gap in the sense of "missing data", that is, a cognate set for which
        a value in a given language is absent.

    Returns
    -------
    consensus : list
        A one-dimensional list.

    Note
    ----
    This consensus method raises an error if the patterns contain incompatible
    columns (non-identical values apart from the gap character in the same
    column). 
    """
    out = []
    for i in range(len(patterns[0])):
        col = [line[i] for line in patterns]
        no_gaps = [x for x in col if x != gap]
        if len(set(no_gaps)) > 1:
            raise ValueError("Your patterns are incompatible")
        out += [no_gaps[0] if no_gaps else gap]
    return out

def compatible_columns(colA, colB, missing='Ø'):
    """Check for column compatibility.

    Parameters
    ----------
    colA, colB = list
        Lists (sequence type) containing a given pattern.
    missing : str (default="Ø")
        A gap in the sense of "missing data", that is, a cognate set for which
        a value in a given language is absent.
    
    Returns
    -------
    score : int
        The score is 0 if there is no conflict, -1 if there are conflicts, and
        n for the number of identical non-empty values in the column if the
        patterns are compatible.
    """

    matches = 0
    for a, b in zip(colA, colB):
        if not missing in [a, b]:
            if a != b:
                return -1
            else:
                matches += 1
    return matches

def nodes_with_condition(alms, ref='cogid', pos='A', prostring=False):
    """
    Calculate the patterns which match a certain condition and those which do
    not match the condition.
    """
    C = {}
    singletons = [k for k, c in iter_rows(alms, ref) if c not in alms.msa]
    for cogid, msa in alms.msa[ref].items():
        if not prostring:
            consensus = get_consensus(msa['alignment'], gaps=True)
            _prostring = prosodic_string(consensus)
            pidx = _prostring.find(pos)
        else:
            fidx = msa['ID'][0]
            consensus = get_consensus(msa['alignment'], gaps=True)
            _prostring = alms[fidx, prostring].split(' + ')
            _prostring = _prostring[alms[fidx, ref].index(cogid)].split(' ')
            if pos in _prostring:
                pidx = _prostring.index(pos)
            else:
                pidx = -1
        if pidx != -1:
            for idx in msa['ID']:
                C[idx] = 1
        else:
            for idx in msa['ID']:
                C[idx] = 0

    for idx in singletons:
        _prostring = prosodic_string(alms[idx, alms._segments])
        if pos in _prostring:
            C[idx] = 1
        else:
            C[idx] = 0
    return C

def compatibility_graph(
        alms, 
        ref='cogid',
        pos='A',
        minimal_taxa=3,
        taxa=None, 
        missing="Ø",
        color=None,
        debug=False,
        prostring=False,
        _filter=False
        ):
    """
    Get a graph of all compatible positions in a set of aligned strings.

    Parameters
    ----------
    alms : ~lingpy.basic.sca.Alignments
        Valid alignment object, containing partial or full alignments.
    ref : str (default="cogid")
        Reference column storing the alignments.
    pos : str (default="A")
        Prosodic position which will be investigated. Follows general handling
        of prostrings in lingpy.
    minimal_taxa : int (default=3)
        Determine the minimal amount of taxonomic units which should reflect a
        given pattern.
    taxa : list (default=None)
        Use this keyword to pass the desired number and order of taxonomic
        units which should be investigated.
    missing : str (default="Ø")
        A gap in the sense of "missing data", that is, a cognate set for which
        a value in a given language is absent.
    color : ~lingpy.data.models.Model (default='color')
        A Model instance of lingpy, the color-model (used for sound coloring)
        as a default.

    Note
    ----
    Compatibility is hereby defined as triadic compatibility: Two nodes are
    compatible, if the patterns contain either
    * one identical position and no non-identical positions (excluding gaps),
      or
    * they are compatible to a third node which contains one identical position
     with both nodes

    """
    if not _filter: _filter=lambda x: x
    graph = nx.Graph()
    stats = [0, 0]

    taxa = taxa or alms.taxa
    color = color or Model('color')
    good_patterns = []
    
    # add the patterns as nodes to the graph
    with pb(desc='COMPATIBILITY GRAPH CREATION', total=len(alms.msa[ref])) as progress:
        for cogid, msa in alms.msa[ref].items():
            progress.update(1)
            _taxa = msa['taxa']
            if len(set(_taxa)) >= minimal_taxa and _filter(msa):
                stats[0] += 1
                if not prostring:
                    consensus = get_consensus(msa['alignment'], gaps=True)
                    _prostring = prosodic_string(consensus)
                    pidx = _prostring.find(pos)
                else:
                    # TODO make this more properly responding here TODO
                    fidx = msa['ID'][0]
                    # wrong line here! we need to take the structure we already
                    # have, make also sure to normalize representation of
                    # structure
                    strucs = [
                            alms[idx, prostring].replace(' ', '').split('+')[
                                alms[idx, ref].index(cogid)]
                                for idx in msa['ID']
                                ]
                    #strucs = [get_simple_structure([x for x in s if x != '-']) for s in msa['seqs']]
                    print(strucs, msa)
                    msa_strucs = [class2tokens(struc, alm) for struc, alm in
                        zip(strucs, msa['alignment'])]
                    consensus = get_consensus(msa['alignment'], gaps=True)
                    _prostring = []
                    for _i in range(len(msa_strucs[0])):
                        row = [x[_i] for x in msa_strucs if x[_i] != '-']
                        _prostring += [row[0]]

                    #_prostring = get_simple_structure(consensus)
                    #if not _prostring:
                    #    _prostring = alms[fidx, prostring].split(' + ')
                    #    _prostring = _prostring[alms[fidx, ref].index(cogid)].split(' ')
                    if len(_prostring) != len(msa['alignment'][0]):
                        print(_prostring)
                        print(consensus)
                        for m in msa['alignment']:
                            print(m)
                        input('not good')

                    if pos in _prostring:
                        pidx = _prostring.index(pos)
                    else:
                        pidx = -1
                if pidx != -1:
                    stats[1] += 1
                    good_patterns += [cogid]
                    reflexes = []
                    for t in taxa:
                        if t not in _taxa:
                            reflexes += [missing]
                        else:
                            reflexes += [msa['alignment'][_taxa.index(t)][pidx]]
                    graph.add_node(str(cogid), column = ' '.join(reflexes),
                            consensus=consensus[pidx], clique=0, cliquesize=0,
                            color = tokens2class(consensus, color)[0])

    # determine compatibility and assign edges
    zero_edges = []
    for nA, nB in combinations(graph.nodes(), r=2):
        colA, colB = (
                graph.node[nA]['column'].split(' '),
                graph.node[nB]['column'].split(' ')
                )
        cc = compatible_columns(colA, colB, missing=missing)
        if cc >= 0:
            graph.add_edge(nA, nB, weight=cc)
            if cc == 0:
                zero_edges += [(nA, nB)]
    if debug: 
        print('[i] CA #1, {0} edges, {1} zero edges.'.format(
            graph.number_of_edges(),
            len(zero_edges)))

    # remove edges where we don't find a triad-compatibility
    remove_edges = []
    for i, (nA, nB) in enumerate(zero_edges):
        neighborsA = {n for n, d in graph.edge[nA].items() if d['weight']}
        neighborsB = {n for n, d in graph.edge[nB].items() if d['weight']}
        common_nodes = neighborsA.intersection(neighborsB)
        triadic_compatibility = False
        for n in common_nodes:
            if graph.edge[n][nA]['weight'] and graph.edge[n][nB]['weight']:
                triadic_compatibility = True
                break
        if not triadic_compatibility:
            remove_edges += [(nA, nB)]
        if not i % 500 and debug:
            print('[i] analyzing {0}'.format(i))
    if debug:
        print('[i] Removing {0} edges.'.format(len(remove_edges)))
    graph.remove_edges_from(remove_edges)

    if debug: 
        print('Patterns in total: {0}\nPatterns with condition: {1}'.format(stats[0], 
            stats[1]))
    return graph

def retrieve_identical_morphemes(alms, segments='segments', ref='cogids'):
    graph = nx.Graph()
    for col in alms.cols:
        keys = alms.get_list(col=col, flat=True)
        concepts = list(map(lambda x: alms[x, 'concept'], keys))
        M = defaultdict(list)
        for k, c in zip(keys, concepts):
            for i, morph in enumerate(slice_segments(alms[k, segments])):
                M[' '.join(morph)] += [(k, i, c)]
        for key, vals in M.items():
            for k, i, c in vals:
                cog = alms[k, ref][i]
                try:
                    graph.node[cog]['language'] += [col]
                    graph.node[cog]['concepts'] += [c]
                except KeyError:
                    graph.add_node(cog, language=[col], concepts=[c])
            for (k1, i1, c1), (k2, i2, c2) in combinations(vals, r=2):
                cog1, cog2 = alms[k1, ref][i1], alms[k2, ref][i2]
                try:
                    graph.edge[cog1][cog2]['language'] += [col]
                    graph.edge[cog1][cog2]['idxs'] += [(str(k1), str(k2))]
                    graph.edge[cog1][cog2]['concepts'] += [(c1, c2)]
                except KeyError:
                    graph.add_edge(cog1, cog2, idxs = [(str(k1), str(k2))],
                        language=[col], concepts=[(c1, c2)])
    for nA, nB, d in graph.edges(data=True):
        if len(d['language']) > 3:
            for l in alms.taxa:
                tmp = ''
                if l in graph.node[nA]['language']:
                    idxA = graph.node[nA]['language'].index(l)
                else:
                    idxA = -1
                if l in graph.node[nB]['language']:
                    idxB = graph.node[nB]['language'].index(l)
                else:
                    idxB = -1       
                lang1 = nA if l in graph.node[nA]['language'] else 0
                lang2 = nB if l in graph.node[nB]['language'] else 0
                cA = graph.node[nA]['concepts'][idxA] if idxA > -1 else "Ø"
                cB = graph.node[nB]['concepts'][idxB] if idxB > -1 else "Ø"
                link = '<->' if l in d['language'] else ''
                print('{0:15}\t{1:5}\t{2:15}\t{3:5}\t{4:5}\t{5:15}'.format(
                    l, lang1, cA, link, lang2, cB))
            print('')
    return graph

def collapse_identical_nodes(graph, column='column'):
    """Collapse those nodes in a graph which are completely identical."""
    D = defaultdict(list)
    for node, data in graph.nodes(data=True):
        D[data[column]] += [node]
    remove = []
    for column, nodes in D.items():
        graph.node[nodes[0]]['nodebunch'] = ','.join([str(x) for x in nodes])
        remove += nodes[1:]
    graph.remove_nodes_from(remove)

def greedy_clique_partition(graph, debug=True, weight_clique=None,
        _approximate=False):
    """Partition the graph in the minimal number of cliques.

    Parameters
    ----------
    graph : networkx.Graph
        Graph object, of the networkx-library. 
    weight_clique : function (default=len)
        A function that handles how the cliques are sorted (essential for the
        algorithm). Defaults to the length of the cliques, but you can pass any
        function that returns an integer and takes a set of nodes as input.

    Returns
    -------
    graph : networkx.Graph
        The graph object, with two additional attributes: one for each clique,
        one for each clique size.

    Note
    ----
    This algorithm is a greedy version of the well-known problem of clique
    partitioning, that is, the problem of finding the minimal number of cliques
    which partition a given graph.
    """
    weight_clique = weight_clique or len
    if _approximate:
        cliques = sorted(
                clique_removal(graph)[1], 
                key=lambda x: weight_clique(x), reverse=True
                )
    else:
        cliques = sorted(nx.find_cliques(graph), key=lambda x: weight_clique(x),
                reverse=True)

    visited = []
    idx = 1
    with pb(desc='GREEDY CLIQUE PARTITION', total=len(cliques)) as progress:
        while cliques:
            progress.update(1)
            next_clique = cliques.pop(0)
            not_visited = [n for n in next_clique if n not in visited]
            nvl = len(not_visited)
            for node in not_visited:
                graph.node[node]['clique'] = idx
                graph.node[node]['cliquesize'] = nvl
            visited += not_visited
            idx += 1
            cliques = [c for c in cliques if len([y for y in c if y not in
                visited]) > 1]
            cliques = sorted(cliques, key=lambda x: weight_clique(x), reverse=True)
            if debug: 
                print('{0} cliques left'.format(len(cliques)))

    # assign new indices to graph for all nodes with clique as "0"
    for node in graph.nodes():
        if graph.node[node]['clique'] == 0: 
            graph.node[node]['clique'] = idx
            graph.node[node]['cliquesize'] = 1
            idx += 1
    return graph

def score_graph(graph, clique='clique'):
    cliques = set([d[clique] for n, d in graph.nodes(data=True)])
    scores = []
    for c in cliques:
        clq = [n for n in graph.nodes() if graph.node[n][clique] == c]
        matrix = _get_matrix(graph, clq)
        if _has_null(matrix):
            scores += len(clq) * [0]
        else:
            score = score_pattern(graph, clq)
            if score < 1:
                scores += len(clq) * [0]
            else:
                scores += [score]
            #scores += [score_pattern(graph, clq)]
    return scores

def get_consensus_patterns(graph, debug=False, missing="Ø"):
    """
    To test the 
    """
    cliques = [x[1]['clique'] for x in graph.nodes(data=True)]
    consensus = {}
    for clique in cliques:
        patterns = [x[1]['column'].split(' ') for x in graph.nodes(data=True) if x[1]['clique'] == clique]
        nodes = ', '.join([x[0] for x in graph.nodes(data=True) if x[1]['clique'] ==
            clique])
        consensus[clique] = (consensus_pattern(patterns), len(patterns), nodes)

    for (clique1, (c1, n1, d1)), (clique2, (c2, n2, d2)) in combinations(consensus.items(), r=2):
        cc = compatible_columns(c1, c2, missing=missing)
        if cc > 0:
            print(clique1, '\t', n1, '\t', '\t'.join(c1))
            print(clique2, '\t', n2, '\t', '\t'.join(c2))
            print('---{0}---'.format(cc))

    return consensus


def find_source(consensus, graph):
    """find the source in a sound change graph for a given subset"""
    alt = reconstruct_from_consensus(consensus, rule='majority')
    sounds = [s for s in consensus if s != 'Ø']
    if len(set(sounds)) == 1:
        return '*'+sounds[0]
    if [s for s in sounds if s not in graph]:
        return alt
    
    subgraph = nx.subgraph(graph, sounds)
    # check for connected graph
    if [n for n in subgraph.nodes() if (
        nx.is_isolate(subgraph, n) and not nx.ancestors(subgraph, n)
        )]:
        return alt
    queue = [subgraph.nodes()[0]]
    sources = []
    visited = []
    while queue:
        node = queue.pop(0)
        visited += [node]
        ancestors = list(nx.ancestors(subgraph, node))
        if not ancestors:
            sources += [node]
        else:
            queue += [ancestor for ancestor in ancestors if ancestor not in
                    visited]
    if not sources:
        return alt

    return '*'+'/'.join(sorted(set(sources)))


def reconstruct_from_consensus(consensus, rule='majority', missing='Ø',
        graph=None):
    """apply different ways to reconstruct from the consensus patterns"""
    if rule == 'majority':
        return sorted([c for c in consensus if c != missing], key=lambda x: consensus.count(x), reverse=True)[0]
    elif rule == 'network':
        if not graph:
            raise ValueError("you must provide a networkx.graph as network")
        return find_source(consensus, graph)


def parse_sound_change_graph(path, source=1, target=2):
    """parse a source-target list of sound changes and make it a graph"""
    csv = csv2list(path, strip_lines=False)
    graph = nx.DiGraph()
    for line in csv[1:]:
        if len(line) > 2 and line[source].strip() and line[source].strip() != line[target].strip():
            graph.add_edge(line[source], line[target])
        else:
            graph.add_node(line[source].strip())
    return graph


def add_proto_slots(alignments, proto_name, structure='structure', threshold=3,
        segments='segments', unknown='?', ref='cogids', alignment='alignment',
        template='imMnNct', consensus=False, rule='majority', graph=None):
    """Add proto-language slots where there are enough reflexes"""

    def create_pattern(alms, idxs, wordlist):
        """helper function extracts the consensus patterns from an alignment"""
        these_taxa = [wordlist[idx, 'doculect'] for idx in idxs]
        out = []
        for i in range(len(alms[0])):
            alm = [row[i] for row in alms]
            out += [[]]
            for t in wordlist.taxa:
                if t in these_taxa:
                    out[-1] += [alm[these_taxa.index(t)]]
                else:
                    out[-1] += ["Ø"]
        return out
    
    concepts = defaultdict(lambda : defaultdict(list))
    for idx, concept, cogids in iter_rows(alignments, 'concept', ref):
        for cogid in cogids:
            concepts[concept][cogid] += [idx]
    plang = {}
    widx = max(alignments)+1
    for concept, cogids in concepts.items():
        pids, pforms, pstrucs, pinfos = [], [], [], []
        for cogid, reflexes in cogids.items():
            # make sure to have attested forms only if there are enough
            # languages, 
            languages = [alignments[r, 'doculect'] for r in reflexes]
            if len(set(languages)) >= threshold:
                pids += [cogid]
                rstrucs, alms, structs = [], [], []
                for idx in reflexes:
                    print(alignments.header, alignments[idx])
                    midx = alignments[idx, ref].index(cogid)
                    rstrucs += [alignments[idx, 'structure'].split(' + ')[midx]]
                    alms += [tokens2morphemes(alignments[idx,
                        alignment])[midx]]
                    
                structs = [class2tokens(s.split(), a) for s, a in zip(rstrucs, alms)]
                pstruc = []
                for line in zip(*structs):
                    pstruc += [[x for x in line if x != '-'][0]]
                pstrucs += [' '.join(pstruc)]
                
                if not consensus:
                    pforms += [' '.join([unknown for x in pstruc])]
                else:
                    pform, pinfo = [], []
                    patterns = create_pattern(alms, reflexes, alignments)
                    for p, ctx in zip(patterns, pstruc):
                        sounds, pattern_info = reconstruct_by_pattern(p, ctx, consensus[ctx],
                                threshold, rule=rule, graph=graph)
                        pform += ['|'.join(sounds)]
                        pinfo += ['<tr><td>'+sounds[0]+'</td><td>'+'</td><td>'.join(
                            pattern_info[0])+'</td></tr>']
                    if '?' in sounds:
                        print(concept)
                        print(cogids)
                        print(cogid)
                        print(reflexes)
                        print(patterns)
                        print(pstruc)
                        print(pform)
                        print(languages)
                        input()
                    pforms += [' '.join(pform)]
                    pinfos += ['<table border=1>'+ ''.join(pinfo)+'</table>']

        if pids:
            print('---', pinfos, pforms, pids, ref)
            plang[widx] = {
                'doculect': proto_name, 
                'concept': concept, 
                segments: ' '.join([x for x in ' + '.join(pforms).split() if x !=
                    '-']),
                structure: ' + '.join(pstrucs),
                alignment: ' + '.join(pforms),
                'patterns': '<br>'.join(pinfos),
                ref: pids
                }
            widx += 1
    D = {}
    D[0] = sorted(
        alignments.header, key=lambda x: alignments.header[x])
    for idx in alignments:
        D[idx] = [alignments[idx, h] for h in D[0]] + ['']
    D[0] += ['patterns']
    for idx in plang:
        D[idx] = [plang[idx].get(h, '') for h in D[0]]
    return Wordlist(D)



def reconstruct_by_pattern(pattern, context, patterns, threshold=3,
        output='string', rule='majority', graph=None, missing="Ø"):
    """Given one pattern for a couple of languages, find the best matching pattern
    
    Note
    ----
    Patterns are a consensus dictionary.
    """
    # we make a simple matching procedure for compatibility to find all
    # patterns with maximal matches
    good_patterns = defaultdict(list)
    less_good_patterns = defaultdict(list)
    for idx, (p, _l, _n) in patterns.items():
        if _l > 1:
            score = compatible_columns(pattern, p, missing=missing)
            if score >= threshold:
                good_patterns[score] += [idx]
        else:
            score = compatible_columns(pattern, p, missing=missing)
            if score >= threshold:
                less_good_patterns[score] += [idx]

    # take the best patterns as a max-argument
    if good_patterns:
        best_patterns = good_patterns[max(good_patterns)]
        # infer proto-forms for the best ones
        protos = defaultdict(list)
        for idx in best_patterns:
            proto_form = reconstruct_from_consensus(
                    patterns[idx][0], rule=rule, graph=graph)
            protos[proto_form] += [idx]
        if output == 'string':
            sorted_protos = sorted(protos, key=lambda x: len(protos[x]),
                    reverse=True)
            for key, vals in protos.items():
                for i, val in enumerate(vals):
                    vals[i] = patterns[idx][0]
            return sorted_protos, [protos[x][0] for x in sorted_protos]

    elif less_good_patterns:
        best_patterns = less_good_patterns[max(less_good_patterns)]
        # infer proto-forms for the best ones
        protos = defaultdict(list)
        for idx in best_patterns:
            proto_form = reconstruct_from_consensus(
                    patterns[idx][0], rule=rule, graph=graph)
            protos[proto_form] += [idx]
        if output == 'string':
            sorted_protos = sorted(protos, key=lambda x: len(protos[x]),
                    reverse=True)
            for key, vals in protos.items():
                for i, val in enumerate(vals):
                    vals[i] = patterns[idx][0]
            return ['!'+x.replace('*', '') for x in sorted_protos], [protos[x][0] for x in sorted_protos]
    else:
        print(pattern, context, good_patterns)
        return '?', '?'


