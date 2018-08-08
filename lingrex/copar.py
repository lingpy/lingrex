# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
from six import text_type
from collections import defaultdict, OrderedDict
from itertools import combinations
from functools import cmp_to_key
import unicodedata

from lingpy.sequence.sound_classes import (
        prosodic_string, tokens2morphemes,
        class2tokens
        )
from lingpy.settings import rc
from lingpy.align.sca import get_consensus, SCA, Alignments
from lingpy.util import pb
from lingpy.compare.partial import Partial, _get_slices
from lingpy.read.qlc import normalize_alignment
from lingpy.algorithm.cython.misc import squareform
from lingpy.algorithm.extra import infomap_clustering
from lingpy import log
from lingpy import basictypes as bt

from lingrex.util import add_structure

def consensus_pattern(patterns, missing="Ø"):
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
        no_gaps = [x for x in col if x != missing]
        if len(set(no_gaps)) > 1:
            raise ValueError("Your patterns are incompatible")
        out += [no_gaps[0] if no_gaps else missing]
    return out

def incompatible_columns(patterns, missing="Ø"):
    """
    Compute whether a pattern has incompatible columns.
    """
    columns = []
    for i in range(len(patterns[0])):
        col = [patterns[j][i] for j in range(len(patterns)) if patterns[j][i]
                != missing]
        if len(set(col)) > 1:
            columns += ['*']
        else:
            columns += ['']
    return columns


def missing_in_pattern(pattern, missing="Ø"):
    """Utility function for quick access of missing data."""
    return pattern.count(missing)


def score_patterns(patterns, missing="Ø", mode='coverage', smooth=1):
    """Function gives a score for the overall number of reflexes.

    Notes
    -----
    This score tells simply to which degree a pattern is filled. It divides the
    number of cells not containing missing data by the number of cells in the
    matrix.
    """
    # we rank the columns by sorting them first
    if mode == 'ranked':
        if len(patterns) == 1:
            return 0
        cols = []
        for i in range(len(patterns[0])):
            cols += [sum([0 if row[i] == missing else 1 for row in patterns])]
        # sort the columns
        ranks, cols = list(range(1, len(cols)+1))[::-1], sorted(cols, reverse=True)
        scores = []
        for rank, col in zip(ranks, cols):
            scores += [rank * col]
        return sum(scores) / sum(ranks) / len(patterns)

    cols = []
    for i in range(len(patterns[0])):
        col = [row[i] for row in patterns]
        cols += [len(patterns) - col.count(missing)]
    if len(patterns) <= smooth:
        return 0
    return (sum(cols) / len(patterns[0])) / len(patterns) # * len(patterns[0]))

def compatible_columns(colA, colB, missing='Ø', gap='-'):
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
    matches, mismatches : tuple
        The score for matches gives zero if there is no conflict but also no
        match. For mismatches it is accordingly. So if you seek for
        compatibility, a mismatch greater 0 means the patterns are not
        compatible.
    """
    matches, mismatches = 0, 0
    for a, b in zip(colA, colB):
        if not missing in [a, b]:
            if a != b:
                mismatches += 1
            else:
                if a == gap:
                    pass
                else:
                    matches += 1
    return matches, mismatches


def density(wordlist, ref='cogid'):
    """Compute the density of a wordlist.

    Notes
    -----
    We define the density of a wordlist by measuring how many words can be
    explained by the same cognate set.
    """
    scores = []
    for concept in wordlist.rows:
        idxs = wordlist.get_list(row=concept, flat=True)
        cogids = [wordlist[idx, ref] for idx in idxs]
        sums = []
        for idx, cogid in zip(idxs, cogids):
            sums += [1 / cogids.count(cogid)]
        score = sum(sums) / len(sums)
        scores += [score]
    return 1-sum(scores) / len(scores)


class CoPaR(Alignments):
    """Correspondence Pattern Recognition class
    
    Parameters
    ----------
    wordlist : ~lingpy.basic.wordlist.Wordlist
        A wordlist object which should have a column for segments and a column
        for cognate sets. Since the class inherits from LingPy's
        Alignments-class, the same kind of data should be submitted.
    ref : str (default="cogid")
        The column which stores the cognate sets.
    segments : str (default="tokens")
        The column which stores the segmented transcriptions.
    alignment : str (default="alignment")
        The column which stores the alignments (or will store the alignments if
        they have not yet been computed).
    """
    def __init__(self, wordlist, minrefs=3, **keywords):
        Alignments.__init__(self, wordlist, **keywords)
        self.ref = keywords['ref']
        self._minrefs=minrefs

    def add_structure(self, model='cv', structure='structure', gap='-'):
        """Add structure to a wordlist.
        
        Parameters
        ----------
        model : str (default='cv')
            Select between "cv" for consonant-vowel(-tone) models, "c" for a
            model in which all segments are assigned to the same class, and
            "nogaps" where all segments are assigned to the same class but
            gappy sites in the alignments are excluded. 
        structure : str (default="structure")
            The name of the column in which the structure information should be
            placed.
        gap : str (default="-")
            The gap symbol to be used for the analysis.
        """
        add_structure(self, model=model, segments=self._segments,
                structure=structure, ref=self._ref, 
                gap=gap)
    
    # helper functions, generic
    def positions_from_alignment(self, alignment, pos):
        """Return positions matching from an alignment."""
        consensus = get_consensus(alignment, gaps=True)
        prostring = prosodic_string(consensus)
        return [i for i, p in enumerate(prostring) if p == pos]
    
    def positions_from_prostrings(self, cogid, indices, alignment, structures, pos):
        """Return positions matching from an alignment and user-defined prosodic strings"""
        if self._mode == 'fuzzy':
            strucs = []
            for idx, struc, alm in zip(indices, structures, alignment):
                pos_ = self[idx, self._ref].index(cogid)
                strucs += [class2tokens(tokens2morphemes(struc)[pos_], alm)]
        else:
            strucs = [class2tokens(struc, alm) for struc, alm in
                    zip(structures, alignment)]
        consensus = get_consensus(alignment, gaps=True)
        prostring = []
        for i in range(len(strucs[0])):
            row = [x[i] for x in strucs if x[i] != '-']
            prostring += [row[0] if row else '+']
        if len(prostring) != len(alignment[0]):
            print(prostring)
            print(consensus)
            print('\n'.join(['\t'.join(alm) for alm in alignment]))
            input('problematic alignment')
        return [i for i, p in enumerate(prostring) if p == pos]

    def reflexes_from_pos(self, position, taxa, current_taxa, alignment,
            missing, irregular):
        reflexes = []
        for t in taxa:
            if t not in current_taxa:
                reflexes += [missing]
            else:
                reflex = alignment[current_taxa.index(t)][position]
                if '/' in reflex:
                    reflex = reflex.split('/')[1]
                elif reflex[0] in irregular:
                    reflex = missing
                reflexes += [reflex]
        return reflexes

    def get_sites(
            self,
            ref=None,
            pos='A',
            minrefs=False,
            missing='Ø',
            structure=False,
            irregular = "!?"
            ):
        """
        Retrieve the alignment sites of interest for initial analysis.

        Parameters
        ----------
        ref : str (default=None)
            Method expects the reference to the cognate identifier column given
            upon initialization.
        pos : str (default='A')
            The position in the prostring representation of the segments in the
            data. As a default, prosodic strings are automatically created and
            only initial consonants are being compared.
        minrefs : int (default=2)
            The minimal number of reflexes per alignment site. 
        irregular : str (default='?!')
            The marker that indicates that a character is known to be irregular
            and has been marked as such in the alignments. 
        structure : str (default=False)
            Indicates the column in which the prosodic structure of the
            sequence can be found. If not provided, it will be automatically
            computed, using LingPy's prostring method.
        """
        ref = ref or self.ref
        patterns, all_patterns = OrderedDict(), OrderedDict()
        taxa = self.cols
        if minrefs:
            self._minrefs = minrefs

        for idx, struc in self.iter_rows(structure):
            self[idx, structure] = self._str_type(struc)
        
        # iterate over all sites in the alignment
        for cogid, msa in pb(sorted(self.msa[ref].items()), 
                desc='CoPaR: get_patterns()',
                total=len(self.msa[ref])):
            # get essential data: taxa, alignment, etc.
            _taxa = [t for t in taxa if t in msa['taxa']]
            _idxs = {t: msa['taxa'].index(t) for t in _taxa}
            _alms = [msa['alignment'][_idxs[t]] for t in _taxa]
            _wlid = [msa['ID'][_idxs[t]] for t in _taxa]
            if len(_taxa) >= minrefs:
                if not structure:
                    positions = self.positions_from_alignment(_alms, pos)
                else:
                    if self._mode == 'fuzzy':
                        _strucs = []
                        for _widx in _wlid:
                            _these_strucs = self[_widx, structure]
                            _strucs += [_these_strucs]
                    else:
                        _strucs = [self[idx, structure] for idx in \
                                _wlid]
                    positions = self.positions_from_prostrings(
                            cogid, _wlid, _alms, _strucs, pos)
                for pidx in positions:
                    reflexes = self.reflexes_from_pos(
                            pidx, taxa, _taxa, _alms, missing, irregular)
                    patterns[cogid, pidx] = reflexes
            for pidx in range(len(_alms[0])):
                reflexes = self.reflexes_from_pos(
                        pidx, taxa, _taxa, _alms, missing, irregular)
                all_patterns[cogid, pidx] = reflexes
        self.patterns = patterns
        self.all_patterns = all_patterns
        self._minrefs = minrefs
        return patterns, all_patterns 

    def sort_patterns(self, patterns=None, missing='Ø', threshold=2, debug=False):
        """
        Pre-sort the patterns to speed up the performance of the cluster algorithm.
        
        Parameters
        ----------
        patterns : dict (default=None)
            The correspondence patterns as identified by the algorithm.
            Defaults to None, meaning that the internally computed patterns
            will be used.
        missing : str (default="Ø")
            The symbol to be used for missing values.
        threshold : int (default=2)
            The threshold of the minimal number of cognate sets per cluster
            which should be regarded as regular.
        
        Notes
        -----
        The sorting algorithm sorts similar patterns straightforwardly into
        initial clusters and then runs through the sorted items, merging all
        compatible patterns into one cluster.
        """
        def sorter(x, y):
            """Sorter arranges patterns to allow for a quicker assembly"""
            if tuple(x[1]) == tuple(y[1]):
                return 0
            mg, mb = compatible_columns(x[1], y[1])
            tA = x[1].count(missing)
            tB = y[1].count(missing)
            if mg > 0:
                if tA < tB:
                    return -1
                if tA > tB:
                    return 1
                return (tuple(x[1]) > tuple(y[1])) - (tuple(x[1]) < tuple(y[1]))
            else:
                return (tuple(x[1]) > tuple(y[1])) - (tuple(x[1]) < tuple(y[1]))

        patterns = patterns or self.patterns
        sorted_patterns = sorted(patterns.items(), key=cmp_to_key(sorter))
        previous = sorted_patterns[0][1]
        out, indices = defaultdict(list), []
        for i, (key, pattern) in enumerate(sorted_patterns):
            match, mism = compatible_columns(pattern, previous)
            if mism > 0 or match == 0:
                for idx in indices:
                    out[tuple(previous)] += [sorted_patterns[idx][0]]
                indices = [i]
                previous = pattern
            else:
                indices += [i]
                previous = consensus_pattern([previous, pattern])
        for idx in indices:
            out[tuple(previous)] += [sorted_patterns[idx][0]]
        
        self.clusters = OrderedDict(**out)
        return out
        
    def cluster_patterns(self, clusters=None, debug=False,
            missing="Ø", match_threshold=1, gap='-', iterations=2,
            score_mode='coverage'):
        """Cluster patterns following an the spirit of Welsh-Powel algorithm.
        
        Parameters
        ----------
        clusters : dict (default=None)
            The clusters are assumed to be computed by the algorithm, so if set
            to None, the pre-computed clusters will be computed.
        missing : str (default="Ø")
            The symbol to be used for missing values.
        gap : str (default="-")
            If not set to "", gaps won't be counted when comparing two
            alignment sites
            for compatibility, in order to avoid problems from highly gappy
            sites.
        iterations : int (default=2)
            To make sure that during the 

        Notes
        -----
        This algorithm follows the spirit of the Welsh-Powell algorithm for
        graph coloring. Since graph coloring is the inverse of clique
        partitioning, we can use the algorithm in the same spirit.

        """
        clusters = clusters or self.clusters
        sorted_clusters = sorted(
                clusters.items(), 
                key=lambda x: (
                    len(x[0]) - x[0].count(missing),
                    score_patterns(
                        [self.patterns[y] for y in x[1]],
                        mode=score_mode
                        ),
                    len(x[1]),
                    ),
                reverse=True)
        out = []
        with pb(desc='CoPaR: cluster_patterns()', total=iterations*len(clusters)) as progress:
            while sorted_clusters:
                plen = len(sorted_clusters)
                (this_cluster, these_vals), remaining_clusters = sorted_clusters[0], sorted_clusters[1:]
                for i in range(iterations):
                    queue = []
                    for next_cluster, next_vals in remaining_clusters:
                        match, mism = compatible_columns(this_cluster, next_cluster,
                                missing=missing, gap=gap)
                        if match >= match_threshold and mism == 0:
                            this_cluster = consensus_pattern([this_cluster,
                                    next_cluster])
                            these_vals += next_vals
                        else:
                            queue += [(next_cluster, next_vals)]
                        remaining_clusters = queue
                out += [(this_cluster, these_vals)]
                sorted_clusters = sorted(
                        queue, 
                        key=lambda x: (
                            len(x[0]) - x[0].count(missing),
                            score_patterns(
                                [self.patterns[y] for y in x[1]],
                                mode=score_mode
                                ),
                            len(x[1])
                            ),
                        reverse=True)
                progress.update(iterations * (plen - len(queue)))
        clusters = {tuple(a): b for a, b in out}
        self.clusters = clusters
        self.ordered_clusters = sorted(clusters, key=lambda x: len(x[1]))
        alls = [c for c in clusters]
        match = 0
        for i, a in enumerate(alls):
            for j, b in enumerate(alls):
                if i < j:
                    ma, mi = compatible_columns(a, b)
                    if ma and not mi:
                        match += 1
        if match:
            log.warn('found {0} clusters which could be further merged'.format(match))
        if debug: print(single, mass, '{0:.2f}'.format(mass / (single + mass)))

    def sites_to_pattern(self, threshold=1, missing="Ø", debug=False):
        """Algorithm assigns alignment sites to patterns.

        Notes
        -----
        We rank according to general compatibility.
        """
        asites = defaultdict(list)
        best_patterns = sorted([c for c, s in self.clusters.items() if len(s) >=
                threshold])
        for consensus in pb(self.clusters, desc='CoPaR: sites_to_pattern()',
                total=len(self.clusters)):
            sites = self.clusters[consensus]
            for cog, pos in sites:
                pattern = self.patterns[cog, pos]
                for consensusB in self.clusters:
                    ma, mi = compatible_columns(pattern, consensusB)
                    if not mi and ma >= threshold:
                        asites[cog, pos] += [(ma, consensusB)]
        self.sites = asites

    def fuzziness(self):

        fuzziness = sum(
                [len(b) for a, b in self.sites.items()]) / len(self.sites)
        return fuzziness

    def refine_sites(self, missing='Ø'):
        if not hasattr(self, 'sites'):
            raise ValueError('run sites_to_pattern first')
        if not hasattr(self, 'clusters'):
            raise ValueError('run cluster detection first')
        new_clusters = defaultdict(list)
        new_sites = {}
        for site, patterns in pb(
                sorted(self.sites.items()),
                desc='COPAR: site to patterns',
                total=len(self.sites)
                ):
            spats = sorted(
                    patterns, 
                    key=lambda x: (
                        x[0], 
                        len(self.clusters[x[1]]), 
                        score_patterns([y for y in self.clusters[x[1]]],
                            mode='coverage')
                        ),
                    reverse=True)
            best_score, best_freq = spats[0][0], len(self.clusters[spats[0][1]])
            selected = []
            while spats:
                spat = spats.pop(0)
                new_score, new_freq = spat[0], len(self.clusters[spat[1]])
                if (best_score, best_freq) == (new_score, new_freq):
                    selected += [spat]
                else:
                    break
            new_sites[site] = selected
            new_clusters[selected[0][1]] += [site]

        self.clusters = OrderedDict()
        for p, sites in new_clusters.items():
            cons = consensus_pattern(
                    [self.patterns[s] for s in sites], missing=missing)
            self.clusters[tuple(cons)] = sites
        self.sites = new_sites



    def similar_patterns(self, mismatches=1, matches=3, debug=False, missing="Ø"):
        """
        Compute clusters by applying a community detection algorithm.

        Parameters
        ----------
        mismatches : int (default=1)
            The maximal amount of incompatible patterns allowed to make an edge
            in the graph.
        matches : int (default=3)
            The minimal amount of mutual matches (no missing characters per
            alignment) which will be accepted per pattern.
        missing : str (default="Ø")
            The symbol to be used for missing values.

        Notes
        -----
        This algorithm uses a community detection approach to cluster the
        patterns. Essentially, it allows for a certain degree of
        incompatibility, given that it is not an explicit solution of the
        clique partitioning problem.
        """
        if not hasattr(self, 'ordered_clusters'):
            raise ValueError("you need to compute clusters first")
        matrix = []
        for c1, c2 in combinations(self.ordered_clusters, r=2):
            match, misms = compatible_columns(c1, c2, missing=missing)
            if misms <= mismatches and match >= matches:
                matrix += [1]
            else:
                matrix += [3]
        matrix = squareform(matrix)

        clr = infomap_clustering(2, matrix, taxa=self.ordered_clusters)
        if debug:
            for key, value in sorted(clr.items(), key=lambda x: len(x[1]),
                    reverse=True):
                if len(value) > 2:
                    print('CLUSTER {0}'.format(key))
                    for pattern in value:
                        print('{0:5}\t'.format(
                            len(self.clusters[pattern]))+'\t'.join(
                                pattern))
                    cols = incompatible_columns(value)
                    print('{0:5}\t'.format('')+'\t'.join(cols))
                    print('---')

        return clr

    def irregular_patterns(self, accepted=2, mismatches=1, matches=1,
            debug=False, missing="Ø", irregular_prefix='!'):
        """
        Try to assign irregular patterns to accepted patterns.

        Parameters
        ----------
        accepted : int (default=2)
            Minimal size of clusters that we regard as regular.

        """
        bad_clusters = [(clr, pts[0]) for clr, pts in self.clusters.items() if len(pts)
                == 1]
        good_clusters = sorted([(clr, pts) for clr, pts in self.clusters.items() if len(pts)
            >= accepted], key=lambda x: len(x[1]), reverse=True)
        new_clusters = {clr: [] for clr, pts in good_clusters}
        irregular_patterns = []
        for clr, ptn in bad_clusters:
            if missing_in_pattern(ptn, missing=missing) <= 2:
                for clrB, pts in good_clusters:
                    match, mism = compatible_columns(clr, clrB)
                    if mism <= matches and match > matches:
                        new_clusters[clrB] += [clr]
                        irregular_patterns += [clr]
                        break
        # re-assign alignments to the data by adding the irregular character
        for key, value in sorted(new_clusters.items(), key=lambda x: len(x[1]),
                reverse=True):
            if len(value) > 0:
                if debug:
                    print('<<<')
                    print('{0:5}\t'.format(len(self.clusters[key]))+'\t'.join(key))
                for i, pattern in enumerate(value):
                    pt = []
                    for lid, (a, b) in enumerate(zip(key, pattern)):
                        if a != b and missing not in [a, b]:
                            pt += [irregular_prefix+b]
                            # assign pattern to the corresponding alignments
                            doc = self.cols[lid]
                            for cogid, position in self.clusters[pattern]:
                                if self._mode == 'fuzzy':
                                    word_indices = self.etd[self.ref][cogid][lid]
                                    if word_indices:
                                        for widx in word_indices:
                                            # get the position in the alignment
                                            alms = tokens2morphemes(self[widx,
                                                'alignment'])
                                            cog_pos = self[widx,
                                                    self.ref].index(cogid)
                                            new_alm = alms[cog_pos]
                                            new_alm[position] = '{0}{1}/{2}'.format(
                                                    irregular_prefix,
                                                    b,
                                                    a)
                                            alms[cog_pos] = new_alm
                                            self[widx, 'alignment'] = ' + '.join([' '.join(x) for x in alms]).split()
                                else:
                                    word_indices = self.etd[self.ref][cogid][lid]
                                    if word_indices:
                                        for widx in word_indices:
                                            alm = self[widx,
                                                    'alignment'].split()
                                            alm[position] = '{0}{1}/{2}'.format(
                                                    irregular_prefix,
                                                    b,
                                                    a)
                                            self[widx, 'alignment'] = ' '.join(alm)
                        else:
                            pt += [b]
                    if debug:
                        print('     \t'+'\t'.join(pt))
                if debug:
                    print('---')
        self.ipatterns = new_clusters
        for pattern, data in [(a, b) for a, b in bad_clusters if a not in 
                irregular_patterns]:
            cogid, position = data
            if self._mode == 'fuzzy':
                for indices in [idx for idx in self.etd[self.ref][cogid] if
                        idx]:
                    for widx in indices:
                        cog_pos = self[widx, self.ref].index(cogid)
                        #print(widx, self[widx, 'alignment'])
                        alms = tokens2morphemes(self[widx, 'alignment'])
                        new_alm = alms[cog_pos]
                        new_alm[position] = '{0}{1}'.format(
                                irregular_prefix, new_alm[position])
                        alms[cog_pos] = new_alm
                        self[widx, 'alignment'] = ' + '.join([
                            ' '.join(x) for x in alms]).split()

        return new_clusters

    def regular_cognates(self, regularity_threshold=0.5, alignment_threshold=0.5,
            correct_cognates=None, ref='cogid', score_mode='coverage'):
        """Evaluate the regularity of cognate sets.

        Notes
        -----
        The regularity is evaluated by pulling out all cognate sets in the data
        and checking how many of their columns are subject to irregularity,
        defined as being below the regularity threshold supplied by the user. For each
        alignment, these irregular cognate sets are summarized, and below the
        alignment threshold (a value between 0 and 1). It requires that
        clusters have been computed. If the users wants to correct those
        cognates by splitting them, this can be defined with a keyword.
        """
        cogs = defaultdict(list)
        if not hasattr(self, 'clusters'):
            raise ValueError('You should computer the clusters first.')
        for key, values in self.clusters.items():
            # assemble patterns
            patterns = [self.patterns[cog, pos] for cog, pos in values]
            regular = 1
            if score_patterns(patterns, mode=score_mode) <= regularity_threshold:
                regular = 0
            for cogid, _ in values:
                cogs[cogid] += [regular] # [score_patterns(patterns, mode=score_mode)]

        regulars = {}
        for cogid, scores in cogs.items():
            regulars[cogid] = sum(scores) / len(scores)
        irregulars = [cogid for cogid in regulars if regulars[cogid] <
                alignment_threshold]
        if correct_cognates:
            assert isinstance(correct_cognates, str)
            D = {}
            new_cogid = max(self.get_etymdict(ref))+1
            for idx, cogid in self.iter_rows(ref):
                if cogid in irregulars:
                    D[idx] = new_cogid
                    new_cogid += 1
                else:
                    D[idx] = cogid
            self.add_entries(correct_cognates, D, lambda x: x)
        return regulars, irregulars

    def load_patterns(self, patterns='patterns', alignment='alignment', 
            structure='structure', missing='Ø'):
        self.id2ptn = OrderedDict()
        self.clusters = OrderedDict()
        self.id2pos = defaultdict(set)
        self.sites = OrderedDict()
        # get the template
        template = [missing for m in self.taxa]
        tidx = {self.cols[i]: i for i in range(self.width)}
        for idx, ptn, alm, struc, doc, cogs in self.iter_rows(patterns, alignment,
                structure, 'doculect', self._ref):
            if self._mode == 'fuzzy':
                ptn = bt.lists(ptn)
                for i in range(len(alm.n)):
                    for j, (p, a) in enumerate(zip(ptn.n[i], alm.n[i])):
                        if not p == '0/n':
                            this_pattern = self.id2ptn.get(p, [t for t in template])
                            if this_pattern[tidx[doc]] == 'Ø':
                                this_pattern[tidx[doc]] = a
                            self.id2ptn[p] = this_pattern
                            self.id2pos[p].add((cogs[i], j))
            else:
                for j, (p, a) in enumerate(zip(ptn, alm)):
                    if not p == '0/n':
                        this_pattern = self.id2ptn.get(p, [t for t in
                            template])
                        if this_pattern[tidx[doc]] == 'Ø':
                            this_pattern[tidx[doc]] = a
                        self.id2ptn[p] = this_pattern
                        self.id2pos[p].add((cogs, j))
                        
        self.ptn2id = {tuple(v): k for k, v in self.id2ptn.items()}
        for k, v in self.id2ptn.items():
            self.clusters[tuple(v)] = list(self.id2pos[k])
            self.id2pos[k] = list(self.id2pos[k])
            for s in self.id2pos[k]:
                self.sites[s] = [(len(self.id2pos[k]), tuple(v))]

    def add_patterns(self, ref="patterns", irregular_patterns=False,
            proto=False):
        """Assign patterns to a new column in the word list.
        """
        if not hasattr(self, 'id2ptn'):
            self.id2ptn = {}
        if not hasattr(self, 'pattern2id'):
            self.ptn2id = {}
        if proto:
            pidx = self.cols.index(proto)
        else:
            pidx = 0

        if irregular_patterns:
            new_clusters = defaultdict(list)
            for reg, iregs in self.ipatterns.items():
                for cogid, position in self.clusters[reg]:
                    new_clusters[reg] += [(cogid, position)]
                for ireg in iregs:
                    for cogid, position in self.clusters[ireg]:
                        new_clusters[reg] += [(cogid, position)]
        else:
            new_clusters = self.clusters
        for pattern, rest in self.clusters.items():
            for cogid, position in rest:
                if (cogid, position) not in new_clusters[pattern]:
                    new_clusters[pattern] += [(cogid, position)]

        P = {idx: bt.lists(['0/n' if x not in  rc('morpheme_separators') else
            '+' for x in self[idx, 'alignment']]) for idx in self}
        for i, (pattern, data) in enumerate(sorted(new_clusters.items(),
            key=lambda x: len(x), reverse=True)):
            pattern_id = '{0}-{1}/{2}'.format(
                    i+1,
                    len(self.clusters[pattern]),
                    pattern[pidx]
                    )
            self.id2ptn[pattern_id] = pattern
            self.ptn2id[pattern] = pattern_id
            for cogid, position in data:
                word_indices = [c for c in self.etd[self.ref][cogid] if c]
                for idxs in word_indices:
                    for idx in idxs:
                        if self._mode == 'fuzzy':
                            #split_patterns = P[idx].n #, cldf=True)
                            pattern_position = self[idx, self.ref].index(cogid)
                            this_pattern = P[idx].n[pattern_position]
                            this_pattern[position] = pattern_id
                            P[idx].change(pattern_position, this_pattern)
                            #split_patterns[pattern_position] = this_pattern
                            ## here we want to insert the Xth element XXX
                            #pattern_position = self[idx, self.ref].index(cogid)
                            #P[idx].n[self[idx, self.ref].index(cogid)]

                            #P[idx] = bt.lists(' + '.join([' '.join(ptn) for ptn in
                            #    split_patterns]))
                        else:
                            P[idx][position] = pattern_id
        self.add_entries(ref, P, lambda x: x)

    def write_patterns(self, filename, proto=False, irregular_patterns=False):
        if proto:
            pidx = self.cols.index(proto)
        else:
            pidx = 0

        if irregular_patterns:
            new_clusters = defaultdict(list)
            for reg, iregs in self.ipatterns.items():
                for cogid, position in self.clusters[reg]:
                    new_clusters[reg] += [(cogid, position)]
                for ireg in iregs:
                    ireg_ = list(ireg)
                    for i, (a, b) in enumerate(zip(reg, ireg)):
                        if a != b and b != 'Ø':
                            ireg_[i] = a+'/'+b
                    ireg_ = tuple(ireg_)
                    for cogid, position in self.clusters[ireg]:
                        new_clusters[ireg_] += [(cogid, position)]
        else:
            new_clusters = self.clusters
        for pattern, rest in self.clusters.items():
            for cogid, position in rest:
                if (cogid, position) not in new_clusters[pattern]:
                    new_clusters[pattern] += [(cogid, position)]
        text = 'ID\tFREQUENCY\t{0}\t{1}\tCOGNATESETS\tCONCEPTS\n'.format(
                self.cols[pidx],
                '\t'.join([c for c in self.cols if c != self.cols[pidx]]))
                
        sound = ''
        idx = 0
        for pattern, entries in sorted(new_clusters.items(), key=lambda x:
                (x[0][pidx], len(x[1])), reverse=True):
            if sound != pattern[pidx]:
                sound = pattern[pidx]
                idx = 0
            concepts = []
            for x, y in entries:
                for entry in self.etd[self.ref][x]:
                    if entry:
                        for value in entry:
                            concepts += [self[value, 'concept']]
            concepts = ' / '.join(sorted(set(concepts)))

            idx += 1
            text += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                    idx, len(entries), pattern[pidx], '\t'.join([
                        p for i, p in enumerate(pattern) if i != pidx]),
                    ', '.join(['{0}:{1}'.format(x, y) for x, y in entries]),
                    concepts)
        with open(filename, 'w') as f:
            f.write(text)

    def predict_words(self, **kw):
        """
        Predict patterns for those cognate sets where we have missing data.
        """
        if not hasattr(self, 'sites'):
            raise ValueError('You need to compute alignment sites first')

        minrefs = kw.get('minrefs', self._minrefs)

        preds = {}
        count = 1
        for cogid, msa in self.msa[self._ref].items():
            missing = [t for t in self.cols if t not in msa['taxa']]
            if len(set(msa['taxa'])) >= minrefs:
                words = [bt.strings('') for m in missing]
                for i, m in enumerate(missing):
                    tidx = self.cols.index(m)
                    for j in range(len(msa['alignment'][0])):
                        try:
                            words[i] += [self.sites[cogid, j][0][1][tidx]]
                        except KeyError:
                            words[i] += ['?']
                if words:
                    preds[cogid] = dict(zip(missing, words))
                    if kw.get('debug'):
                        for m, w in zip(missing, words):
                            #print(count, cogid, self[msa['ID'][0], 'concept'], m, ' '.join(w))
                            count += 1
        return preds

    def stats(self, threshold=2):
        sums = sum([len(x) for x in self.clusters.values() if len(x) >=
            threshold])
        oris = sum([len(x) for x in self.clusters.values()])
        return sums / oris


                
                

