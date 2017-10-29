# *-* coding: utf-8 *-*
from __future__ import print_function, division, unicode_literals
from six import text_type
from collections import defaultdict
from itertools import combinations
from functools import cmp_to_key

from lingpy.sequence.sound_classes import (
        prosodic_string, tokens2morphemes,
        class2tokens
        )
from lingpy.align.sca import get_consensus, SCA, Alignments
from lingpy.util import pb
from lingpy.compare.partial import Partial, _get_slices
from lingpy.read.qlc import normalize_alignment


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
                matches += 1
    return matches, mismatches

class CorPatR(Alignments):
    def __init__(self, wordlist, **keywords):
        Alignments.__init__(self, wordlist, **keywords)
    
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
                pos = self[idx, self._ref].index(cogid)
                strucs += [class2tokens(tokens2morphemes(struc)[pos], alm)]
        else:
            strucs = [class2tokens(struc.split(' '), alm) for struc, alm in
                    zip(structures, alignment)]
        consensus = get_consensus(alignment, gaps=True)
        prostring = []
        for i in range(len(strucs[0])):
            row = [x[i] for x in strucs if x[i] != '-']
            prostring += [row[0]]
        if len(prostring) != len(alignment[0]):
            print(prostring)
            print(consensus)
            print('\n'.join(['\t'.join(alm) for alm in alignmet]))
            input('problematic alignment')
        return [i for i, p in enumerate(prostring) if p == pos]

    def reflexes_from_pos(self, position, taxa, current_taxa, alignment,
            missing):
        reflexes = []
        for t in taxa:
            if t not in current_taxa:
                reflexes += [missing]
            else:
                reflexes += [alignment[current_taxa.index(t)][position]]
        return reflexes

    def get_patterns(
            self,
            ref='cogid',
            pos='A',
            minimal_taxa=3,
            taxa=None,
            missing='Ø',
            debug=False,
            prostring=False,
            irregular = "!",
            ):
        """
        """
        patterns, all_patterns = {}, {}
        taxa = taxa or self.cols
        
        # iterate over all sites in the alignment
        with pb(desc='COMPATIBILITY PATTERN EXTRACTION',
                total=len(self.msa[ref])) as progress:
            for cogid, msa in self.msa[ref].items():
                progress.update(1)
                # get essential data: taxa, alignment, etc.
                _taxa = [t for t in taxa if t in msa['taxa']]
                _idxs = {t: msa['taxa'].index(t) for t in _taxa}
                _alms = [msa['alignment'][_idxs[t]] for t in _taxa]
                _wlid = [msa['ID'][_idxs[t]] for t in _taxa]
                if len(_taxa) >= minimal_taxa:
                    if not prostring:
                        positions = self.positions_from_alignment(_alms, pos)
                    else:
                        _strucs = [self[idx, prostring] for idx in \
                                _wlid]
                        positions = self.positions_from_prostrings(
                                cogid, _wlid, _alms, _strucs, pos)
                    for pidx in positions:
                        reflexes = self.reflexes_from_pos(
                                pidx, taxa, _taxa, _alms, missing)
                        patterns[cogid, pidx] = reflexes
                for pidx in range(len(_alms[0])):
                    reflexes = self.reflexes_from_pos(
                            pidx, taxa, _taxa, _alms, missing)
                    all_patterns[cogid, pidx] = reflexes
        self.patterns = patterns
        self.all_patterns = all_patterns
        return patterns, all_patterns 

    def sort_patterns(self, patterns=None, missing='Ø', threshold=2, debug=False):
        def sorter(x, y):
            """Sorter arranges patterns to allow for a quicker assembly"""
            cc = compatible_columns(x[1], y[1])
            lA = len([i for i in x[1] if i == missing])
            lB = len([i for i in y[1] if i == missing])
            tA = x[1].count(missing)
            tB = y[1].count(missing)
            if cc == -1:
                return -1
            if cc == 0:
                if tA < tB:
                    return -1
                return 1
            if lA < lB:
                return -1
            elif lA > lB:
                return 1
            return 0
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
        
        single, mass = 0, 0
        for key, vals in out.items():
            if len(vals) >= threshold:
                if debug: print(str(len(vals))+ '\t ITMS\t'+ '\t'.join(key))
                print(vals)
                for v in vals:
                    mass += 1
                    if debug: print('---\t'+'\t'.join(
                        patterns[v])
                        )
                if debug: print('===')
            else:
                single += len(vals)
        self.clusters = out
        return out
        
    def cluster_patterns(self, clusters=None, threshold=2, debug=False):
        """Cluster patterns following an the spirit of Welsh-Powel algorithm."""
        clusters = clusters or self.clusters
        sorted_clusters = sorted(clusters.items(), key=lambda x: len(x[1]),
                reverse=True)
        out = []
        while sorted_clusters:
            (this_cluster, these_vals), remaining_clusters = sorted_clusters[0], sorted_clusters[1:]
            queue = []
            for next_cluster, next_vals in remaining_clusters:
                match, mism = compatible_columns(this_cluster, next_cluster)
                if match != 0 and mism == 0:
                    this_cluster = consensus_pattern([this_cluster,
                            next_cluster])
                    these_vals += next_vals
                else:
                    queue += [(next_cluster, next_vals)]
            out += [(this_cluster, these_vals)]
            sorted_clusters = sorted(queue, key=lambda x: len(x[1]),
                    reverse=True)
        clusters = {tuple(a): b for a, b in out}
        single, mass = 0, 0
        for key, vals in sorted(clusters.items(), key=lambda x: len(x[1])):
            if len(vals) >= threshold:
                if debug: print(str(len(vals))+ '\t'+ '\t'.join(key))
                for v in vals:
                    mass += 1
                    if debug: print('---\t'+'\t'.join(
                        self.patterns[v])
                        )
                if debug: print('===')
            else:
                single += len(vals)
        self.clusters = clusters
        print(single, mass)
