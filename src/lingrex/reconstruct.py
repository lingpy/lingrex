"""
Module provides methods for linguistic reconstruction.
"""
import itertools
import collections

from lingpy.align.sca import Alignments, get_consensus
from lingpy.sequence.sound_classes import prosodic_string, class2tokens
from lingpy.align.multiple import Multiple
from lingpy.align.pairwise import edit_dist, nw_align
from lingpy.evaluate.acd import _get_bcubed_score as get_bcubed_score
from lingpy.align.sca import normalize_alignment
import networkx as nx
from networkx.algorithms.clique import find_cliques
from lingpy import log

from lingrex.util import clean_sound, ungap, alm2tok


class CorPaRClassifier(object):

    def __init__(self, minrefs=2, missing=0, threshold=1):
        self.G = nx.Graph()
        self.missing = 0
        self.threshold = threshold

    def compatible(self, ptA, ptB):
        """
        Check for compatibility of two patterns.
        """
        res = {True: 0, False: 0}
        for a, b in zip(ptA, ptB):
            if a and b:
                res[a == b] += 1
        return res[True], res[False]

    def consensus(self, nodes):
        """
        Create a consensus pattern of multiple alignment sites.
        """
        cons = []
        for i in range(len(nodes[0])):
            nocons = True
            for node in nodes:
                if node[i] != self.missing:
                    cons += [node[i]]
                    nocons = False
                    break
            if nocons:
                cons += [self.missing]
        return tuple(cons)

    def fit(self, X, y):
        """
        Train the prediction of data in y with data in X.

        :param X: Two-dimensional array with observations.
        :param y: One-dimensional array with results.
        """
        # get identical patterns
        P = collections.defaultdict(list)
        for i, row in enumerate(X):
            P[tuple(row + [y[i]])] += [i]
        # make graph
        for (pA, vA), (pB, vB) in itertools.combinations(P.items(), r=2):
            match_, mismatch = self.compatible(pA, pB)
            if not mismatch and match_ >= self.threshold:
                if pA not in self.G:
                    self.G.add_node(pA, freq=len(vA))
                if pB not in self.G:
                    self.G.add_node(pB, freq=len(vB))
                self.G.add_edge(pA, pB, weight=match_)
        self.patterns = collections.defaultdict(collections.Counter)
        self.lookup = collections.defaultdict(collections.Counter)
        # get cliques
        for nodes in find_cliques(self.G):
            cons = self.consensus(list(nodes))
            self.patterns[cons[:-1]][cons[-1]] = len(nodes)
            for node in nodes:
                self.lookup[node[:-1]][cons[:-1]] += len(nodes)
        self.predictions = {
            ptn: counts.most_common(1)[0][0] for ptn, counts in self.patterns.items()}
        for ptn, counts in self.lookup.items():
            self.predictions[ptn] = self.predictions[counts.most_common(1)[0][0]]

        # make index of data points for quick search based on attested data
        self.ptnlkp = collections.defaultdict(list)
        for ptn in self.patterns:
            for i in range(len(ptn)):
                if ptn[i] != self.missing:
                    self.ptnlkp[i, ptn[i]] += [ptn]

    def predict(self, matrix):
        out = []
        for row in matrix:
            ptn = tuple(row)
            if ptn in self.predictions:
                out.append(self.predictions[ptn])
            else:
                candidates = collections.Counter()
                for i in range(len(ptn) - 1):
                    if ptn[i] != self.missing:
                        for ptnB in self.ptnlkp[i, ptn[i]]:
                            if ptnB not in candidates:
                                match_, mismatch = self.compatible(ptn, ptnB)
                                if match_ and not mismatch:
                                    candidates[ptnB] = match_ + len(ptn)
                                elif match_ - mismatch:
                                    candidates[ptnB] = match_ - mismatch
                if candidates:
                    self.predictions[tuple(row)] = self.predictions[candidates.most_common(1)[0][0]]
                    out += [self.predictions[tuple(row)]]
                else:
                    out += [self.missing]
        return out


class ReconstructionBase(Alignments):
    """
    Basic class for the phonological reconstruction.
    """
    def __init__(
            self, infile, target=None, ref="cogids", fuzzy=True,
            transcription="form", missing="Ø", gap="-"):
        Alignments.__init__(self, infile, fuzzy=fuzzy, ref=ref, transcription=transcription)
        self.target = target
        self.missing = missing
        self.gap = gap
        self.languages = [t for t in self.cols if t != target]
        self.target = target
        self.tgtidx = self.cols.index(target)
        self.lngidx = {t: self.cols.index(t) for t in self.languages}

    def iter_sequences(self, aligned=False):
        """
        Iterate over aligned or unaligned sequences with or without the target \
                sequence.
        """
        seq_ref = self._alignments if aligned else self._segments
        for cogid, idxs in self.etd[self._ref].items():
            if idxs[self.tgtidx]:
                if self._mode == "fuzzy":
                    target = self[idxs[self.tgtidx][0], seq_ref].n[
                        self[idxs[self.tgtidx][0], self._ref].index(cogid)]
                else:
                    target = self[idxs[self.tgtidx][0], seq_ref]
                alignment, languages = [], []
                for j, lng in enumerate(self.languages):
                    lidx = self.lngidx[lng]
                    if idxs[lidx]:
                        languages += [lng]
                        idx = idxs[lidx][0]
                        if self._mode == "fuzzy":
                            alm = self[idx, seq_ref].n[self[idx, self._ref].index(cogid)]
                        else:
                            alm = self[idx, seq_ref]
                        alignment.append([clean_sound(x) for x in alm])
                alignment.append([clean_sound(x) for x in target])
                if aligned:
                    alignment = normalize_alignment(alignment)
                languages.append(self.target)
                yield cogid, alignment, languages


class OneHot(object):
    """
    Create a one-hot-encoder from a matrix.
    """

    def __init__(self, matrix):
        self.vals = []
        for i in range(len(matrix[0])):
            cols = [row[i] for row in matrix]
            self.vals += [sorted(set(cols)) + ["?"]]

    def __call__(self, matrix):
        out = [[] for row in matrix]
        for i, vals in enumerate(self.vals):
            for j in range(len(matrix)):
                template = [0 for k in vals]
                try:
                    template[matrix[j][i]] = 1
                except IndexError:
                    template[-1] = 1
                out[j] += template
        return out


def transform_alignment(seqs,
                        languages,
                        all_languages,
                        align=True,
                        training=True,
                        missing="Ø",
                        gap="-",
                        startend=False,
                        prosody=False,
                        position=False,
                        firstlast=False):
    """
    Basic alignment function used for phonological reconstruction.
    """
    if align:
        seqs = [[s for s in seq if s != gap] for seq in seqs]
        msa = Multiple([[s for s in seq if s != gap] for seq in seqs])
        msa.prog_align()
        alms = [alm for alm in msa.alm_matrix]
    else:
        seqs = [[s for s in seq if s != gap] for seq in seqs]
        alms = normalize_alignment([s for s in seqs])
    if training:
        alms = ungap(alms, languages, languages[-1])
        these_seqs = seqs[:-1]
    else:
        these_seqs = seqs
    matrix = [[missing for x in all_languages] for y in alms[0]]
    for i in range(len(alms[0])):
        for j, lng in enumerate(languages):
            lidx = all_languages.index(lng)
            snd = clean_sound(alms[j][i])
            matrix[i][lidx] = snd
    if position:
        for i in range(len(matrix)):
            matrix[i] += [i]
    if startend:
        matrix[0] += [0]
        for i in range(1, len(matrix) - 1):
            matrix[i] += [1]
        if len(matrix) > 1:
            matrix[-1] += [2]
    if prosody:
        for i, c in enumerate(
                get_consensus(
                    [class2tokens(prosodic_string(seqs[j], _output="CcV"), alms[j])
                     for j in range(len(these_seqs))],
                    gaps=True)):
            matrix[i] += [c]
    if firstlast:
        if training:
            all_seqs = len(all_languages) - 1
        else:
            all_seqs = len(all_languages)
        for i, row in enumerate(matrix):
            for j in range(all_seqs):
                matrix[i] += [matrix[0][j], matrix[-1][j]]

    # for debugging
    for row in matrix:
        assert len(row) == len(matrix[0])
    return matrix


class PatternReconstructor(ReconstructionBase):
    """
    Automatic reconstruction with correspondence patterns.
    """

    def fit(self, clf=None, onehot=False, func=None, aligned=False):
        """
        Fit a classifier to the data.

        :param clf: a classifier with a predict function.
        """
        self.patterns = collections.defaultdict(lambda: collections.defaultdict(list))
        self.occurrences = collections.defaultdict(list)
        self.func = func or transform_alignment

        for cogid, alignment, languages in self.iter_sequences():
            if len(alignment) >= 2:
                matrix = self.func(
                    alignment,
                    languages,
                    self.languages + [self.target],
                    training=True)
                for i, row in enumerate(matrix):
                    ptn = tuple(row[:len(self.languages)] + row[len(self.languages) + 1:])
                    self.patterns[ptn][row[len(self.languages)]] += [
                        (cogid, i)]
                    for j, lng in enumerate(self.languages):
                        if row[j] not in [self.missing]:
                            self.occurrences[lng, j, row[j]] += [(cogid, i)]
                    for j in range(len(self.languages) + 1, len(row)):
                        self.occurrences["feature-{0}".format(j - 1), j - 1, row[j]] += [(cogid, i)]

        self.snd2idx = {(i, self.missing): 0 for i in range(len(matrix[0]))}
        for i in range(len(matrix[0])):
            self.snd2idx[i, self.gap] = 1

        idxtracker = {i: 2 for i in range(len(matrix[0]))}
        for lng, lidx, sound in self.occurrences:
            last_idx = idxtracker[lidx]
            if (lidx, sound) not in self.snd2idx:
                self.snd2idx[lidx, sound] = last_idx
                idxtracker[lidx] += 1

        self.tgt2idx = {}
        idx = 1
        for pattern in self.patterns:
            for sound in self.patterns[pattern]:
                if sound not in self.tgt2idx:
                    self.tgt2idx[sound] = idx
                    idx += 1

        self.matrix = []
        self.solutions = []
        for pattern, sounds in self.patterns.items():
            for sound, vals in sounds.items():
                tidx = self.tgt2idx[sound]
                row = []
                for i in range(len(pattern)):
                    sidx = self.snd2idx[i, pattern[i]]
                    row += [sidx]
                for cogid, idx in vals:
                    self.matrix += [row]
                    self.solutions += [tidx]
        self.dim = len(self.matrix[0])
        if clf is not None:
            self.clf = clf
        else:
            self.clf = CorPaRClassifier()
        log.info("fitting classifier")
        if onehot:
            self.onehot = OneHot(self.matrix)
            self.clf.fit(self.onehot(self.matrix), self.solutions)
        else:
            self.clf.fit(self.matrix, self.solutions)
        self.idx2tgt = {v: k for k, v in self.tgt2idx.items()}
        log.info("fitted the classifier")

    def predict(
            self, alignment, languages, unknown="?", onehot=False,
            desegment=True):
        """
        Predict a word form from an alignment.

        :param desegment: Return the form without gaps and ungapped tokens.
        """
        matrix = self.func(alignment, languages, self.languages, training=False)
        for row in matrix:
            assert len(row) == self.dim
        new_matrix = [[0 for char in row] for row in matrix]
        for i, row in enumerate(matrix):
            for j, char in enumerate(row):
                new_matrix[i][j] = self.snd2idx.get((j, char), 0)
        if hasattr(self, "onehot"):
            new_matrix = self.onehot(new_matrix)
        out = [self.idx2tgt.get(idx, unknown) for idx in self.clf.predict(new_matrix)]
        return alm2tok(out) if desegment else out


def eval_by_dist(data, func=None, **kw):
    """
    Evaluate by measuring distances between sequences.

    :param data: List of tuples with prediction and attested sequence.
    :param func: Alignment function (defaults to edit distance)

    :note: Defaults to the unnormalized edit distance.
    """
    func = func or edit_dist
    scores = []
    for seqA, seqB in data:
        if not seqA:
            seqA = ["?"]
        if not seqB:
            seqB = ["?"]
        scores += [func(seqA, seqB, **kw)]
    return sum(scores) / len(scores)


def eval_by_bcubes(data, func=None, **kw):
    """
    Evaluate by measuring B-Cubed F-scores.

    :param data: List of tuples with prediction and attested sequence.
    :param func: Alignment function (defaults to Needleman-Wunsch)
    """
    numsA, numsB = {"": 0}, {"": 0}
    func = func or nw_align
    almsA, almsB = [], []
    for seqA, seqB in data:
        if not seqA:
            seqA = ["?"]
        if not seqB:
            seqB = ["?"]
        almA, almB, score = func(seqA, seqB, **kw)
        for a, b in zip(almA, almB):
            if a not in numsA:
                numsA[a] = max(numsA.values()) + 1
            if b not in numsB:
                numsB[b] = max(numsB.values()) + 1
            almsA += [numsA[a]]
            almsB += [numsB[b]]
    p, r = get_bcubed_score(almsA, almsB), get_bcubed_score(almsB, almsA)
    return 2 * (p * r) / (p + r)
