import math
import pathlib
import itertools
import collections

from lingpy.sequence.sound_classes import class2tokens
from lingpy.settings import rc
from lingpy.align.sca import get_consensus, Alignments
from lingpy.util import pb
from lingpy import log
from lingpy import basictypes as bt

import networkx as nx


def consensus_pattern(patterns, missing="Ø"):
    """
    Return consensus pattern of multiple patterns.

    :param patterns: list of patterns
    :param missing: the character used to represent missing values

    .. note:: This consensus method raises an error if the patterns contain incompatible
       columns (non-identical values apart from the missing data character in the same
       column).
    """
    out = []
    for i in range(len(patterns[0])):
        col = [line[i] for line in patterns]
        no_gaps = [x for x in col if x != missing]
        if len(set(no_gaps)) > 1:
            raise ValueError("Your patterns are incompatible")
        out += [no_gaps[0] if no_gaps else missing]
    return tuple(out)


def incompatible_columns(patterns, missing="Ø"):
    """
    Compute whether a pattern has incompatible columns.
    """
    columns = []
    for i in range(len(patterns[0])):
        col = [
            patterns[j][i] for j in range(len(patterns)) if patterns[j][i] != missing
        ]
        columns.append("*" if len(set(col)) > 1 else "")
    return columns


def score_patterns(patterns, missing="Ø", mode="coverage"):
    """
    Function gives a score for the overall number of reflexes.

    .. note:: This score tells simply to which degree a pattern is filled. It divides the
       number of cells not containing missing data by the number of cells in the
       matrix.
    """
    # return -1 if the patterns are not compatible
    for i in range(len(patterns[0])):
        if len(set([row[i] for row in patterns if row[i] != missing])) > 1:
            return -1
    if len(patterns) <= 1:
        return -1

    if mode not in ["ranked", "pairs", "squared", "coverage"]:
        raise ValueError("you must select an appropriate mode")

    # we rank the columns by sorting them first
    if mode == "ranked":
        cols = []
        for i in range(len(patterns[0])):
            cols += [sum([0 if row[i] == missing else 1 for row in patterns])]
        # sort the columns
        ranks, cols = list(range(1, len(cols) + 1))[::-1], sorted(cols, reverse=True)
        scores = []
        for rank, col in zip(ranks, cols):
            scores += [rank * col]
        return sum(scores) / sum(ranks) / len(patterns)

    if mode == "squared":
        psize = len(patterns[0])
        scores = [((psize - row.count(missing)) / psize) ** 2 for row in patterns]
        return sum(scores) / len(scores)

    if mode == "pairs":

        # count the number of pairs in the data
        pairs = 0
        covered = 0
        m, n = len(patterns[0]), len(patterns)
        for i in range(n):
            vals = m - patterns[i].count(missing)
            pairs += (vals ** 2 - vals) / 2
        for i in range(m):
            vals = n - [p[i] for p in patterns].count(missing)
            pairs += (vals ** 2 - vals) / 2
            if vals != 0:
                covered += 1
        return ((pairs / n) / covered) / m

    if mode == "coverage":
        cols = []
        for i in range(len(patterns[0])):
            col = [row[i] for row in patterns]
            cols += [len(patterns) - col.count(missing)]
        return (sum(cols) / len(patterns[0])) / len(patterns)  # * len(patterns[0]))


def compatible_columns(colA, colB, missing="Ø", gap="-"):
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
        if missing not in [a, b]:
            if a != b:
                mismatches += 1
            else:
                if a != gap:
                    matches += 1
    return matches, mismatches


def density(wordlist, ref="cogid"):
    """Compute the density of a wordlist.


    XXX TOTO: better move to util or elsewhere

    Notes
    -----
    We define the density of a wordlist by measuring how many words can be
    explained by the same cognate set.
    """
    scores = []
    for concept in wordlist.rows:
        idxs = wordlist.get_list(row=concept, flat=True)
        cogids = [wordlist[idx, ref] for idx in idxs]
        sums = [1 / cogids.count(cogid) for idx, cogid in zip(idxs, cogids)]
        scores.append(sum(sums) / len(sums))
    return 1 - sum(scores) / len(scores)


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

    def __init__(
        self,
        wordlist,
        minrefs=3,
        ref="cogids",
        structure="structure",
        missing="Ø",
        gap="-",
        irregular="!?",
        **keywords
    ):
        Alignments.__init__(self, wordlist, ref=ref, **keywords)
        self.ref = ref
        self._structure = structure
        self.minrefs = minrefs
        self.missing = missing
        self.gap = gap
        self.irregular = irregular
        if structure not in self.columns:
            raise ValueError("no column {0} for structure was found".format(structure))

    def positions_from_prostrings(self, cogid, indices, alignment, structures):
        """
        Return positions matching from an alignment and user-defined prosodic strings
        """
        if self._mode == "fuzzy":
            strucs = []
            for idx, struc, alm in zip(indices, structures, alignment):
                pos_ = self[idx, self._ref].index(cogid)
                strucs += [class2tokens(struc.n[pos_], alm)]
        else:
            strucs = [
                class2tokens(struc, alm) for struc, alm in zip(structures, alignment)
            ]
        get_consensus(alignment, gaps=True)
        prostring = []
        for i in range(len(strucs[0])):
            row = [x[i] for x in strucs if x[i] != "-"]
            prostring += [row[0] if row else "+"]
        return [(i, p) for i, p in enumerate(prostring)]

    def reflexes_from_pos(
        self, position, taxa, current_taxa, alignment, missing, irregular
    ):
        reflexes = []
        for t in taxa:
            if t not in current_taxa:
                reflexes += [missing]
            else:
                reflex = alignment[current_taxa.index(t)][position]
                if "/" in reflex:
                    reflex = reflex.split("/")[1]
                elif reflex[0] in irregular:
                    reflex = missing
                reflexes += [reflex]
        return reflexes

    def _check(self):
        """
        Check for problematic patterns in the data.
        """
        errors = []
        for idx, struc, alm in self.iter_rows(self._structure, self._alignment):
            self[idx, self._structure] = self._str_type(struc)
            self[idx, self._alignment] = self._str_type(alm)
            if not len(self[idx, self._structure]) == len(
                [x for x in self[idx, self._alignment] if x != "-"]
            ):
                print(
                    idx,
                    self[idx, self._structure],
                    "|",
                    self[idx, self._alignment],
                    "|",
                    self[idx, "tokens"],
                )
                log.warning("alignment and structure do not match in {0}".format(idx))
                errors += [idx]
        return errors

    def get_sites(self):
        """
        Retrieve the alignment sites of interest for initial analysis.
        """
        sites, all_sites, taxa = (
            collections.OrderedDict(),
            collections.OrderedDict(),
            self.cols,
        )
        errors = self._check()
        if errors:
            raise ValueError("found {0} problems in the data".format(len(errors)))

        # iterate over all sites in the alignment
        visited = []
        for cogid, msa in pb(
            sorted(self.msa[self.ref].items()),
            desc="CoPaR: get_patterns()",
            total=len(self.msa[self.ref]),
        ):
            # get essential data: taxa, alignment, etc.
            _taxa = [t for t in taxa if t in msa["taxa"]]
            _idxs = {t: msa["taxa"].index(t) for t in _taxa}
            _alms = [msa["alignment"][_idxs[t]] for t in _taxa]
            _wlid = [msa["ID"][_idxs[t]] for t in _taxa]

            # store visited entries
            visited += msa["ID"]
            if len(_taxa) >= self.minrefs:
                if self._mode == "fuzzy":
                    _strucs = []
                    for _widx in _wlid:
                        _these_strucs = self[_widx, self._structure]
                        _strucs += [_these_strucs]
                else:
                    _strucs = [self[idx, self._structure] for idx in _wlid]
                positions = self.positions_from_prostrings(cogid, _wlid, _alms, _strucs)
                for pidx, pos in positions:
                    reflexes = self.reflexes_from_pos(
                        pidx, taxa, _taxa, _alms, self.missing, self.irregular
                    )
                    sites[cogid, pidx] = [pos, tuple(reflexes)]
            for pidx in range(len(_alms[0])):
                reflexes = self.reflexes_from_pos(
                    pidx, taxa, _taxa, _alms, self.missing, self.irregular
                )
                all_sites[cogid, pidx] = reflexes

        # add non-visited segments
        for idx in [i for i in self if i not in visited]:
            if self._mode == "fuzzy":
                for tt, ss, cogid in zip(
                    self[idx, self._segments].n,
                    self[idx, self._structure].n,
                    self[idx, self._ref],
                ):
                    for i, (t, s) in enumerate(zip(tt, ss)):
                        all_sites[cogid, i] = [
                            self.missing if tax != self[idx][self._colIdx] else t
                            for tax in self.cols
                        ]
            else:
                for i, (t, s) in enumerate(
                    zip(self[idx, self._segments], self[idx, self._structure])
                ):
                    all_sites[self[idx, self.ref], i] = [
                        self.missing if tax != self[idx][self._colIdx] else t
                        for tax in self.cols
                    ]

        self.sites = sites
        self.all_sites = all_sites

    def cluster_sites(self, match_threshold=1, score_mode="pairs"):
        """Cluster alignment sites using greedy clique cover.
        :param match_threshold: The threshold of matches for accepting two
            compatible columns.
        :param score_mode: select between "pairs", "coverage"

        .. note:: This algorithm follows the spirit of the Welsh-Powell algorithm for
           graph coloring. Since graph coloring is the inverse of clique
           partitioning, we can use the algorithm in the same spirit.

        """
        if not hasattr(self, "clusters"):
            self.clusters = collections.defaultdict(list)
            for (cogid, idx), (pos, ptn) in self.sites.items():
                self.clusters[pos, ptn] += [(cogid, idx)]
        clusters = self.clusters
        while True:
            prog = 0
            with pb(
                desc="CoPaR: cluster_sites()", total=len(self.clusters)
            ) as progress:
                sorted_clusters = sorted(
                    clusters.items(),
                    key=lambda x: (
                        score_patterns(
                            [self.sites[y][1] for y in x[1]], mode=score_mode
                        ),
                        len(x[1]),
                    ),
                    reverse=True,
                )
                out = []
                while sorted_clusters:
                    ((this_pos, this_cluster), these_vals), remaining_clusters = (
                        sorted_clusters[0],
                        sorted_clusters[1:],
                    )
                    queue = []
                    for (next_pos, next_cluster), next_vals in remaining_clusters:
                        match, mism = compatible_columns(
                            this_cluster,
                            next_cluster,
                            missing=self.missing,
                            gap=self.gap,
                        )
                        if (
                            this_pos == next_pos
                            and match >= match_threshold  # noqa: W503
                            and mism == 0  # noqa: W503
                        ):
                            this_cluster = consensus_pattern(
                                [this_cluster, next_cluster]
                            )
                            these_vals += next_vals
                        else:
                            queue += [((next_pos, next_cluster), next_vals)]
                    sorted_clusters = queue
                    out += [((this_pos, this_cluster), these_vals)]
                    progress.update(len(self.sites) - len(queue) - prog)
                    prog = len(self.sites) - len(queue)
                clusters = {tuple(a): b for a, b in out}
                alls = [c for c in clusters]
                match = 0
                for i, (_a, a) in enumerate(alls):
                    for j, (_b, b) in enumerate(alls):
                        if i < j and _a == _b:
                            ma, mi = compatible_columns(
                                a, b, missing=self.missing, gap=self.gap
                            )
                            if ma and not mi:
                                match += 1
                if not match:
                    break
                else:
                    log.warning(
                        "iterating, since {0} clusters can further be merged".format(
                            match
                        )
                    )
        self.clusters = clusters
        self.ordered_clusters = sorted(clusters, key=lambda x: len(x[1]))

    def sites_to_pattern(self, threshold=1):
        """Algorithm assigns alignment sites to patterns.

        Notes
        -----
        We rank according to general compatibility.
        """
        asites = collections.defaultdict(list)
        for consensus in pb(
            self.clusters, desc="CoPaR: sites_to_pattern()", total=len(self.clusters)
        ):
            sites = self.clusters[consensus]
            for cog, pos in sites:
                struc, pattern = self.sites[cog, pos]
                for strucB, consensusB in self.clusters:
                    ma, mi = compatible_columns(pattern, consensusB)
                    if struc == strucB and not mi and ma >= threshold:
                        asites[cog, pos] += [(ma, struc, consensusB)]
        self.patterns = asites

    def fuzziness(self):
        return sum([len(b) for a, b in self.patterns.items()]) / len(self.patterns)

    def irregular_patterns(self, accepted=2, matches=1, irregular_prefix="!"):
        """
        Try to assign irregular patterns to accepted patterns.

        Parameters
        ----------
        accepted : int (default=2)
            Minimal size of clusters that we regard as regular.

        """
        bad_clusters = [
            (clr, pts[0]) for clr, pts in self.clusters.items() if len(pts) == 1
        ]
        good_clusters = sorted(
            [(clr, pts) for clr, pts in self.clusters.items() if len(pts) >= accepted],
            key=lambda x: len(x[1]),
            reverse=True,
        )
        new_clusters = {clr: [] for clr, pts in good_clusters}
        irregular_patterns = []
        for clr, ptn in bad_clusters:
            if ptn.count(self.missing) <= 2:
                for clrB, pts in good_clusters:
                    match, mism = compatible_columns(clr[1], clrB[1])
                    if mism <= matches and match > matches:
                        new_clusters[clrB] += [clr]
                        irregular_patterns += [clr]
                        break
        # re-assign alignments to the data by adding the irregular character
        for key, value in sorted(
            new_clusters.items(), key=lambda x: len(x[1]), reverse=True
        ):
            if len(value) > 0:
                for i, pattern in enumerate(value):
                    pt = []
                    for lid, (a, b) in enumerate(zip(key[1], pattern[1])):
                        if a != b and self.missing not in [a, b]:
                            pt += [irregular_prefix + b]
                            # assign pattern to the corresponding alignments
                            for cogid, position in self.clusters[pattern]:
                                if self._mode == "fuzzy":
                                    word_indices = self.etd[self.ref][cogid][lid]
                                    if word_indices:
                                        for widx in word_indices:
                                            # get the position in the alignment
                                            alms = self[widx, self._alignment].n
                                            cog_pos = self[widx, self.ref].index(cogid)
                                            new_alm = alms[cog_pos]
                                            new_alm[position] = "{0}{1}/{2}".format(
                                                irregular_prefix, b, a
                                            )
                                            alms[cog_pos] = new_alm
                                            self[
                                                widx, self._alignment
                                            ] = self._str_type(
                                                " + ".join(
                                                    [" ".join(x) for x in alms]
                                                ).split()
                                            )
                                else:
                                    word_indices = self.etd[self.ref][cogid][lid]
                                    if word_indices:
                                        for widx in word_indices:
                                            alm = self._str_type(
                                                self[widx, self._alignment]
                                            )
                                            alm[position] = "{0}{1}/{2}".format(
                                                irregular_prefix, b, a
                                            )
                                            self[
                                                widx, self._alignment
                                            ] = self._str_type(" ".join(alm))
                        else:
                            pt += [b]

        self.ipatterns = new_clusters
        for pattern, data in [
            (a, b) for a, b in bad_clusters if a not in irregular_patterns
        ]:
            cogid, position = data
            if self._mode == "fuzzy":
                for indices in [idx for idx in self.etd[self.ref][cogid] if idx]:
                    for widx in indices:
                        cog_pos = self[widx, self.ref].index(cogid)
                        alms = self[widx, self._alignment].n
                        new_alm = alms[cog_pos]
                        new_alm[position] = "{0}{1}".format(
                            irregular_prefix, new_alm[position]
                        )
                        alms[cog_pos] = new_alm
                        self[widx, self._alignment] = self._str_type(
                            " + ".join([" ".join(x) for x in alms]).split()
                        )

        return new_clusters

    def load_patterns(self, patterns="patterns"):
        self.id2ptn = collections.OrderedDict()
        self.clusters = collections.OrderedDict()
        self.id2pos = collections.defaultdict(set)
        self.sites = collections.OrderedDict()
        # get the template
        template = [self.missing for m in self.cols]
        tidx = {self.cols[i]: i for i in range(self.width)}
        for idx, ptn, alm, struc, doc, cogs in self.iter_rows(
            patterns, self._alignment, self._structure, "doculect", self._ref
        ):
            if self._mode == "fuzzy":
                ptn = bt.lists(ptn)
                for i in range(len(alm.n)):
                    for j, (p, a) in enumerate(zip(ptn.n[i], alm.n[i])):
                        if not p == "0/n":
                            this_pattern = self.id2ptn.get(p, [t for t in template])
                            if this_pattern[tidx[doc]] == "Ø":
                                this_pattern[tidx[doc]] = a
                            self.id2ptn[p] = this_pattern
                            self.id2pos[p].add((cogs[i], j))
            else:
                for j, (p, a) in enumerate(zip(ptn, alm)):
                    if not p == "0/n":
                        this_pattern = self.id2ptn.get(p, [t for t in template])
                        if this_pattern[tidx[doc]] == "Ø":
                            this_pattern[tidx[doc]] = a
                        self.id2ptn[p] = this_pattern
                        self.id2pos[p].add((cogs, j))

        self.ptn2id = {tuple(v): k for k, v in self.id2ptn.items()}
        for k, v in self.id2ptn.items():
            self.clusters[tuple(v)] = list(self.id2pos[k])
            self.id2pos[k] = list(self.id2pos[k])
            for s in self.id2pos[k]:
                self.sites[s] = [(len(self.id2pos[k]), tuple(v))]

    def add_patterns(
        self, ref="patterns", irregular_patterns=False, proto=False, override=True
    ):
        """Assign patterns to a new column in the word list."""
        if not hasattr(self, "id2ptn"):
            self.id2ptn = {}
        if not hasattr(self, "pattern2id"):
            self.ptn2id = {}
        if proto:
            pidx = self.cols.index(proto)
        else:
            pidx = 0

        if irregular_patterns:
            new_clusters = collections.defaultdict(list)
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

        P = {
            idx: bt.lists(
                [
                    "0/n" if x not in rc("morpheme_separators") else "+"
                    for x in self[idx, self._alignment]
                ]
            )
            for idx in self
        }
        for i, ((struc, pattern), data) in enumerate(
            sorted(new_clusters.items(), key=lambda x: len(x), reverse=True)
        ):
            pattern_id = "{0}-{1}/{2}".format(
                i + 1, len(self.clusters[struc, pattern]), pattern[pidx]
            )
            self.id2ptn[pattern_id] = pattern
            self.ptn2id[pattern] = pattern_id
            for cogid, position in data:
                word_indices = [c for c in self.etd[self.ref][cogid] if c]
                for idxs in word_indices:
                    for idx in idxs:
                        if self._mode == "fuzzy":
                            pattern_position = self[idx, self.ref].index(cogid)
                            this_pattern = P[idx].n[pattern_position]
                            try:
                                this_pattern[position] = pattern_id
                                P[idx].change(pattern_position, this_pattern)
                            except:  # noqa: E722
                                log.warning("error in {0}".format(cogid))

                        else:
                            P[idx][position] = pattern_id
        self.add_entries(ref, P, lambda x: x, override=override)

    def write_patterns(self, filename, proto=False, irregular_patterns=False):
        if proto:
            pidx = self.cols.index(proto)
        else:
            pidx = 0

        if not hasattr(self, "id2ptn"):
            raise ValueError("You should run CoPaR.add_patterns first!")

        if irregular_patterns:
            new_clusters = collections.defaultdict(list)
            for (pos, reg), iregs in self.ipatterns.items():
                for cogid, position in self.clusters[pos, reg]:
                    new_clusters[pos, reg] += [(cogid, position)]
                for _, ireg in iregs:
                    ireg_ = list(ireg)
                    print(ireg_)
                    for i, (a, b) in enumerate(zip(reg, ireg)):
                        print(i, a, b)
                        if a != b and b != self.missing:
                            ireg_[i] = a + "/" + b
                    ireg_ = tuple(ireg_)
                    self.ptn2id[ireg_] = self.ptn2id[reg]
                    for cogid, position in self.clusters[pos, ireg]:
                        new_clusters[pos, ireg_] += [(cogid, position)]
        else:
            new_clusters = self.clusters
        for (struc, pattern), rest in self.clusters.items():
            for cogid, position in rest:
                if (cogid, position) not in new_clusters[struc, pattern]:
                    new_clusters[struc, pattern] += [(cogid, position)]
        text = "ID\tSTRUCTURE\tFREQUENCY\t{0}\t{1}\tCOGNATESETS\tCONCEPTS\n".format(
            self.cols[pidx], "\t".join([c for c in self.cols if c != self.cols[pidx]])
        )

        sound = ""
        idx = 0
        for (struc, pattern), entries in sorted(
            new_clusters.items(),
            key=lambda x: (x[0][0], x[0][1][pidx], len(x[1])),
            reverse=True,
        ):
            if sound != pattern[pidx]:
                sound = pattern[pidx]
                idx = 0
            concepts = []
            for x, y in entries:
                for entry in self.etd[self.ref][x]:
                    if entry:
                        for value in entry:
                            concepts += [self[value, "concept"]]
            concepts = " / ".join(sorted(set(concepts)))

            idx += 1
            text += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(
                self.ptn2id[pattern].split("/")[0],
                struc,
                len(entries),
                pattern[pidx],
                "\t".join([p for i, p in enumerate(pattern) if i != pidx]),
                ", ".join(["{0}:{1}".format(x, y) for x, y in entries]),
                concepts,
            )
        pathlib.Path(filename).write_text(text, encoding="utf8")

    def purity(self):
        """
        Compute the purity of the cluster analysis.

        .. note:: The purity is here interpreted as the degree to which
           patterns are filled with non-missing values. In this sense, it
           indicates to which degree information is computed and to which
           degree information is already provided by the data itself.
        """

        def get_purity(patterns):
            all_sums = []
            for i in range(len(patterns[0])):
                col = [line[i] for line in patterns]
                subset = set(col)
                sums = []
                for itm in subset:
                    if itm != self.missing:
                        sums += [col.count(itm) ** 2]
                if sums:
                    sums = math.sqrt(sum(sums)) / len(col)
                else:
                    sums = 0
                all_sums += [sums]
            return sum(all_sums) / len(all_sums)

        graph = self.get_cluster_graph()
        purities = []
        for node, data in graph.nodes(data=True):
            patterns = []
            for neighbor in graph[node]:
                patterns += [graph.nodes[neighbor]["pattern"].split()]
            if patterns:
                purities += [get_purity(patterns)]
            else:
                purities += [0]
        return sum(purities) / len(purities)

    def get_cluster_graph(self):
        """
        Compute a graph of the clusters.

        .. note:: In the cluster graph, the sites in the alignments are the
           nodes and the edges are drawn between nodes assigned to the same
           pattern.
        """

        graph = nx.Graph()
        for (pos, ptn), sites in self.clusters.items():
            for site in sites:
                graph.add_node(
                    "{0[0]}-{0[1]}".format(site),
                    pattern=" ".join(ptn),
                    site=" ".join(self.sites[site][1]),
                )

        for ((s1, p1), ptn1), ((s2, p2), ptn2) in itertools.combinations(
            self.sites.items(), r=2
        ):
            if ptn1[0] == ptn2[0]:
                m, mm = compatible_columns(ptn1[1], ptn2[1])
                if m and not mm:
                    graph.add_edge("{0}-{1}".format(s1, p1), "{0}-{1}".format(s2, p2))
        return graph

    def upper_bound(self):
        """
        Compute upper bound for clique partitioning following Bhasker 1991.
        """
        degs = {s: 0 for s in self.sites}
        sings = {s: 0 for s in self.sites}
        for (nA, (posA, ptnA)), (nB, (posB, ptnB)) in itertools.combinations(
            self.sites.items(), r=2
        ):
            if posA == posB:
                m, n = compatible_columns(ptnA, ptnB)
                if n > 0:
                    degs[nA] += 1
                    degs[nB] += 1
                else:
                    sings[nA] += 1
                    sings[nB] += 1
            else:
                degs[nA] += 1
                degs[nB] += 1

        return max([b for a, b in degs.items() if sings[a] > 0])

    def predict_words(self, **kw):
        """
        Predict patterns for those cognate sets where we have missing data.

        .. note::

           Purity (one of the return values) measures how well a given sound
           for a given site is reflected by one single sound (rather than
           multiple patterns pointing to different sounds) for a given
           doculect. It may be seen as a control case for the purity of a given
           prediction: if there are many alternative possibilities, this means
           that there is more uncertainty regarding the reconstructions or
           predictions.

        """
        if not hasattr(self, "sites"):
            raise ValueError("You need to compute alignment sites first")

        minrefs = self.minrefs
        missing = self.missing
        samples = kw.get("samples", 3)

        # pre-analyse the data to get for each site the best patterns in ranked
        # form
        ranked_sites = {}
        ranked_clusters = sorted(
            [(s, p, len(f)) for (s, p), f in self.clusters.items()],
            key=lambda x: x[2],
            reverse=True,
        )
        for (cogid, pos), ptns in self.patterns.items():
            struc, ptn = self.sites[cogid, pos]
            missings = [i for i in range(self.width) if ptn[i] == missing]
            if (struc, ptn) in self.clusters:
                ranked_sites[cogid, pos] = [
                    (len(self.clusters[struc, ptn]), struc, ptn)
                ]
            else:
                ranked_sites[cogid, pos] = [(1, struc, ptn)]
            for strucB, ptnB, freq in ranked_clusters:
                m, mm = compatible_columns(ptn, ptnB)
                if struc == strucB and m >= 1 and mm == 0:
                    if len(missings) > len(
                        [ptnB[i] for i in missings if ptnB[i] == missing]
                    ):
                        ranked_sites[cogid, pos] += [(freq, strucB, ptnB)]

        purity = {site: {} for site in ranked_sites}

        preds = {}
        for cogid, msa in self.msa[self._ref].items():
            missings = [t for t in self.cols if t not in msa["taxa"]]
            if len(set(msa["taxa"])) >= minrefs:
                words = [bt.strings("") for m in missings]
                for i, m in enumerate(missings):
                    tidx = self.cols.index(m)
                    for j in range(len(msa["alignment"][0])):
                        segments = collections.defaultdict(int)
                        sidx = 0
                        if (cogid, j) in ranked_sites:
                            while True:
                                this_segment = ranked_sites[cogid, j][sidx][2][tidx]
                                score = ranked_sites[cogid, j][sidx][0]
                                if this_segment != missing:
                                    segments[this_segment] += score
                                sidx += 1
                                if sidx == len(ranked_sites[cogid, j]):
                                    break

                        if not (cogid, j) in purity:
                            purity[cogid, j] = {}

                        if not segments:
                            words[i] += ["Ø"]
                            purity[cogid, j][m] = 0
                        else:
                            purity[cogid, j][m] = math.sqrt(
                                sum(
                                    [
                                        (s / sum(segments.values())) ** 2
                                        for s in segments.values()
                                    ]
                                )
                            )
                            words[i] += [
                                "|".join(
                                    [
                                        s
                                        for s in sorted(
                                            segments,
                                            key=lambda x: segments[x],
                                            reverse=True,
                                        )
                                    ][:samples]
                                )
                            ]
                if words:
                    preds[cogid] = dict(zip(missings, words))

        pudity = {doc: [] for doc in self.cols}
        for site, docs in purity.items():
            for doc in docs:
                pudity[doc] += [purity[site][doc]]
        for doc, purs in pudity.items():
            if purs:
                pudity[doc] = sum(purs) / len(purs)
            else:
                pudity[doc] = 0

        return preds, purity, pudity
