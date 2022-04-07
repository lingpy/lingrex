"""
Basic code for borrowing detection.
"""
import itertools
import collections

from lingpy import Pairwise
from lingpy.compare.partial import Partial
from lingpy.compare.lexstat import LexStat

import networkx as nx

from lingpy.util import pb


def internal_cognates(
    wordlist,
    family="family",
    partial=True,
    method="lexstat",
    runs=10000,
    threshold=0.50,
    smooth=1,
    ratio=(2, 1),
    vscale=0.5,
    restricted_chars="_",
    modes=[("global", -1, 0.5), ("overlap", -1, 0.5)],
    ref="autocogids",
    cluster_method="upgma",
    model="sca",
):
    """
    Cluster the data into cognate sets, but only inside each family.

    :param family: name of the column in which language family information can
       be found (defaults="family")
    """
    families = {wordlist[k, family] for k in wordlist}

    # split data into parts
    D = {k: {} for k in sorted(families)}
    for idx, fam in wordlist.iter_rows(family):
        D[fam][idx] = [cell for cell in wordlist[idx]]

    gcogid = 0
    G = {}
    for fam, data in D.items():
        data[0] = [h for h in wordlist.columns]
        if partial:
            lex = Partial(data, model=model)
            if method == "lexstat":
                lex.get_partial_scorer(
                    runs=runs,
                    smooth=smooth,
                    ratio=ratio,
                    vscale=vscale,
                    restricted_chars=restricted_chars,
                    modes=modes,
                )
            lex.partial_cluster(
                ref=ref,
                method=method,
                cluster_method=cluster_method,
                threshold=threshold,
            )
        else:
            lex = LexStat(data, model=model)
            if method == "lexstat":
                lex.get_scorer(
                    runs=runs,
                    smooth=smooth,
                    ratio=ratio,
                    vscale=vscale,
                    restricted_chars=restricted_chars,
                    modes=modes,
                )
            lex.cluster(
                ref=ref,
                method=method,
                cluster_method=cluster_method,
                threshold=threshold,
            )

        # prepare global cognate indicies
        if partial:
            C = {idx: len(lex[idx, ref]) * [0] for idx in lex}
            etd = lex.get_etymdict(ref=ref)
            for cogid, idxs in etd.items():
                for idx_ in idxs:
                    if idx_:
                        for idx in idx_:
                            cogids = lex[idx, ref]
                            C[idx][cogids.index(cogid)] = cogid + gcogid
        else:
            C = {idx: 0 for idx in lex}
            etd = lex.get_etymdict(ref=ref)
            for cogid, idxs in etd.items():
                for idx_ in idxs:
                    if idx_:
                        for idx in idx_:
                            C[idx] = cogid + gcogid
        for idx in lex:
            G[idx] = C[idx]
        gcogid += max(etd) + 1

    renumber = {}
    cogid = 1
    if partial:
        for idx, vals in G.items():
            f = wordlist[idx, family]
            new_cogids = []
            for v in vals:
                if (f, v) in renumber:
                    new_cogids += [renumber[f, v]]
                else:
                    renumber[f, v] = cogid
                    new_cogids += [cogid]
                    cogid += 1
            G[idx] = new_cogids
    else:
        for idx, val in G.items():
            f = wordlist[idx, family]
            if (f, val) not in renumber:
                renumber[f, val] = cogid
                cogid += 1
            G[idx] = renumber[f, val]

    wordlist.add_entries(ref, G, lambda x: x)


def external_cognates(
    wordlist,
    cognates="autocogid",
    ref="autoborid",
    threshold=0.3,
    segments="tokens",
    gop=-1,
    family="family",
    doculect="doculect",
    concept="concept",
    align_mode="overlap",
):
    """
    Compute language-external cognates and assign them to cognate sets.

    :param cognates: The column which holds previously calculated cognates.
    :param ref: The column which will store the new borrowing identifiers.
    :param family: The column storing family information.
    :param doculect: The column storing doculect information.
    """

    B = {}
    borid = 1
    # iterate over the concepts
    for concept in pb(wordlist.rows):
        idxs = wordlist.get_list(row=concept, flat=True)
        for idx in idxs:
            B[idx] = 0
        taxa = [wordlist[idx, doculect] for idx in idxs]
        famis = [wordlist[idx, family] for idx in idxs]
        if len(set(famis)) > 1:
            G = nx.Graph()
            tokens = [wordlist[idx, segments] for idx in idxs]
            cogids = [wordlist[idx, cognates] for idx in idxs]

            # assemble cogids to groups
            groups = collections.defaultdict(list)
            for i, d, t, c in zip(idxs, taxa, tokens, cogids):
                groups[c] += [(i, d, t)]

            for group, items in groups.items():
                G.add_node(
                    str(group),
                    concept=concept,
                    taxa=", ".join([t[1] for t in items]),
                    idxs=", ".join([str(t[0]) for t in items]),
                    family=wordlist[[t[0] for t in items][0], family],
                )

            # compare groups
            for (gA, iA), (gB, iB) in itertools.combinations(list(groups.items()), r=2):
                if G.nodes[str(gA)]["family"] != G.nodes[str(gB)]["family"]:
                    wpairs = [(a[2], b[2]) for a, b in itertools.product(iA, iB)]

                    pairs = Pairwise(wpairs)
                    pairs.align(distance=True, gop=gop, mode=align_mode)
                    dst = []
                    for i, p in enumerate(pairs._alignments):
                        dst += [p[2]]

                    dst = sum(dst) / len(dst)
                    if dst <= threshold:
                        G.add_edge(str(gA), str(gB), distance=dst)

            # components
            for i, comp in enumerate(nx.connected_components(G)):
                if len(comp) > 1:
                    table = []
                    for cogid in comp:
                        idxs = [int(x) for x in G.nodes[cogid]["idxs"].split(", ")]
                        for idx in idxs:
                            table += [
                                [
                                    wordlist[idx, doculect],
                                    wordlist[idx, concept],
                                    str(wordlist[idx, segments]),
                                    wordlist[idx, family],
                                    cogid,
                                ]
                            ]
                        for idx in idxs:
                            B[idx] = borid
                    borid += 1
    wordlist.add_entries(ref, B, lambda x: x)
