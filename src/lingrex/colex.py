"""
Functions for partial colexification manipulations.
"""
import collections


def find_bad_internal_alignments(alignments, ref="cogids"):
    """
    Helper function discards wrongly assigned cross-semantic cognates.

    .. note:: The function essentially iterates over the alignments and picks
        out those in which the same language has the same cognate ID, and if
        the alignment itself differs, it assigns it a new cognate ID. It
        presupposes that the data has not been analyzed in search for
        cross-semantic cognates.
    """
    newIDs = {}

    def get_all_indices(lst):
        idxs = collections.defaultdict(list)
        for i, l in enumerate(lst):
            idxs[l] += [i]
        return idxs

    new_cogid = max(alignments.msa[ref]) + 1
    for cogid, msa in alignments.msa[ref].items():
        idxs = [i for t, i in get_all_indices(msa["taxa"]).items() if len(i) > 1]
        for idx in idxs:
            tups = [tuple(msa["alignment"][x]) for x in idx]
            if len(set(tups)) > 1:
                bestt = sorted(tups, key=lambda x: tups.count(x), reverse=True)[0]
                for x in idx:
                    if tuple(msa["alignment"][x]) != bestt:
                        newIDs[msa["ID"][x]] = (cogid, new_cogid)
                        new_cogid += 1

    for idx, (cogid, new_cogid) in newIDs.items():
        this_idx = alignments[idx, ref].index(cogid)
        alignments[idx, ref][this_idx] = new_cogid


def expand_alignment(msa, taxa, missing="Ø"):
    """
    Expand an alignment by adding a symbol for missing taxa.
    """
    out = []
    for taxon in taxa:
        if taxon in msa["taxa"]:
            tidx = msa["taxa"].index(taxon)
            out.append(
                [x.split("/")[1] if "/" in x else x for x in msa["alignment"][tidx]]
            )
        else:
            out.append(len(msa["alignment"][0]) * [missing])
    return out


def compatible(msa1, msa2, missing="Ø", gap="-"):
    """Compare two alignments and check whether they colexify."""
    matches = 0
    for line1, line2 in zip(msa1, msa2):
        if [x for x in line1 if x != gap] == [
            x for x in line2 if x != gap
        ] and missing not in line1 + line2:
            matches += 1
        else:
            if list(set(line1))[0] != missing and list(set(line2))[0] != missing:
                return False
    return matches


def merge_alignments(almA, almB, missing="Ø", gap="-"):
    """
    Merge two alignments which are compatible.
    """
    out = []
    missing_taxa = []
    for k, (a, b) in enumerate(zip(almA, almB)):
        if (
            len(set(a)) == 1
            and list(set(a))[0] == missing  # noqa: W503
            and len(set(b)) == 1  # noqa: W503
            and list(set(b))[0] == missing  # noqa: W503
        ):
            missing_taxa += [k]
    i, j = 0, 0
    while i < len(almA[0]) and j < len(almB[0]):
        colA, colB = [row[i] for row in almA], [row[j] for row in almB]
        if colA == colB:
            out += [colA]
            i += 1
            j += 1
        else:
            col = []
            for a, b in zip(colA, colB):
                if a == gap and a != b and b != missing:
                    ncol = []
                    for k, c in enumerate(colA):
                        if c == missing and k not in missing_taxa:
                            ncol += [gap]
                        else:
                            ncol += [c]
                    out += [ncol]
                    i += 1
                    col = []
                    break
                if b == gap and a != b and a != missing:
                    ncol = []
                    for k, c in enumerate(colB):
                        if c == missing and k not in missing_taxa:
                            ncol += [gap]
                        else:
                            ncol += [c]
                    out += [ncol]
                    j += 1
                    col = []
                    break

                col.append(b if a == missing else a)
            if col:
                out += [col]
                i += 1
                j += 1
    if i < len(almA[0]):
        ncol = []
        for k, c in enumerate([row[i] for row in almA]):
            if c == missing and k not in missing_taxa:
                ncol += [gap]
            else:
                ncol += [c]
        out += [ncol]
    elif j < len(almB[0]):
        ncol = []
        for k, c in enumerate([row[j] for row in almB]):
            if c == missing and k not in missing_taxa:
                ncol += [gap]
            else:
                ncol += [c]
        out += [ncol]

    nalm = []
    for i in range(len(out[0])):
        nalm += [[row[i] for row in out]]
    return nalm


def find_colexified_alignments(
    alignments, cognates="cogids", missing="Ø", ref="crossids"
):
    """
    Identify identical alignments in a dataset and label them as homophones.
    """

    queue = []
    for cogid, msa in sorted(
        alignments.msa[cognates].items(),
        key=lambda x: len(set(x[1]["taxa"])),
        reverse=True,
    ):
        queue += [(cogid, expand_alignment(msa, alignments.taxa, missing=missing))]

    merged = {}

    while queue:
        this_cogid, this_msa = queue.pop(0)
        deletes = []
        merged[this_cogid] = this_cogid
        for i, (other_cogid, other_msa) in enumerate(queue):
            if compatible(this_msa, other_msa) >= 1:
                this_msa = merge_alignments(this_msa, other_msa)
                merged[other_cogid] = this_cogid
                deletes += [i]

        for i in deletes[::-1]:
            del queue[i]

    # assemble the clusters now
    if alignments._mode == "fuzzy":
        alignments.add_entries(ref, cognates, lambda x: [merged.get(y, y) for y in x])
    else:
        alignments.add_entries(ref, cognates, lambda x: merged.get(x, x))
