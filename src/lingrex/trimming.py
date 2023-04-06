"""
Trimming functionalities in lingrex.
"""
import random
import collections
from lingpy.sequence.sound_classes import token2class


def revert(alms):
    return [[row[i] for row in alms] for i in range(len(alms[0]))]        


def get_skeleton(alms, gap="-"):
    return [token2class(([c for c in col if c != gap] or ["+"])[0], "cv") for col in revert(alms)]


def apply_trim(alms, idxs):
    """
    Basic trimming function, based on a selection of indices that are trimmed.
    """
    return [[row[i] for i in range(len(row)) if i not in idxs] for row in alms]


def subsequence_of(source, target):
    """
    Check if source is a subsequence of target.
    """
    q_1, q_2 = list(target), list(source)
    while q_1:
        s = q_1.pop(0)
        if q_2 and q_2[0] == s:
            q_2.pop(0)
        elif q_2:
            pass
        else:
            break
    if not q_2:
        return True
    return False


def consecutive_gaps(seq, gap="-"):
    """
    Return consecutive gaps in line.
    """
    start, end = [0], [len(seq)]
    for i in range(len(seq)):
        if seq[i] == gap:
            gapped = True
        else:
            gapped = False
        if gapped:
            start += [i + 1]
        if not gapped:
            break
    for i in range(len(seq) - 1, 0, -1):
        if seq[i] == "-":
            gapped = True
        else:
            gapped = False
        if gapped:
            end += [i]
        if not gapped:
            break
    return start[:-1], end[::-1][:-1]


def gap_profile(alms, gap="-"):
    """
    Return a profile of the gap-ration per column.
    """
    return [col.count(gap) / len(col) for col in revert(alms)]


def trim_by_gap(alms, threshold=0.5, skeletons=("CV", "VC"), gap="-", exclude="_+"):
    """
    Trim alignment sites by gaps.

    :parameter alms: Alignment sites of a cognate set.
    :type alms: list
    :parameter threshold: Threshold by which sites with gaps should be trimmed.
        Defaults to '0.5'.
    :type threshold: int
    :param skeletons: Tuple of syllable-skeletons that should be preserved
        for further processing. Defaults to '("CV", "VC")'.
    :type skeletons: tuple
    :parameter gap: String that codes gaps in alignment sites. Defaults to '-'.
    :type gap: string
    :param exclude: Sequence of strings that should be excluded from further processing,
        e.g. morpheme boundaries. Defaults to '_+'.
    :param exclude: str
    :return: Indices of trimmed strings.
    :rtype: set
    """
    skeleton = get_skeleton(alms, gap=gap)
    profile = gap_profile(alms, gap)
    # exclude markers
    idxs = [i for i, c in enumerate(skeleton) if c in exclude]
    # order by gap weights
    sorted_profile = sorted(enumerate(profile), key=lambda x: x[1], reverse=True)
    while sorted_profile:
        idx, score = sorted_profile.pop(0)
        if score >= threshold:
            current_skeleton = "".join([skeleton[i] for i in
                                        range(len(skeleton)) if i not in idxs+[idx]])
            if any([subsequence_of(s, current_skeleton) for s in skeletons]):
                idxs += [idx]
            else:
                break
        else:
            break
    return sorted(set(idxs))


def trim_by_core(
        alms, threshold=0.5, skeletons=("CV", "VC"), gap="-",
        exclude="_+"
        ):
    """
    Trim alignment sites by gaps, preserving a core of sites.

    :parameter alms: Alignment sites of a cognate set.
    :type alms: list
    :parameter threshold: Threshold by which sites with gaps should be trimmed.
        Defaults to '0.5'.
    :type threshold: int
    :param skeletons: Tuple of syllable-skeletons that should be preserved
        for further processing. Defaults to '("CV", "VC")'.
    :type skeletons: tuple
    :parameter gap: String that codes gaps in alignment sites. Defaults to '-'.
    :type gap: string
    :return: Indices of trimmed strings.
    :rtype: set
    """
    cons = [gap if c >= threshold else "S" for c in gap_profile(alms, gap=gap)]
    skeleton = get_skeleton(alms, gap=gap)
    left, right = consecutive_gaps(cons)

    idxs = [i for i, c in enumerate(skeleton) if c in exclude]

    # order by first
    sorted_indices = right[::-1] + left
    while sorted_indices:
        idx = sorted_indices.pop(0)
        current_skeleton = "".join([skeleton[i] for i in
                                    range(len(skeleton)) if i not in idxs+[idx]])
        if any([subsequence_of(s, current_skeleton) for s in skeletons]):
            idxs += [idx]
        else:
            break
    return sorted(set(idxs))


def trim_random(
        alms,
        func=None,
        threshold=0.5,
        skeletons=("CV", "VC"),
        gap="-",
        exclude="_+",
        ):
    """
    For a base trim function, return a random version with a similar CV distribution.

    :parameter alms: Alignment sites of a cognate set.
    :type alms: list
    :parameter func: Trimming function that should be applied. Defaults to 'None'.
    :type func: function
    :parameter threshold: Threshold by which sites with gaps should be trimmed.
        Defaults to '0.5'.
    :type threshold: int
    :param skeletons: Tuple of syllable-skeletons that should be preserved
        for further processing. Defaults to '("CV", "VC")'.
    :type skeletons: tuple
    :parameter gap: String that codes gaps in alignment sites. Defaults to '-'.
    :type gap: string
    :return: Indices of trimmed strings.
    :rtype: set
    """
    func = func or trim_by_gap
    reference = apply_trim(
            alms,
            func(
                alms, threshold=threshold, skeletons=skeletons, gap=gap,
                exclude=exclude)
            )
    reference_skeleton = get_skeleton(reference)
    # create a freq dict of ref skel
    rs_freqs = collections.defaultdict(int)
    for c in reference_skeleton:
        rs_freqs[c] += 1
    # get a dictionary of indices by position
    indices = collections.defaultdict(list)
    for i, c in enumerate(get_skeleton(alms)):
        indices[c] += [i]
    # random sample indices to be retained
    retain = []
    for c, _ in rs_freqs.items():
        retain += random.sample(indices[c], rs_freqs[c])
    to_trim = [i for i in range(len(alms[0])) if i not in retain]
    return to_trim
