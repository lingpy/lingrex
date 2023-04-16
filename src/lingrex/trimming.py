"""
Functionality to trim alignments by removing sites.
"""
import random
import typing
import functools
import itertools
import collections

from lingpy.sequence.sound_classes import token2class

from lingrex.util import subsequence_of

__all__ = ['GAP', 'Site', 'Sites', 'prep_alignments']
GAP = '-'


class Site(list):
    """
    A site in an alignment is a "column", i.e. a list of the n-th sound in the aligned words.
    """
    def gap_ratio(self, gap: str = GAP) -> float:
        return self.count(gap) / len(self)

    def first_sound(self, gap=GAP):
        for s in itertools.dropwhile(lambda c: c == gap, self):
            return s

    def soundclass(self, gap: str = GAP) -> str:
        return token2class(self.first_sound(gap=gap) or "+", "cv")


class Sites(list):
    """
    A Sites object represents an alignment in the orthogonal view, i.e. listing columns rather
    than rows.
    """
    def __init__(self, alms: typing.List[typing.Union[str, typing.List[str]]], gap: str = GAP):
        """
        :parameter alms: List of aligned sequences.
        :parameter gap: String that codes gaps in alignment sites.
        """
        self.gap = gap
        super().__init__(Site([row[i] for row in alms]) for i in range(len(alms[0])))

    @property
    def gap_ratios(self) -> typing.List[float]:
        return [s.gap_ratio(gap=self.gap) for s in self]

    @property
    def soundclasses(self) -> typing.List[str]:
        return [s.soundclass(gap=self.gap) for s in self]

    def trimmed(self, idxs: typing.Iterable[int]) -> 'Sites':
        """
        Trim by removing the sites specified by index in `idxs`.
        """
        for idx in sorted(set(idxs), reverse=True):
            del self[idx]
        return self

    def to_alignment(self) -> typing.List[typing.List[str]]:
        return [[s[i] for s in self] for i in range(len(self[0]))]

    def trimmed_by_candidates(self,
                              candidate_indices: typing.List[int],
                              skeletons: typing.Iterable[str] = ("CV", "VC"),
                              exclude="_+") -> 'Sites':
        """
        Trim by removing sites specified by index in `candidate_indices` as long as this leaves an
        alignment containing at least one of the cv-patterns from `skeletons`.
        """
        skeleton = list(enumerate(self.soundclasses))
        idxs = {i for i, c in skeleton if c in exclude}  # Exclude markers.
        for idx in candidate_indices:
            current_skeleton = [c for i, c in skeleton if i not in idxs | {idx}]
            if any(subsequence_of(s, current_skeleton) for s in skeletons):
                # Trimming this site leaves a "big enough" remainder.
                idxs.add(idx)
            else:
                break
        return self.trimmed(idxs)

    def trimmed_by_gap(self,
                       threshold: float = 0.5,
                       skeletons: typing.Iterable[str] = ("CV", "VC"),
                       exclude="_+") -> 'Sites':
        """
        Trim alignment sites by gaps.

        :parameter threshold: Threshold for gap ratio by which sites should be trimmed.
        :param skeletons: Iterable of syllable-skeletons at least one of which should be preserved \
        for further processing.
        :param exclude: Sequence of strings that should be excluded from further processing,
            e.g. morpheme boundaries. Defaults to '_+'.
        """
        candidates = [  # Indices of sites with big enough gap ratio ordered by decreasing ratio.
            idx for idx, score
            in sorted(enumerate(self.gap_ratios), key=lambda x: x[1], reverse=True)
            if score >= threshold]
        return self.trimmed_by_candidates(candidates, skeletons=skeletons, exclude=exclude)

    def trimmed_by_core(self,
                        threshold: float = 0.5,
                        skeletons: typing.Iterable[str] = ("CV", "VC"),
                        exclude="_+") -> 'Sites':
        """
        Trim alignment sites by gaps, preserving a core of sites.

        :parameter threshold: Threshold by which sites with gaps should be trimmed.
        :param skeletons: Tuple of syllable-skeletons that should be preserved
            for further processing. Defaults to '("CV", "VC")'.
        :parameter gap: String that codes gaps in alignment sites. Defaults to '-'.
        """
        gap_or_not = [self.gap if ratio >= threshold else "S" for ratio in self.gap_ratios]
        takewhile_gap = functools.partial(itertools.takewhile, lambda c: c[1] == self.gap)
        leading_gaps = [i for i, _ in takewhile_gap(enumerate(gap_or_not))]
        trailing_gaps = [
            len(gap_or_not) - 1 - i for i, _ in takewhile_gap(enumerate(reversed(gap_or_not)))]
        return self.trimmed_by_candidates(
            trailing_gaps + leading_gaps, skeletons=skeletons, exclude=exclude)

    def trimmed_random(self,
                       method: typing.Union[str, typing.Callable] = 'trimmed_by_gap',
                       threshold: float = 0.5,
                       skeletons: typing.Iterable[str] = ("CV", "VC"),
                       exclude="_+") -> 'Sites':
        """
        For a base trim function, return a random version with a similar CV distribution.

        :parameter method: Trimming function that should be applied to compute the CV distribution.\
        Specified as name of a suitable method of `Sites`, or as callable.
        :parameter threshold: Threshold by which sites with gaps should be trimmed.
        :param skeletons: Tuple of syllable-skeletons that should be preserved
            for further processing. Defaults to '("CV", "VC")'.
        """
        if isinstance(method, str):
            method = getattr(Sites, method)
        reference_skeleton = method(
            Sites(self.to_alignment(), gap=self.gap),
            threshold=threshold,
            skeletons=skeletons,
            exclude=exclude).soundclasses
        # create a freq dict of ref skel
        rs_freqs = collections.Counter(reference_skeleton)
        # get a dictionary of indices by position
        indices = {  # soundclass mapped to list of indices in cv template.
            sc: [i[0] for i in items] for sc, items in itertools.groupby(
                sorted(enumerate(self.soundclasses), key=lambda ii: ii[1]),
                lambda ii: ii[1])}
        # random sample indices to be retained
        retain = [random.sample(indices[c], rs_freqs[c]) for c, _ in rs_freqs.items()]
        retain = set(itertools.chain(*retain))
        return self.trimmed([i for i in range(len(self)) if i not in retain])


def prep_alignments(aligned_wl, skeletons=("CV", "VC"), ref="cogid"):
    """"
    Preparing the alignments assures that the structure is correctly
    added to the wordlist.

    :param wordlist: A lingpy Alignments.
    :type wordlist: :class:lingpy.Alignments
    :param skeletons: Tuple of syllable-skeletons that should be preserved
        for further processing. Defaults to '("CV", "VC")'.
    :type skeletons: tuple
    :param ref: The column which stores the cognate sets, defaults to 'cogid'
    :type ref: str
    :return: Pre-processed alignments.
    :rtype: :class:lingpy.Alignments
    """
    whitelist = []
    for _, msa in aligned_wl.msa[ref].items():
        skel = Sites(msa["alignment"]).soundclasses
        if any([subsequence_of(s, skel) for s in skeletons]):
            whitelist += msa["ID"]
    aligned_wl.add_entries(
        "structure", "tokens", lambda x: " ".join(Sites([c for c in x]).soundclasses))
    dct = {0: aligned_wl.columns}
    for idx in whitelist:
        dct[idx] = aligned_wl[idx]
    return aligned_wl.__class__(dct, transcription="form")
