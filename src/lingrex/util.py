"""
Utility functions for the lingrex package.
"""
import math
import pathlib

from lingpy import tokens2class, prosodic_string
from lingpy.align.sca import get_consensus
from lingpy import basictypes as bt
from lingpy.sequence.ngrams import get_n_ngrams


def lingrex_path(*comps):
    return str(pathlib.Path(__file__).parent.joinpath(*comps))


def bleu_score(word, reference, n=4, weights=None, trim=True):
    """
    Compute the BLEU score for predicted word and reference.

    :param word: the predicted word
    :param reference: the predicted reference
    :param n: the order of ngrams
    :param weights: list of weights, should be the same size as n
    :param trim: bool, decide to trim n-grams or not
    """

    if not weights:
        weights = [1 / n for x in range(n)]

    scores = []
    for i in range(1, n + 1):
        new_wrd = list(get_n_ngrams(word, i))
        new_ref = list(get_n_ngrams(reference, i))
        if trim and i > 1:
            new_wrd = new_wrd[i - 1:-(i - 1)]
            new_ref = new_ref[i - 1:-(i - 1)]

        clipped, divide = [], []
        for itm in set(new_wrd):
            clipped += [new_ref.count(itm)]
            divide += [new_wrd.count(itm)]
        scores += [sum(clipped) / sum(divide)]

    # calculate arithmetic mean
    out_score = 1
    for weight, score in zip(weights, scores):
        out_score = out_score * (score ** weight)

    bp = 1 if len(word) > len(reference) else math.e ** (1 - (len(reference) / len(word)))
    return bp * (out_score ** (1 / sum(weights)))


def clean_sound(sound):
    """
    Get rid of "a/b" notation for sound segments.
    """
    return ".".join([s.split('/')[1] if "/" in s else s for s in sound.split('.')])


def alm2tok(seq, gap="-"):
    """Turn an alignment into a sequence."""
    return [clean_sound(x) for x in unjoin(seq) if x != gap]


def unjoin(seq):
    """
    Turn segments joined by a dot into unjoined segments.
    """
    out = []
    for itm in seq:
        out += itm.split('.')
    return out


def ungap(alignment, languages, proto):
    """
    Trim an MSA to remove all gaps in the target sequence.

    :examples:
      >>> ungap([['a', 'b'], ['x', '-'], ['y', '-']], ['proto', 'l1', 'l2'], 'proto')
      ... [['a.b'], ['x'], ['y']]
      >>> ungap([['a', 'b'], ['x', '-'], ['y', 'h']], ['proto', 'l1', 'l2'], 'proto')
      ... [['a', 'b'], ['x', '-'], ['y', 'h']]
    """
    pidxs = [i for i, taxon in enumerate(languages) if taxon == proto]
    merges = []
    for i in range(len(alignment[0])):  # go through the rows of the alignment ...
        col = [row[i] for row in alignment]
        # ... looking for gap-only alignments (in non-proto languages):
        if {site for j, site in enumerate(col) if j not in pidxs} == {'-'}:
            merges += [i]
    if not merges:
        return alignment
    new_alms = []
    for i, row in enumerate(alignment):
        new_alm, mergeit, started = [], False, True
        for j, cell in enumerate(row):
            if j in merges or mergeit:
                mergeit = False
                if not started:  # j != 0:
                    if cell != "-":
                        new_alm[-1] += '.' + cell if new_alm[-1] else cell
                else:
                    mergeit = True
                    new_alm.append("" if cell == "-" else cell)
            else:
                started = False
                new_alm.append(cell)
        new_alms.append([cell or "-" for cell in new_alm])
    return new_alms


def add_structure(
        wordlist, model="cv", segments="tokens", structure="structure", ref="cogid", gap="-"):
    """Add structure to a wordlist to make sure correspondence patterns can be
    inferred"""
    if model not in ["cv", "c", "CcV", "ps", "nogap"]:
        raise ValueError("[i] you need to select a valid model")
    D = {}
    if model == "cv":
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = " ".join(tokens2class(tks, "cv")).lower()

    if model == "c":
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = (
                " ".join(tokens2class(tks, "cv"))
                .lower()
                .replace("v", "c")
                .replace("t", "c")
            )
    if model == "nogap":
        assert hasattr(wordlist, "msa")
        for cogid, msa in wordlist.msa[ref].items():
            cons = [
                "c" if c != gap else gap
                for c in get_consensus(msa["alignment"], gaps=True)
            ]
            for idx, alm in zip(msa["ID"], msa["alignment"]):
                struc = []
                for a, b in zip(cons, alm):
                    if b != "-":
                        struc += [a]
                D[idx] = " ".join(struc)
        for idx, tks in wordlist.iter_rows(segments):
            if idx not in D:
                D[idx] = " ".join(["c" if c != "+" else c for c in tks])
    if model == "CcV":
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = " ".join(
                list(prosodic_string(tks, _output="CcV").replace("_", "+"))
            )
    if model == "ps":
        for idx, tks in wordlist.iter_rows(segments):
            D[idx] = " ".join(list(prosodic_string(tks)))

    if hasattr(wordlist, "_mode") and wordlist._mode == "fuzzy":
        struc_ = bt.lists
    else:
        struc_ = bt.strings
    wordlist.add_entries(structure, D, lambda x: struc_(x))
