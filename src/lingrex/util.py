"""
Utility functions for the lingrex package.
"""
import os
from lingpy.sequence.sound_classes import tokens2class, prosodic_string
from lingpy import *
from lingpy.align.sca import get_consensus
from lingpy import basictypes as bt


def lingrex_path(*comps):
    """
    Our data-path in CLICS.
    """
    return os.path.join(os.path.dirname(__file__), os.pardir, *comps)


def add_structure(
    wordlist, model="cv", segments="tokens", structure="structure", ref="cogid", gap="-"
):
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
