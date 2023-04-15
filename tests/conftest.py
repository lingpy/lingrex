import pathlib

import pytest


@pytest.fixture
def data():
    return pathlib.Path(__file__).parent / 'data'


@pytest.fixture
def wl_with_alignments():
    return {
        0: ["doculect", "concept", "form", "tokens", "alignment", "cogid"],
        1: ["A", "one", "atawu", "ata+wu", "a t a w u", 1],
        2: ["B", "one", "a_twu", "a_twu", "a t - w u", 1],
        3: ["C", "one", "tawu", "tawu", "- t a w u", 1],
        4: ["D", "one", "tefu", "tefu", "- t e f u", 1],
        5: ["A", "two", "satu", "satu", "s a t u", 2],
        6: ["A", "two", "seram", "seram", "s e r a m", 2]
    }
