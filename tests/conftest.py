import pathlib

import pytest


@pytest.fixture
def data():
    return pathlib.Path(__file__).parent / 'data'


@pytest.fixture
def wl_with_alignments():
    return {
        0: ["doculect", "concept", "form", "tokens", "alignment", "cogid"],
        1: ["A", "one", "atawu", "a t a w u", "a t a w u", 1],
        2: ["B", "one", "a_twu", "a _ t w u", "a t - w u", 1],
        3: ["C", "one", "tawu", "t a w u", "- t a w u", 1],
        4: ["D", "one", "tefu", "tʲ e f u", "- t e f u", 1],
        5: ["A", "two", "satu", "s a t u", "s a t u", 2],
        6: ["A", "two", "seram", "s e r a m", "s e r a m", 2]
    }
