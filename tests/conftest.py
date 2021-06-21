import pathlib

import pytest


@pytest.fixture
def data():
    return pathlib.Path(__file__).parent / 'data'
