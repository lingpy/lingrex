import shutil
import pathlib
import subprocess

import pytest


@pytest.fixture
def clean_dir(tmp_path):  # pragma: no cover
    def _clean_dir(d):
        shutil.copytree(pathlib.Path(__file__).parent / 'workflows' / d, tmp_path / d)
        return tmp_path / d
    return _clean_dir


def _run(wd, *cmds):  # pragma: no cover
    for cmd in cmds:
        try:
            subprocess.check_call(cmd, cwd=wd, shell=True)
        except subprocess.CalledProcessError as e:  # pragma: no cover
            print(e)
            print(e.output)
            raise


@pytest.mark.workflow
def test_bodt(clean_dir):  # pragma: no cover
    _run(
        clean_dir('bodt-2019'),
        'python predict.py',
        'python test-prediction.py bodt-khobwa-cleaned.tsv -r 0.5',
    )


@pytest.mark.workflow
def test_list(clean_dir):  # pragma: no cover
    _run(
        clean_dir('list-2019'),
        'python general.py',
        'python predict.py data/burmish-240-8.tsv -r 0.75 --runs 2',
        'python predict.py data/chinese-623-14.tsv -r 0.75 --runs 2',
        'python predict.py data/polynesian-210-10.tsv -r 0.75 --runs 2',
        'python predict.py data/japanese-200-10.tsv -c crossid -r 0.75 --runs 2',
    )


@pytest.mark.workflow
def test_wu(clean_dir):  # pragma: no cover
    _run(
        clean_dir('wu-2020'),
        'python 4_crosssemantic.py',
        'python 5_correspondence.py',
    )
