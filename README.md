# LingRex: Linguistic Reconstruction with LingPy

[![Build Status](https://github.com/lingpy/lingrex/workflows/tests/badge.svg)](https://github.com/lingpy/lingrex/actions?query=workflow%3Atests)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.1544943.svg)](https://doi.org/10.5281/zenodo.1544943)
[![PyPI version](https://badge.fury.io/py/lingrex.png)](https://badge.fury.io/py/lingrex)

LingRex offers the code needed for the automatic inference of sound correspondence patterns as described in the following paper:

> List, J.-M. (2019): Automatic inference of sound correspondence patterns across multiple languages. Computational Linguistics 45.1. 137-161. [DOI: 10.1162/coli_a_00344](https://doi.org/10.1162/coli_a_00344)

To test this workflow, please check the workflow code example in [`tests/workflows/list-2019`](tests/workflows/list-2019).

LingRex offers also the code needed for a baseline algorithm for automatic word prediction or automatic phonological reconstruction in a supervised fashion.

> List, J.-M. and R. Forkel and N. W. Hill (2022): A New Framework for Fast Automated Phonological Reconstruction Using Trimmed Alignments and Sound Correspondence Patterns. Proceedings of the 3rd International Workshop on Computational Approaches to Historical Language Change (LChange 2022). Dublin. Ireland. https://aclanthology.org/2022.lchange-1.9

This algorithm is also used as a baseline for a Shared Task on the Prediction of Cognate Reflexes (https://sigtyp.github.io/st2022.html), organized as part of the SIGTYP Workshop at NAACL 2022.

> List, J.-M., E. Vylomova, R. Forkel, N. Hill, and R. Cotterell (2022): The SIGTYP shared task on the prediction of cognate reflexes. In: Proceedings of the 4th Workshop on Computational Typology and Multilingual NLP. Association for Computational Linguistics 52-62. https://aclanthology.org/2022.sigtyp-1.7

Methods for the handling of partial cognates were introduced in a study by Wu and List (2023):

> Wu, M.-S. and J.-M. List (2023): Annotating cognates in phylogenetic studies of South-East Asian languages. Language Dynamics and Change. https://doi.org/10.1163/22105832-bja10023

Methods for the trimming of phonetic alignments were introduced in a study by Blum and List (2023):

> Blum, F. and J.-M. List (2023): Trimming phonetic alignments improves the inference of sound correspondence patterns from multilingual wordlists. In: Proceedings of the 5th Workshop on Computational Typology and Multilingual NLP. Association for Computational Linguistics 52-64. https://aclanthology.org/2023.sigtyp-1.6.pdf

Methods for the handling and creation of fuzzy / uncertain phonological reconstructions were introduced in a study by List et al. (forthcoming):

> List, J.-M.; Hill, N. W.; Blum, F.; and Forkel, R. (forthcoming): A New Framework for the Representation and Computation of Uncertainty in Phonological Reconstruction. To appear in:  Proceedings of the 4th Workshop on Computational Approaches to Historical Language Change.

When using this package in your research, please make sure to quote the respective papers, depending on the algorithms you use, and quote the software package as follows:

> List, J.-M. and R. Forkel (2023): LingRex: Linguistic Reconstruction with LingPy. [Computer software, Version 1.4.0]. With contributions by Frederic Blum and Mei-Shin Wu. Leipzig: Max Planck Institute for Evolutionary Anthropology. https://pypi.org/project/lingrex

Since this software package itself makes use of LingPy's alignment algorithms, you should also quote the LingPy package itself.

> List, J.-M. and R. Forkel (2023): LingPy. A Python library for quantitative tasks in historical linguistics. Version 2.6.10. Max Planck Institute for Evolutionary Anthropology: Leipzig. https://lingpy.org

## Installation

Install the package via `pip`:

```shell
pip install lingrex
```

## Further Examples

The borrowing detection algorithm implemented in LingRex is introduced in the
paper:

> List, J.-M. and R. Forkel (2021): Automated identification of borrowings in multilingual wordlists [version 1; peer review: 3 approved, 1 approved with reservations]. Open Research Europe 1.79. 1-11. [DOI: 10.12688/openreseurope.13843.1](https://doi.org/10.12688/openreseurope.13843.1)

If you use this algorithm, please cite LingRex and this paper.

In addition to the paper in which the correspondence pattern inference algorithm was first introduced, LingRex also offers the code to compute the workflow described in the following paper:

> Wu, M.-S., N. Schweikhard, T. Bodt, N. Hill, and J.-M. List (2020): Computer-Assisted Language Comparison. State of the Art. Journal of Open Humanities Data 6.2. 1-14. [DOI: 10.5334/johd.12](https://doi.org/10.5334/johd.12)

To test this workflow, please check the workflow code example in `tests/workflows/wu-2020`. 

If you use this workflow in your work, please quote this paper as well.

In addition, our experiment (with T. Bodt) on predicting words with the help of sound correspondence patterns also made use of the LingRex package.

> Bodt, T. and J.-M. List (2021): Reflex prediction. A case study of Western Kho-Bwa. Diachronica 0.0. 1-38. [DOI: 10.1075/dia.20009.bod](https://doi.org/10.1075/dia.20009.bod)

To test this workflow, please check the workflow code example in `tests/workflows/bodt-2019`.
