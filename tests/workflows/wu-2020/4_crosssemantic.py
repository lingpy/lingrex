from lingpy import *
from lingrex.colex import find_colexified_alignments, find_bad_internal_alignments
from lingrex.align import template_alignment
from sys import argv

if 'all' in argv:
    fname='A_Chen_'
else:
    fname='D_Chen_'


alms = Alignments(fname+'aligned.tsv', ref='cogids')
print('[i] search for bad internal alignments')
find_bad_internal_alignments(alms)

print('[i] search for colexified alignments')
find_colexified_alignments(
        alms,
        cognates='cogids',
        segments='tokens',
        ref='crossids'
        )

# re-align the data
print('[i] re-align the data')
template_alignment(alms,
                   ref='crossids',
                   template='imnct+imnct+imnct+imnct+imnct+imnct',
                   structure = 'structure',
                   fuzzy=True,
                   segments='tokens')

alms.output('tsv', filename=fname+'crossids', prettify=False)
