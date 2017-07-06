from pyburmish import burmish_path
from lingrex.colex import find_colexified_alignments
from lingpy import *

alms = Alignments(burmish_path('dumps', 'alignments.tsv'), ref='cogids',
        segments='tokens')
merged = find_colexified_alignments(alms)
alms.output('tsv', filename='alignments-merged.tsv')
