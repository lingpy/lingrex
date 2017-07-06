from lingrex.util import add_structure, align_by_structure, renumber_partials
from pyburmish import burmish_path
from lingpy import *

wl = Wordlist(burmish_path('dumps', 'burmish.tsv'))
renumber_partials(wl, 'cogids')
add_structure(wl)
count = 1
for idx, tokens, structures in iter_rows(wl, 'tokens', 'structure'):
    if len(structures.split(' ')) != len(tokens):
        print('\t'.join(tokens))
        print('\t'.join(structures.split(' ')))
        print(count)
        print(' ')
        count += 1
align_by_structure(wl, segments='tokens')
wl.output('tsv', filename=burmish_path('dumps', 'alignments'))
