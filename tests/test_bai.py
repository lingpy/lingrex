from pycddb.dataset import Dataset
from lingpy import *
from lingrex.cpa import *
from lingrex.util import lingrex_path, align_by_structure
from lingrex.colex import *

print('start')
alms = Alignments(Dataset('Allen2007').words, ref='cogids',
        segments='segments')
alms.add_entries('tokens', 'segments', lambda x: x)
align_by_structure(alms, template='imnNcCt', segments='tokens', ref='cogids',
        structure='structure', alignment='alignment')

# add alignments
find_colexified_alignments(alms, cognates='cogids', ref='crossids')
alms.add_alignments(ref='crossids', fuzzy=True)

print('added alignments')
# open the file to which we write the patterns which we find
f = open('bai-consensus.tsv', 'w')
f.write('ID\tCONTEXT\t'+'\t'.join(alms.taxa)+'\tFREQUENCY\tINFERRED\tPROTO\tCOGNATES\n')

print('opened consensus')
consensus = {}
for c in 'imnNcCt':
    graph = compatibility_graph(alms, ref='crossids', prostring='structure',
            pos=c,
            minimal_taxa = 3)
    graph = greedy_clique_partition(graph, debug=False)
    consensus[c] = get_consensus_patterns(graph)
    for i, (a, b, cc) in consensus[c].items():
        f.write('\t'.join([str(i),
            c]+a+[str(b)]+[reconstruct_from_consensus(a,
                graph=False,
                rule='majority')]+['']+ [cc])+'\n')
print('adding proto-forms')
# add proto-slots for all words that are sufficiently reflected
wl = add_proto_slots(alms, 'Proto-Bai', segments='tokens',
        consensus=consensus, rule='majority',
        template='imnNcCt', ref='crossids',
        graph=False)
wl.output('tsv', filename='proto-bai')
