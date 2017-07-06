from lingrex.util import lingrex_path, renumber_partials, save_network
from pyburmish import burmish_path
from lingrex.cpa import (
        compatibility_graph, greedy_clique_partition,
        get_consensus_patterns, reconstruct_from_consensus, add_proto_slots,
        parse_sound_change_graph)
from lingpy import *
from lingrex.colex import find_colexified_alignments


# first step, load the alignments in the data
alms = Alignments(burmish_path('dumps', 'alignments.tsv'), ref='cogids',
        segments='tokens')
# second step, merge the alignments
merged = find_colexified_alignments(alms, cognates='cogids', ref='crossids')

alms.add_alignments(ref='crossids', fuzzy=True)

# prepare the sound-change graph
sc_graph = parse_sound_change_graph(burmish_path('reconstruction/change-directions.tsv'))

# open the file to which we write the patterns which we find
f = open(burmish_path('reconstruction', 'consensus.tsv'), 'w')
f.write('ID\tCONTEXT\t'+'\t'.join(alms.taxa)+'\tFREQUENCY\tINFERRED\tPROTO\tCOGNATES\n')

# write the sound change network to file to make it possible to look at it
save_network('sounds.gml', sc_graph)


print(alms.header)
# get the consensus strings for the different patterns
consensus = {}
for c in 'imMnNct':
    graph = compatibility_graph(alms, ref='crossids', prostring='structure',
            pos=c,
            minimal_taxa = 3)
    graph = greedy_clique_partition(graph, debug=False)
    consensus[c] = get_consensus_patterns(graph)
    for i, (a, b, cc) in consensus[c].items():
        f.write('\t'.join([str(i),
            c]+a+[str(b)]+[reconstruct_from_consensus(a,
                graph=sc_graph,
                rule='network')]+['']+ [cc])+'\n')
input()
# add proto-slots for all words that are sufficiently reflected
wl = add_proto_slots(alms, 'Proto-Burmish', segments='tokens',
        consensus=consensus, rule='network',
        graph=sc_graph)

# write data to file
wl.output('tsv', filename=burmish_path('dumps', 'proto-burmish'))


