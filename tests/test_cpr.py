from lingrex.cpr import *
from lingpy import *


cpu = CorPatR('GER.tsv', ref='cogid')
cpu.align()
D = {}
for idx, tks in cpu.iter_rows('tokens'):
    structure = []
    for t in tks:
        if t not in '_#':
            structure += ['c']
        else:
            structure += ['+']
    D[idx] = ' '.join(structure)
cpu.add_entries('structure', D, lambda x: x)
patterns, allpats = cpu.get_patterns()
patterns, allpats = cpu.get_patterns(pos='c', prostring='structure')
cpu.sort_patterns(patterns, threshold=2, debug=True)
cpu.cluster_patterns(debug=True)
