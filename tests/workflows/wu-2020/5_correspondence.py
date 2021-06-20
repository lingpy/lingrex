from lingrex.copar import CoPaR
from sys import argv

if 'all' in argv:
    fname='A_Chen_'
else:
    fname='D_Chen_'

cop = CoPaR(
        fname+'crossids.tsv',
        ref='crossids',
        fuzzy=True,
        segments='tokens',
        minrefs=3,
        structure="structure"
        )
cop.get_sites()
cop.cluster_sites()
cop.sites_to_pattern()
cop.add_patterns()
cop.write_patterns(fname+'all_patterns.tsv')
cop.output('tsv', filename=fname+'patterns', prettify=False)

# statistics 
sps=['i','m','n','c','t']

total_correspondence_sets = len(cop.clusters)
print('{0}: {1}'.format('The total sound correspondence cluster sets', total_correspondence_sets))

print('The number of regular correspondence sets in each position')
for sp in sps:
    t = [x[1] for x, y in cop.clusters.items() if len(y)>1 and x[0] ==sp]
    print('{0}: {1}'.format(sp, len(t)))

print('The number of singletons in each position ')
for sp in sps:
    t = [x[1] for x, y in cop.clusters.items() if len(y)==1 and x[0] ==sp]
    print('{0}: {1}'.format(sp, len(t)))
