from lingpy import *
from lingrex.copar import *
from glob import glob
from tabulate import tabulate
import numpy as np

data = [
        ('burmish-240-8', 'crossids'),
        ('chinese-623-14', 'crossids'),
        ('polynesian-210-10', 'crossids'),
        ('japanese-200-10', 'crossid')
        ]

table = [[
    'dataset',
    'sites',
    'patterns',
    'singletons',
    'coverage',
    'irregulars',
    'purity'
    ]]
for f, c in data:
    name = f.split('-')[0]
    print(name)
    cp = CoPaR('data/'+f+'.tsv', ref=c, fuzzy=(c=='crossids'),
            transcription="segments", segments='segments',
            minrefs=2, structure="structure")
    cp.get_sites()
    cp.cluster_sites()
    cp.sites_to_pattern()
    cp.add_patterns()
    singletons = len([a for a in cp.clusters.items() if len(a[1]) == 1])
    cp.irregular_patterns()
    iregs = sum([len(a) for a in cp.ipatterns.values()])
    table += [[
        name,
        len(cp.sites),
        len(cp.clusters),
        singletons,
        '{0:.2f}'.format(
            (len(cp.sites)-singletons) / len(cp.sites)),
        iregs,
        cp.purity()
        ]]
    cp.output('tsv', filename='results/out-'+name)
print(tabulate(table, headers='firstrow', tablefmt='latex')) 
