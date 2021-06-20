from lingpy import *
from sys import argv
from lingrex.copar import CoPaR
from sys import argv

cp = CoPaR('bodt-khobwa-cleaned.tsv', ref='crossids', fuzzy=True, 
       minrefs=2, structure='structure', transcription="tokens")

# make function to extract correspondence patterns
cp.get_sites()
cp.cluster_sites()
cp.sites_to_pattern()

preds, purity, pudity = cp.predict_words()
goods = 0
with open('predictions-automatic.tsv', 'w') as f:
    f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(
        'NUMBER', 'GOOD_PREDICTION', 'COGNATESET', 'LANGUAGE', 'CONCEPT', 'MORPHEME', 'WORD1',
        'WORD2', 'WORD3'
        ))
    num = 1
    for key, vals in sorted(preds.items(), key=lambda x: x[0]):
        # get the morphemes
        idx = cp.msa['crossids'][key]['ID'][0]
        cidx = cp[idx, 'crossids'].index(key)
        try:
            morph = cp[idx, 'morphemes'][cidx]
        except:
            morph = '?'
        for doc in vals:
            val1 = ' '.join([x.split('|')[0] for x in vals[doc]])
            if "Ø" in val1:
                no = '?'
            else:
                no = ''
                goods += 1
            val2 = ' '.join(['|'.join(x.split('|')[0:2]) for x in vals[doc]])
            val3 = ' '.join(vals[doc])

            f.write('\t'.join([str(num), no, str(key), doc, cp[idx, 'concept'],
                morph, val1, val2, val3])+'\n')
            num += 1
print('useful predictions', goods)
