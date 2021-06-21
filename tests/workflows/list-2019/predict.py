import collections

from lingpy import *
from lingrex.copar import *
from sys import argv
import random
from lingpy.basictypes import *
from tabulate import tabulate
from lingpy.compare.sanity import average_coverage
import codecs


class CustomCoPaR(CoPaR):
    def stats(self, score_mode='pairs'):
        rest = 0
        if self._mode == 'fuzzy':
            for idx, cogids, alms in self.iter_rows(self._ref, self._alignment):
                for alm, cogid in zip(alms, cogids):
                    if cogid not in self.msa[self._ref]:
                        rest += len(alm)
                    else:
                        pass
        else:
            for idx, cogid, alm in self.iter_rows(self._ref, self._alignment):
                if cogid not in self.msa[self._ref]:
                    rest += len(alm)
        scores = [0 for i in range(rest)]
        for (p, ptn), sites in self.clusters.items():
            scores += len(sites) * [score_patterns(
                [
                    self.sites[site][1] for site in sites
                ], mode=score_mode)]
        return sum(scores) / len(scores)


def run_experiments(
        f, 
        ref, 
        ratio, 
        subset=None, 
        runs=100, 
        verbose=False, 
        fuzzy=True, 
        samples=1, 
        noout=False,
        score_mode='pairs'
        ):
    
    if not noout:
        outfile = codecs.open(
                'results/'+f.split('/')[-1][:-4]+'-'+str(int(ratio*100+0.5))+'.txt', 
                'w', 'utf-8')
        outfile.write('\t'.join([
                'accuracy', 'proportion', 'density', 'fuzziness', 'coverage',
                'purity', 'sounds', 'missing', 'csetsize', 'clusters', 'props',
                'patterns', 'predicted', 'predictable', 'removed', 'regular',
                'purityx'])+'\n')

    cpb = CustomCoPaR(f, ref=ref, fuzzy=fuzzy, split_on_tones=False,
            segments='segments', minrefs=2, structure="structure",
            transcription="segments")
    
    if not noout:
        inout = codecs.open(
                'results/'+f.split('/')[-1][:-4]+'-individual-'+str(int(ratio*100+0.5))+'.tsv', 
                'w', 'utf-8')
        inout.write('\t'.join(['run', 'doculect','accuracy', 'purity', 'words', 'sounds'])+'\n')
    
    # define the scores
    all_scores = []
    all_samples = set()
    all_pscores = {d: [] for d in cpb.cols}
    all_pud = {d: [] for d in cpb.cols}
    all_words = {d: [] for d in cpb.cols}
    all_sounds = {d: [] for d in cpb.cols}
    for key, msa in cpb.msa[ref].items():
        for alm, t in zip(msa['alignment'], msa['taxa']):
            all_samples.add((key, ' '.join(alm), t))
   
    for run in range(runs):    
        remove_idxs = random.sample(list(all_samples), int(len(all_samples)*ratio+0.5))
        D = {0: cpb.columns}
        for idx, cogid, alm, tax, tokens, structures in cpb.iter_rows(
                ref, 'alignment', 'doculect', 'segments', 'structure'):
            if fuzzy:
                cogids, alms, toks, strucs = [], [], [], []
                for c, a, t, s in zip(cogid, lists(alm).n, lists(tokens).n,
                        lists(structures).n):
                    if (c, str(a), tax) not in remove_idxs:
                        cogids += [c]
                        alms += [str(a)]
                        toks += [str(t)]
                        strucs += [str(s)]
                if not cogids:
                    pass
                else:
                    D[idx] = cpb[idx]
                    D[idx][cpb.header[ref]] = ints(cogids)
                    D[idx][cpb.header['segments']] = ' + '.join(toks)
                    D[idx][cpb.header['structure']] = ' + '.join(strucs)
                    D[idx][cpb.header['alignment']] = ' + '.join(alms)
            else:
                if (cogid, str(alm), tax) in remove_idxs:
                    pass
                else:
                    D[idx] = cpb[idx]
        
        cp = CustomCoPaR(D, ref=ref, fuzzy=fuzzy, split_on_tones=False,
                segments='segments', minrefs=2, structure="structure",
                transcription="segments")
        if 'l' in argv: 
            cp.load_patterns()
        else:
            cp.get_sites()
            cp.cluster_sites(score_mode=score_mode)
            cp.sites_to_pattern()

        # compute size of alphabets
        sounds = {d: collections.defaultdict(int) for d in cp.cols}
        for idx, doc, tks in cp.iter_rows('doculect', 'segments'):
            for t in tks:
                if t != '+':
                    sounds[doc][t.split('/')[1] if '/' in t else t] += 1
        ave = sum([len(s) for s in sounds.values()]) / cp.width

        # good words
        our_sample = {}
        for cogid, alm, doc in remove_idxs:
            our_sample[cogid, doc] = strings(alm)
        pscores = {d: [] for d in cp.cols}
        
        regs = sum([len(a[1]) for a in cp.clusters.items() if len(a[1]) > 1]) / len(cp.sites)
                
        predicted, purity, pudity = cp.predict_words(minrefs=2, samples=samples)
        scores = []
        unknown, all_segs, predictable, cogsize = 0, 0, 0, 0
        for k, v in predicted.items():
            for doc in v:
                if (k, doc) in our_sample and (doc == subset or not subset):
                    predictable += 1
                    cogsize += len(cp.msa[ref][k]['ID'])

                    # check for different alignments
                    msaA = cp.msa[ref][k]
                    msaB = cpb.msa[ref][k]
                    if len(msaA['alignment'][0]) != len(msaB['alignment'][0]):
                        # carve out the taxa which are still existent to find which
                        # column to delete
                        new_alm = [msaB['alignment'][i] for i in
                            range(len(msaB['alignment'])) if msaB['taxa'][i] in \
                                    msaA['taxa']]
                        almA, almB = [], []
                        for i in range(len(msaA['alignment'][0])):
                            almA += [tuple([line[i] for line in msaA['alignment']])]
                        for i in range(len(msaB['alignment'][0])):
                            almB += [tuple([line[i] for line in new_alm])]
                        out = []
                        for i, col in enumerate(almB):
                            if col not in almA:
                                out += [i]
                    else:
                        out = []

                    wA, wB = v[doc], our_sample[k, doc]
                    ms = 0
                    wB = strings([x for i, x in enumerate(wB) if i not in out]) 
                    for a, b in zip(wA, wB):
                        b = b.split('/')[1] if '/' in b else b
                        a = a.split('|')
                        for i, a_ in enumerate(a):
                            if b == a_:
                                ms += 1 * (1/(i+1))
                        if a[0] == 'Ø':
                            unknown += 1
                        all_segs += 1
        
                    score = ms / len(wA)
                    pscores[doc] += [score]
                    if verbose: 
                        print('{0:5}\t{1:15}\t{2:20}\t{3:20}\t{4:.2f}\t{5}'.format(
                            str(k), doc, str(wA), str(wB), score, len(set(msaA['taxa']))))
                    if verbose and score != 1.0:
                        purs = []
                        for i, elm in enumerate(wA):
                            if (k, i) in purity:
                                purs += ['{0:.2f}'.format(purity[k, i][doc])]
                            else:
                                purs += ['?']
                                print((cogid, i) in cp.sites)
                                print([_s for _s in cp.sites if _s[0] == cogid],
                                        cogid)
                        print('<---')
                        print('\t'.join([x for x in wA]))
                        print('\t'.join([x for x in wB]))
                        print('\t'.join(purs))
                        print('--->')
                    scores += [score]
        ubound = cp.upper_bound()
        all_scores += [(
            sum(scores) / len(scores),         
            len(cp) / len(cpb),
            density(cp, ref=ref), 
            cp.fuzziness(),
            cp.stats(score_mode=score_mode), 
            sum(pudity.values()) / len(pudity.values()), 
            ave, 
            unknown/all_segs,
            cogsize / predictable, 
            len(cp.clusters), 
            len(cp.clusters) / ubound,
            len(cp.sites),
            predictable / len(remove_idxs),
            predictable,
            len(remove_idxs),
            regs,
            cp.purity()
            )]
        if verbose:
            print('{0:.2f}'.format(all_scores[-1][0]))
        
        cov = cp.coverage()
        for p in pscores:
            all_pscores[p] += [sum(pscores[p]) / len(pscores[p])]
            all_pud[p] += [pudity[p]]
            all_words[p] += [cov[p]]
            all_sounds[p] += [len(sounds[p])]
            
            if not noout:
                inout.write('\t'.join([
                    str(run+1),
                    p, 
                    str(all_pscores[p][-1]), 
                    str(pudity[p]), 
                    str(cov[p]),
                    str(len(sounds[p]))
                    ])+'\n')
        if not noout:
            outfile.write(str(run+1)+'\t'+'\t'.join(['{0:.4f}'.format(x) for x in
                all_scores[-1]])+'\n')
        print('{0:.2f} / {1:.2f}'.format(sum(scores) / len(scores), len(cp) /
            len(cpb)))

    
    new_scores = [[
            'accuracy', 'proportion', 'density', 'fuzziness', 'coverage',
            'purity', 'sounds', 'missing', 'csetsize', 'clusters', 'props',
            'patterns', 'predicted', 'predictable', 'removed', 'regs', 'purityx']]
    new_scores += [[
        round(sum([x[0] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[1] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[2] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[3] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[4] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[5] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[6] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[7] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[8] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[9] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[10] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[11] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[12] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[13] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[14] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[15] for x in all_scores]) / len(all_scores), 4),
        round(sum([x[16] for x in all_scores]) / len(all_scores), 4),
            ]]
    if not noout:
        outfile.close()
        inout.close()


    if noout:
        print(tabulate(new_scores, headers='firstrow'))
        
    return purity, pudity, sounds, cp

if __name__ == '__main__':
    from sys import argv

    # defaults
    f = argv[1]
    ref = 'crossids'
    ratio = 0.5
    proto = None
    verbose = False
    runs = 100
    samples = 1
    noout = False
    
    # parse arguments
    if '-r' in argv:
        ratio = float(argv[argv.index('-r')+1])
    if '-c' in argv:
        ref = argv[argv.index('-c')+1]
    if '-v' in argv or '--verbose' in argv:
        verbose = True
    if '--runs' in argv:
        runs = int(argv[argv.index('--runs')+1])
    if ref in ['crossid', 'cogid']:
        fuzzy = False
    else:
        fuzzy = True
    if '--samples' in argv:
        samples = int(argv[argv.index('--samples')+1])
    if '--noout' in argv:
        noout = True

    if '--seed' in argv:
        random.seed(1)
        
    p1, p2, p3, cop = run_experiments(
            f, 
            ref, 
            ratio, 
            fuzzy=fuzzy, 
            verbose=verbose,
            runs=runs, 
            samples=samples, 
            noout=noout,
            )
    if verbose:
        cop.add_patterns()
        cop.output(
                'tsv',
                filename='results/'+f.split('/')[1].split('-')[0]+str(int(100*ratio+0.5)),
                ignore='all',
                prettify=False
                )

    
