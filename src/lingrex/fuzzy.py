"""Create fuzzy reconstructions."""
from lingrex.reconstruct import PatternReconstructor
import random
from lingpy.util import pb as progressbar 
import lingpy


def ntile(words, n=5, gap="-", missing="Ø"):
    """
    Represent aligned words in form of n-tiles.
    """
    if len(words) == 1:
        return ' '.join([x for x in words[0] if x != gap])

    # start counting the occurrences
    cols = []
    for i in range(len(words[0])):
        col = [line[i] for line in words]
        cols += [col]
    
    ntile = len(words) / n 

    sounds = []
    for col in cols:
        col = [x for x in col if x != missing]
        if not col:
            sounds += ['?']
        else:
            ntile = len(col) / n
            dist = {}
            sounds += [[]]
            for s in set(col):
                dist[s] = int(col.count(s) / ntile + 0.5)
            for s, t in sorted(dist.items(), key=lambda x: x[1], reverse=True):
                for i in range(t):
                    sounds[-1] += [s]
            iterated = 0
            while len(sounds[-1]) < n:
                sounds[-1] += sounds[-1]
                iterated += 1
                if iterated >= n:
                    sounds[-1] += n * ["Ø"]
            sounds[-1] = sorted(sounds[-1][:n], key=lambda x:
                    sounds[-1].count(x), reverse=True)
            sounds[-1] = '|'.join(sounds[-1])

    return ' '.join([s for s in sounds if s.split('|').count(gap) !=
        len(s.split('|'))-1])


class FuzzyReconstructor:

    def __init__(
            self, infile, target, ref="cogid", fuzzy=False,
            transcription="form"):
        
        if isinstance(infile, (str, dict)):
            wordlist = lingpy.align.sca.Alignments(infile, ref=ref,
                    transcription=transcription)
        elif isinstance(
                infile, 
                (
                    lingpy.align.sca.Alignments, 
                    lingpy.basic.wordlist.Wordlist
                    )):
            wordlist = infile
        else:
            raise ValueError("Argument for infile must be a string or a wordlist.")
        self.wordlist = wordlist
        self.target = target or self.wordlist.cols[0]
        self.ref = ref
        self.fuzzy = fuzzy

    def random_splits(self, splits=10, retain=0.9):

        idxs = [
                idx for idx in self.wordlist if self.wordlist[idx, "doculect"] \
                        != self.target]
        tidxs = self.wordlist.get_list(col=self.target, flat=True)
        cogids = [self.wordlist[idx, self.ref] for idx in tidxs]

        self.samples = []
        for i in range(splits):
            self.samples += [random.sample(idxs, int(retain*len(idxs)+0.5))]

        self.wordlists = {}
        for i, sample in enumerate(self.samples):
            D = {0: [c for c in self.wordlist.columns]}
            for idx in sample:
                D[idx] = [self.wordlist[idx, c] for c in D[0]]
            selected_cogids = [self.wordlist[idx, self.ref] for idx in sample]
            for cogid, tidx in zip(cogids, tidxs):
                if cogid in selected_cogids:
                    D[tidx] = [self.wordlist[tidx, c] for c in D[0]]
            self.wordlists[i] = PatternReconstructor(
                    D,
                    ref=self.ref,
                    target=self.target,
                    fuzzy=self.fuzzy
                    )


    def fit_samples(self, clf, onehot=False, func=None, aligned=False, pb=False):
        pb = progressbar if pb else lambda x, desc: x
        for i, wordlist in pb(self.wordlists.items(), desc="fitting data"):
            wordlist.fit(clf=clf(), onehot=onehot, func=func, aligned=aligned)


    def predict(
            self, alignment, languages, desegment=True, orchar="¦",
            scorechar=":", output="percentiles"):
        words = []
        for i, wordlist in self.wordlists.items():
            word = wordlist.predict(alignment, languages,
                    desegment=False)
            words += [word]
        # transform to dictionary
        counts = {i: [] for i in range(len(words[0]))}
        for word in words:
            for i, sound in enumerate(word):
                counts[i] += [sound]
        # get percentiles
        if output in ["percentiles", "wp"]:
            out = []
            for i, sounds in sorted(counts.items(), key=lambda x: x[0]):
                distinct = {s: sounds.count(s)/len(sounds) for s in set(sounds)}
                distinct_s = ["{0}{1}{2}".format(k, scorechar, int(100*v+0.5)) for k, v in
                    sorted(distinct.items(), key=lambda x: x[1], reverse=True)]
                out += [orchar.join(distinct_s)]
            if output == "percentiles":
                return out
            return words, out
        elif output == "words":
            return words



