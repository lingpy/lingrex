from lingpy import *
import random
from itertools import combinations
from collections import defaultdict
from lingpy.sequence.generate import MCPhon
from lingpy.sequence.sound_classes import token2class

def bad_alignments(wordlist, alignment='alignment', ref='cogid'):
    """Create bad alignments for a given wordlist."""
    A = {}
    etd = wordlist.get_etymdict(ref=ref)
    for key, values in etd.items():
        idxs = []
        for v in values:
            if v:
                idxs += v
        reflexes = [wordlist[idx, 'tokens'] for idx in wordlist]


def seed_deviations(wordlist, events=50, alignment='alignment', ref='cogid',
        model='sca', structure='structure', pos='c', segments='tokens'):
    """Seed problematic sounds to the data that should not occur there."""
    # select all alignments
    alms = []
    for key, msa in wordlist.msa[ref].items():
        idxs = msa['ID']
        if len(idxs) > 3:
            alms += idxs
    selected = random.sample(alms, events)
    # prepare sound derivation
    changer = defaultdict(list)
    for key, value in Model(model).converter.items():
        changer[value] += [key]
    D = {idx: '' for idx in wordlist}
    for idx in selected:
        # get active indices
        tks = [h for h in wordlist[idx, segments]]
        struc = wordlist[idx, structure].split()
        active = [i for i in range(len(tks)) if struc[i] == pos]
        chosen = random.choice(active)
        sc = token2class(tks[chosen], model)
        new_sound = random.choice([x for x in changer[sc] if x != tks[chosen]]
                or ['?'])
        tks[chosen] = new_sound
        alm = class2tokens(tks, wordlist[idx, alignment])
        D[idx] = '{0}//{1}//{2}'.format(alm.index(new_sound), tks[chosen], new_sound)

        wordlist[idx, alignment] = ' '.join(alm)
    for idx, alm in wordlist.iter_rows(alignment):
        wordlist[idx, alignment] = alm.split()
    wordlist.add_entries('original_form', D, lambda x: x)

def seed_borrowings(wordlist, events=50, pairs=10, segments='tokens',
        ref='cogid'):
    """Seed fake borrowings in a dataset.

    Parameters
    ----------
    wordlist: ~lingpy.basic.wordlist.Wordlist
        Wordlist object, which should have values for segments and cognate
        identifiers.
    events: int
        The presumed number of borrowing events.
    pairs: int
        The number of language pairs among which this should be tested.
    """

    # get the language pairs
    doc_ids = list(range(wordlist.width))
    doc_pairs = list(combinations(doc_ids, r=2)) + [(y, x) for x, y in
            combinations(doc_ids, r=2)]
    if len(doc_pairs) < pairs:
        raise ValueError(
                'Your pairs must be less than there are possible' 
                ' doculect combinations')

    pair_sample = random.sample(doc_pairs, pairs)

    # get cognate sets for all pairs and just append the indices (from -> to)
    base_list = []
    for idxA, idxB in pair_sample:
        docA, docB = (
                wordlist.get_list(col=wordlist.cols[idxA]),
                wordlist.get_list(col=wordlist.cols[idxB]))
        for i, j in zip(docA, docB):
            if i != 0:
                base_list += [(i, j)]
    # now sample a number of transfers, and insert them in a new wordlist
    if events > len(base_list):
        raise ValueError('The number of events is larger than the number' 
                ' of choices')
    # list of potential transfers
    transfers = random.sample(base_list, events)
    # list of transfers including multiple ones (where we'll shuffle)
    transferred = defaultdict(list)
    for i, j in transfers:
        transferred[j] += [i]
    # random choice one from the list only which word is given the preference
    for i, j in transferred.items():
        if i != 0:
            transferred[i] = random.choice(j)
    
    # prepare cloning
    D = {0: [c for c in wordlist.columns] + ['seeds', 'original_form']}
    for idx, tks in wordlist.iter_rows(segments):
        row = [h for h in wordlist[idx]]
        if idx in transferred:
            # get the source
            source = wordlist[transferred[idx], segments]
            row += [transferred[idx]]
            row += [' '.join(tks)]
            row[wordlist.header[segments]] = source
            row[wordlist.header[ref]] = wordlist[transferred[idx], ref]
            row[-2] = transferred[idx]
        else:
            row += ['', '']
        D[idx] = row

    # return stuff
    return Wordlist(D)

def seed_neologs(wordlist, events=50, doculects=2, segments='tokens',
        ref='cogid', no_singletons=False):
    # get the wordlists
    selected = random.sample(wordlist.cols, doculects)
    
    # check for singletons
    singletons = []
    if no_singletons:
        for cogid, refs in wordlist.get_etymdict(ref=ref).items():
            if len([x for x in refs if x]) == 1:
                singletons += [x for x in refs if x][0]

    # make markov chains
    chains = {}
    all_indices = []
    for doc in selected:
        idxs = wordlist.get_list(col=doc, flat=True)
        tks = [wordlist[idx, segments] for idx in idxs]
        ipa = [' '.join(tk) for tk in tks]
        chains[doc] = MCPhon(ipa)
        all_indices += [idx for idx in idxs if idx not in singletons]
        
    sel_idxs = random.sample(all_indices, events) if len(all_indices) >= events \
            else all_indices
    new_words = {}
    for idx in sel_idxs:
        new_word = chains[wordlist[idx, 'doculect']].get_string(new=True,
                tokens=True)
        new_words[idx] = [x[1] for x in new_word]
    # prepare cloning
    D = {0: [c for c in wordlist.columns] + ['original_form']}
    for idx, tks in wordlist.iter_rows(segments):
        row = [h for h in wordlist[idx]]
        if idx in new_words:
            # get the source
            row += [' '.join(tks)]
            row[wordlist.header[segments]] = new_words[idx]
        else:
            row += ['']
        D[idx] = row
    return Wordlist(D)


        
