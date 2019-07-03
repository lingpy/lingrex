from lingpy import *
from collections import defaultdict

def find_bad_internal_alignments(alignments, ref='cogids', segments='tokens',
        alignment='alignment', transcription='ipa'):
    newIDs = {}
    def get_all_indices(lst):
        idxs = defaultdict(list)
        for i, l in enumerate(lst):
            idxs[l] += [i]
        return idxs
    
    new_cogid = max(alignments.msa[ref])
    for cogid, msa in alignments.msa[ref].items():
        idxs = [i for t, i in get_all_indices(msa['taxa']).items() if len(i) >
                1]
        for idx in idxs:
            for i in idx[1:]:
                newIDs[msa['ID'][i]] = (cogid, new_cogid)
                new_cogid += 1
    for idx, (cogid, new_cogid) in newIDs.items():
        this_idx = alignments[idx, ref].index(cogid)
        alignments[idx, ref][this_idx] = new_cogid
    return Alignments(alignments, ref=ref, fuzzy=alignments._mode,
            segments=segments, alignment=alignment, 
            transcription=transcription)


def expand_alignment(msa, taxa, missing="Ø"):
    """Expand an alignment by adding a symbol for missing taxa"""
    out = []
    for taxon in taxa:
        if taxon in msa['taxa']:
            tidx = msa['taxa'].index(taxon)
            out += [[x.split('/')[1] if '/' in x else x for x in
                msa['alignment'][tidx]]]
        else:
            out += [missing]
    return out

def compatible(msa1, msa2, missing="Ø"):
    """Compare two alignments and check whether they colexify."""
    matches = 0
    for line1, line2 in zip(msa1, msa2):
        if line1 == line2 and missing not in [line1, line2]:
            matches += 1
        else:
            if missing in [line1, line2]:
                pass
            else:
                return False
    return matches


def merge_alignments(msa1, msa2, missing="Ø"):
    out = []
    for line1, line2 in zip(msa1, msa2):
        if line1 == line2:
            out += [line1]
        elif line1 == missing:
            out += [line2]
        elif line2 == missing:
            out += [line1]
        else:
            print(line1, line2)
            raise ValueError("Alignments cannot be merged.")
    return out


def find_colexified_alignments(alignments, cognates='cogids', segments='tokens',
        missing="Ø", ref='crossids'):
    """Identify identical alignments in a dataset and label them as homophones"""
    
    queue = []
    for cogid, msa in sorted(alignments.msa[cognates].items(), key=lambda x:
            len(set(x[1]['taxa']))):
        queue += [(cogid, expand_alignment(msa, alignments.taxa,
            missing=missing))]

    merged = {}

    while queue:
        this_cogid, this_msa = queue.pop(0)
        deletes = []
        merged[this_cogid] = this_cogid
        for i, (other_cogid, other_msa) in enumerate(queue):
            if compatible(this_msa, other_msa) > 1:
                this_msa = merge_alignments(this_msa, other_msa)
                merged[other_cogid] = this_cogid
                deletes += [i]
                
        for i in deletes[::-1]:
            del queue[i]
    
    # assemble the clusters now
    if alignments._mode == 'fuzzy':
        alignments.add_entries(ref, cognates, lambda x: [merged.get(y, y) for y in x])
    else:
        alignments.add_entries(ref, cognates, lambda x: merged.get(x, x))


