from lingpy import *

def expand_alignment(msa, taxa, missing="Ø"):
    """Expand an alignment by adding a symbol for missing taxa"""
    out = []
    for taxon in taxa:
        if taxon in msa['taxa']:
            tidx = msa['taxa'].index(taxon)
            out += [msa['alignment'][tidx]]
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
    for cogid, msa in alignments.msa[cognates].items():
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
    alignments.add_entries(ref, cognates, lambda x: [merged.get(y, y) for y in x])


