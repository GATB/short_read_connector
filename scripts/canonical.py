

    
comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def rev_comp(seq):
    return (''.join(comp[c] for c in seq))[::-1]
    

def canonical_representation (kmer):
    rc=rev_comp(kmer)
    if(rc<kmer): return rc
    return kmer
