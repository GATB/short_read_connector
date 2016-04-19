query_set="data/c1.fasta.gz"
kmer_size=31
abundance_min=2
out_dsk="solid_kmers"


EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Count kmers using dsk
$EDIR/thirdparty/dsk/bin/macosx/dsk -file ${query_set} -kmer-size ${kmer_size} -abundance-min ${abundance_min} -out ${out_dsk} 

