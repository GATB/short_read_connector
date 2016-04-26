bank_set="data/humch1_00096_reads.fasta.gz"
query_set="data/humch1_00096_reads.fasta.gz"
# bank_set="data/c1.fasta.gz"
# query_set="data/c1.fasta.gz"
result_file="associated_reads.txt"
kmer_size=31
abundance_min=2
fingerprint_size=8
kmer_threshold=50
out_dsk="solid_kmers"


EDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# Count kmers using dsk
$EDIR/thirdparty/dsk/bin/macosx/dsk -file ${bank_set} -kmer-size ${kmer_size} -abundance-min ${abundance_min} -out ${out_dsk} 

$EDIR/build/tools/kmer_quasi_indexer/kmer_quasi_indexer -graph ${out_dsk}  -bank ${bank_set} -query ${query_set} -out ${result_file} -kmer_threshold ${kmer_threshold} -fingerprint_size {fingerprint_size}
