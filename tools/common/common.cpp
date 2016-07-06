//
//  common.cpp
//  short read connector
//
//  Created by Pierre Peterlongo on 05/07/16.
//  Copyright (c) 2016 Pierre Peterlongo. All rights reserved.
//

#include "common.hpp"


static int NT2int(char nt){
	return (nt>>1)&3;
}


bool repeated_kmers(Kmer<KMER_SPAN(1)>::ModelCanonical& model, Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator& itKmer){
    std::unordered_set<string> kmer_seen;
    for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
        string string_kmer = model.toString(itKmer->value());
        if (kmer_seen.find(string_kmer) != kmer_seen.end()) {
            return true;
        }
        kmer_seen.insert(string_kmer);
//        std::cout << "kmer " << model.toString(itKmer->value()) << ",  value " << itKmer->value() << std::endl;
    }
    return false;
}


bool valid_sequence(Sequence& seq, const int kmer_size){
	size_t lenseq =seq.getDataSize();
    if (kmer_size>lenseq){return false;} //BUG HERE WE SHOULD NOT NEED THIS LINE (THE KMER ITERATOR CREATES FALSE KMERS WHEN SIZE OF THE SEQUENCE IS LOWER THAN K)
    
	const char* data = seq.getDataBuffer();
	int DUSTSCORE[64]={0}; // all tri-nucleotides
    
	if (data[0]!='A' && data[0]!='C' && data[0]!='G' && data[0]!='T')  { return false; }
	if (data[1]!='A' && data[1]!='C' && data[1]!='G' && data[1]!='T')  { return false; }
    
	for (int j=2; j<lenseq; ++j){
		++DUSTSCORE[NT2int(data[j-2])*16 + NT2int(data[j-1])*4 + NT2int(data[j])];
		if (data[j]!='A' && data[j]!='C' && data[j]!='G' && data[j]!='T')  { return false; }
	}
	int m,s=0;
    
	for (int i=0; i<64; ++i)
	{
		m = DUSTSCORE[i];
		s  += (m*(m-1))/2;
	}
    
	return s<((lenseq-2)/4 * (lenseq-6)/4)/2;
}

