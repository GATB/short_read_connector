//
//  common.cpp
//  short read connector
//
//  Created by Pierre Peterlongo on 05/07/16.
//  Copyright (c) 2016 Pierre Peterlongo. All rights reserved.
//

#include "common.hpp"

using namespace std;


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

template<class T> T min(T a, T b, T c) {
    return a < b ? std::min(a, c) : std::min(b, c);
}

size_t edit_distance(const string& A, const string& B)
{
    
    
        
    
    
    // LINEAR MEMORY
    
    size_t NA = A.size();
    size_t NB = B.size();
    
    vector<size_t> current(NB+1);
    vector<size_t> previous(NB+1);
    
    for (size_t b = 0; b <= NB; ++b)
        previous[b] = b;
    
    
    for (size_t a = 1; a <= NA; ++a){
        current[0]=a;
        for (size_t b = 1; b <= NB; ++b){
            size_t x = previous[b] + 1;
            size_t y = current[b-1] + 1;
            size_t z = previous[b-1] + (A[a-1] == B[b-1] ? 0 : 1);
            current[b] = min(x,y,z);
        }
        previous.swap(current);
    }
    
//    if(previous[NB]<2) cout<<A<<endl<<B<<endl<<previous[NB]<<endl;
    return previous[NB];
    
//    
//    
//    
//    
//    
//    vector<vector<size_t>> M(NA + 1, vector<size_t>(NB + 1));
//    
//    for (size_t a = 0; a <= NA; ++a)
//        M[a][0] = a;
//    
//    for (size_t b = 0; b <= NB; ++b)
//        M[0][b] = b;
//    
//    for (size_t a = 1; a <= NA; ++a)
//        for (size_t b = 1; b <= NB; ++b)
//        {
//            size_t x = M[a-1][b] + 1;
//            size_t y = M[a][b-1] + 1;
//            size_t z = M[a-1][b-1] + (A[a-1] == B[b-1] ? 0 : 1);
//            M[a][b] = min(x,y,z);
//        }
//    
//
//    return M[A.size()][B.size()];
}