//
//  common.h
//  short read connector
//
//  Created by Pierre Peterlongo on 05/07/16.
//  Copyright (c) 2016 Pierre Peterlongo. All rights reserved.
//

#ifndef short_read_connector_common_h
#define short_read_connector_common_h
#include <gatb/gatb_core.hpp>
#include <unordered_set>


static int NT2int(char nt){
	return (nt>>1)&3;
}

static int max_populated_window(const vector<bool> populated, const int windows_size, const int start=0, int stop=-1){
    
    //        cout<< "mpw  "<<start<<" "<<stop<<endl; //DEBUG
    if (stop==-1) stop=populated.size();
    //        cout<< "mpw2 "<<start<<" "<<stop<<endl; //DEBUG
    const int last_excluded_starting_window_position = stop-windows_size+1>=0?stop-windows_size+1:0;
    //        cout<<"leswp "<<last_excluded_starting_window_position<<endl; //DEBUG
    int max = 0;
    int number_populated=0;
    
    for(int i=start;i<windows_size && i<stop;i++){
        
        if (populated[i]) {
            number_populated++;
        }
    }
    max=number_populated;
    
    for (int pos=start+1;pos<last_excluded_starting_window_position;pos+=1){
        if (populated[pos-1]){number_populated--;}
        if (populated[pos+windows_size-1]){number_populated++;}
        if (number_populated>max) max=number_populated;
    }
//    cout<<max<<" "<<windows_size<<endl; //DEBUG
    return max;
}





static bool repeated_kmers(Kmer<KMER_SPAN(1)>::ModelCanonical& model, Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator& itKmer){
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

static bool is_high_complexity(Sequence& seq, const int kmer_size){


    
	size_t lenseq =seq.getDataSize();
    if (kmer_size>lenseq){return false;} //BUG HERE WE SHOULD NOT NEED THIS LINE (THE KMER ITERATOR CREATES FALSE KMERS WHEN SIZE OF THE SEQUENCE IS LOWER THAN K)
    
	const char* data = seq.getDataBuffer();
	int DUSTSCORE[64]={0}; // all tri-nucleotides
    
//	if (data[0]!='A' && data[0]!='C' && data[0]!='G' && data[0]!='T')  { return false; }
//	if (data[1]!='A' && data[1]!='C' && data[1]!='G' && data[1]!='T')  { return false; }
    
	for (int j=2; j<lenseq; ++j){
		++DUSTSCORE[NT2int(data[j-2])*16 + NT2int(data[j-1])*4 + NT2int(data[j])];
//		if (data[j]!='A' && data[j]!='C' && data[j]!='G' && data[j]!='T')  { return false; }
	}
	int m,s=0;
    
	for (int i=0; i<64; ++i)
	{
		m = DUSTSCORE[i];
		s  += (m*(m-1))/2;
	}
    
	return s<((lenseq-2)/4 * (lenseq-6)/4)/2;
}


#endif
