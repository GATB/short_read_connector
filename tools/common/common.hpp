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
#include <cassert>
#include <string>
#include <vector>



bool repeated_kmers(Kmer<KMER_SPAN(1)>::ModelCanonical& model, Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator& itKmer);
bool valid_sequence(Sequence& seq, const int kmer_size);
size_t edit_distance(const string& A, const string& B);

#endif
