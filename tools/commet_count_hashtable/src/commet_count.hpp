/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _TOOL_commet_count_HPP_
#define _TOOL_commet_count_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <gatb/system/impl/SystemInfoCommon.hpp>// for having the memory

#include "../../../thirdparty/IteratorKmerH5/IteratorKmerH5.hpp"
#include "../../../thirdparty/quasi_dictionnary/src/quasidictionnary.h"
/********************************************************************************/


////////////////////////////////////////////////////////////////////////////////
//
// THIS FILE IS AUTOMATICALLY GENERATED...
//
// THIS IS A SIMPLE EXAMPLE HOW TO USE THE Tool CLASS. IF YOU WANT MORE FEATURES,
// YOU CAN HAVE A LOOK AT THE ToyTool SNIPPET HERE:
//
//      http://gatb-core.gforge.inria.fr/snippets_tools.html
//
////////////////////////////////////////////////////////////////////////////////

class commet_count : public Tool
{
	
private:
//	quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char > quasiDico;
	std::unordered_map<u_int64_t, unsigned char> hashDico;
	u_int64_t nbSolidKmers;
	int kmer_size;
	static const size_t span = KMER_SPAN(1);

public:

    // Constructor
	commet_count ();

    // Actual job done by the tool is here
    void execute ();

    void create_and_fill_quasi_dictionary (int fingerprint_size, const int nbCores);


    void parse_query_sequences(int threshold, const int nbCores);
};

/********************************************************************************/

#endif /* _TOOL_commet_count_HPP_ */
