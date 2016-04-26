//! [snippet1]

#include <kmer_quasi_indexer.hpp>

using namespace std;

/********************************************************************************/

// We define some constant strings for names of command line parameters
//static const char* STR_FOO = "-foo";

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
kmer_quasi_indexer::kmer_quasi_indexer ()  : Tool ("kmer_quasi_indexer")
{
    // We add some custom arguments for command line interface
    getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
//    getParser()->push_back (new OptionOneParam (STR_VERBOSE,   "verbosity (0:no display, 1: display kmers, 2: display distrib",  false, "0"));
    getParser()->push_back (new OptionOneParam (STR_URI_INPUT, "bank input",    true));


}



/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void kmer_quasi_indexer::create_quasi_dictionary (){
	const int display = getInput()->getInt (STR_VERBOSE);
	// We get a handle on the HDF5 storage object.
	// Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
	auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (getInput()->getStr(STR_URI_GRAPH)));
	// We get the group for dsk
	Group& dskGroup = storage->getGroup("dsk");
	kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());

	// We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");



	double lg2 = log(2);
	float NBITS_PER_KMER = log (16*kmer_size*(lg2*lg2))/(lg2*lg2);
	NBITS_PER_KMER = 12;
	/** We get the number of solid kmers. */
	const u_int64_t solidFileSize = solidKmers.getNbItems();
//	u_int64_t estimatedBloomSize = (u_int64_t) ((double)solidFileSize * NBITS_PER_KMER);
//	if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }
//	int nb_hash				= 7;


	nbSolidKmers = solidKmers.getNbItems();



	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());

//	int i=0;
//    cout << "---------------------------------------------------- " << endl;
//
//
//	for(auto  &element: iteratorOnKmers){
//		cout<<"first  "<<element<<endl;
//	}
//
//    cout << "---------------------------------------------------- " << endl;
//	i=0;
//	for(auto  &element : iteratorOnKmers){
//		i++;
//		if (i<10) {cout<<"next 10: "<<element<<endl;}
//		else break;
//	}
//    cout << "---------------------------------------------------- " << endl;


	int fingerprint_size 	= 8; //TODO: must be a parameter and set to zero if query == reference bank

//	quasiDico = quasiDictionnary<IteratorKmerH5, std::vector< list<u_int32_t> > (nbSolidKmers, iteratorOnKmers, all_list.begin(), fingerprint_size, sizeof(list<u_int32_t>));
	quasiDico = quasiDictionnaryKeyGeneric<IteratorKmerH5Wrapper, u_int32_t> (nbSolidKmers, iteratorOnKmers, fingerprint_size);



}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void kmer_quasi_indexer::fill_quasi_dictionary (){
	bool verbose = getInput()->get(STR_VERBOSE) != 0;
	bool exists;
	u_int32_t nbSequences = 0;
    IBank* bank = Bank::open (getInput()->getStr(STR_URI_INPUT));
    cout<<kmer_size<<"-mers from banks "<<getInput()->getStr(STR_URI_INPUT)<<endl;
    LOCAL (bank);
    // We declare a canonical model with a given span size.
    Kmer<span>::ModelCanonical model (kmer_size);
    // We declare an iterator on a given sequence.
    Kmer<span>::ModelCanonical::Iterator itKmer (model);
    // We create an iterator over this bank.
    ProgressIterator<Sequence> itSeq (*bank);
    // We loop over sequences.
    for (itSeq.first(); !itSeq.isDone(); itSeq.next())
    {
        // We set the data from which we want to extract kmers.
        itKmer.setData (itSeq->getData());
        // We iterate the kmers.
        for (itKmer.first(); !itKmer.isDone(); itKmer.next())
        {
//            if (verbose)  {  cout <<itKmer->value() << " " <<model.toString (itKmer->value()) << " " << itKmer->value().getVal()<<endl;  }
//            std::cout << "forward value  is: " << itKmer->forward()                    << std::endl;
//            std::cout << "forward string is: " << model.toString(itKmer->forward())    << std::endl;
//            std::cout << "revcomp value  is: " << itKmer->revcomp()                    << std::endl;
//            std::cout << "revcomp string is: " << model.toString(itKmer->revcomp())    << std::endl;
//            std::cout << "used strand is   : " << toString(itKmer->strand())           << std::endl;


        	// Adding the read id to the list of ids associated to this kmer.
        	// note that the kmer may not exist in the dictionnary if it was under the solidity threshold.
        	// in this case, nothing is done
            quasiDico.set_value(itKmer->value().getVal(), nbSequences);
        }
        //  We increase the sequences counter.
        nbSequences++;
    }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void kmer_quasi_indexer::execute ()
{
    // We can do here anything we want.
    // For further information about the Tool class, please have a look
    // on the ToyTool snippet  (http://gatb-core.gforge.inria.fr/snippets_tools.html)

    // We gather some statistics.
    getInfo()->add (1, "input");
//    getInfo()->add (2, STR_FOO,  "%d",  getInput()->getInt(STR_FOO));
    getInfo()->add (1, &LibraryInfo::getInfo());
    create_quasi_dictionary();
    fill_quasi_dictionary();
}

