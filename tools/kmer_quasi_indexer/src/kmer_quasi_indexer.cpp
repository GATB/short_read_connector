//! [snippet1]

#include <kmer_quasi_indexer.hpp>

using namespace std;

/********************************************************************************/

// We define some constant strings for names of command line parameters
static const char* STR_FOO = "-foo";

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
    getParser()->push_front (new OptionOneParam (STR_FOO, "my option",  false, "1"));
    getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
    getParser()->push_back (new OptionOneParam (STR_VERBOSE,   "verbosity (0:no display, 1: display kmers, 2: display distrib",  false, "0"));

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
	// We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");


	int kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());

	double lg2 = log(2);
	float NBITS_PER_KMER = log (16*kmer_size*(lg2*lg2))/(lg2*lg2);
	NBITS_PER_KMER = 12;
	/** We get the number of solid kmers. */
	const u_int64_t solidFileSize = solidKmers.getNbItems();
	u_int64_t estimatedBloomSize = (u_int64_t) ((double)solidFileSize * NBITS_PER_KMER);
	if (estimatedBloomSize ==0 ) { estimatedBloomSize = 1000; }
	int nb_hash				= 7;
	int fingerprint_size 	= 8;


	nbSolidKmers = solidKmers.getNbItems();

	std::vector< list<u_int32_t> > all_lists(nbSolidKmers);


	IteratorKmerH5Wrapper iteratorOnKmers(solidKmers.iterator(), kmer_size);


	//BUG ICI:
	int i=0;
	for(auto & element : iteratorOnKmers){
		cout<<element<<endl;
		i++;
		if (i>10) break;
	}

    cout << "---------------------------------------------------- " << endl;

	i=0;
	for(auto & element : iteratorOnKmers){
		cout<<element<<endl;
		i++;
		if (i>10) break;
	}

    cout << "---------------------------------------------------- " << endl;

//	quasiDico = quasiDictionnary<IteratorKmerH5, std::vector< list<u_int32_t> > (nbSolidKmers, iteratorOnKmers, all_list.begin(), fingerprint_size, sizeof(list<u_int32_t>));
	quasiDico = quasiDictionnary<IteratorKmerH5Wrapper, std::vector< list<u_int32_t> >::iterator > (nbSolidKmers, iteratorOnKmers, fingerprint_size, sizeof(list<u_int32_t>));

	quasiDico.createGenericValues<list<u_int32_t> >();
	          //quasiDictionnary                                                (u_int64_t nelement, RangeKeyOnly& itKey, RangeKeyValue& it, const int fingerprint_size, const int value_size, double gammaFactor=1, int nthreads=1)


	// Ici ou dans une autre fonction: parcourir les reads. Pour les kmers présents dans le quasidictionnary: ajouter l'id du read dans leur liste de reads (valeur du quasi dictionary).

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
    getInfo()->add (2, STR_FOO,  "%d",  getInput()->getInt(STR_FOO));
    getInfo()->add (1, &LibraryInfo::getInfo());
    create_quasi_dictionary();
}

