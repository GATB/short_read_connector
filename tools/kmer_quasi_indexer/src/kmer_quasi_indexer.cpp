//! [snippet1]

#include <kmer_quasi_indexer.hpp>

using namespace std;

/********************************************************************************/

// We define some constant strings for names of command line parameters
//static const char* STR_FOO = "-foo";
static const char* STR_URI_BANK_INPUT = "-bank";
static const char* STR_URI_QUERY_INPUT = "-query";
static const char* STR_FINGERPRINT = "-fingerprint_size";
static const char* STR_THRESHOLD = "-kmer_threshold";
static const char* STR_OUT_FILE = "-out";

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
    getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
    getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT, "query input",    true));
    getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
    getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal number of shared kmers for considering 2 reads as similar",    false, "10"));
    getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "8"));


}



/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void kmer_quasi_indexer::create_quasi_dictionary (int fingerprint_size){
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
	bool exists;
	u_int32_t nbSequences = 0;
    IBank* bank = Bank::open (getInput()->getStr(STR_URI_BANK_INPUT));
    cout<<"Index "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;
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
void kmer_quasi_indexer::parse_query_sequences (int threshold){

	u_int32_t read_id = 0;
    IBank* bank = Bank::open (getInput()->getStr(STR_URI_QUERY_INPUT));
    cout<<"Query "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;

    ofstream  outFile;
    outFile.open (getInput()->getStr(STR_OUT_FILE));

    outFile << "#query_read_id [target_read_id number_shared_kmers]*"<<endl;

    bool exists;
	vector<u_int32_t> associated_read_ids;

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

  	  vector<u_int32_t>similar_read_ids;



        // We set the data from which we want to extract kmers.
        itKmer.setData (itSeq->getData());
        // We iterate the kmers.
        for (itKmer.first(); !itKmer.isDone(); itKmer.next())
        {
        	quasiDico.get_value(itKmer->value().getVal(),exists,associated_read_ids);
        	similar_read_ids.insert(similar_read_ids.end(),associated_read_ids.begin(),associated_read_ids.end());

//        	cout<<"read "<<nbSequences<<" is associated to "<<endl;
//			cout<<"FOR THIS KMER"<<endl;
//        	for (auto &element:associated_read_ids){
//        		cout<<" "<<element<<endl;
//        		if (element==nbSequences) continue;
//        	}
//			cout<<"FOR ALL KMER"<<endl;
//        	for (auto &element:similar_read_ids){
//        		cout<<" "<<element<<endl;
//        		if (element==nbSequences) continue;
//        	}

        }
        sort(similar_read_ids.begin(),similar_read_ids.end());
        outFile<<read_id<<" ";
        u_int32_t prev_id=-1;
        int nb_shared=0;
        for (auto &element:similar_read_ids){
        	if(prev_id != element){
        		if(prev_id!=-1){
        			if(nb_shared>threshold)
        				outFile<<"["<<prev_id<<" "<<nb_shared<<"] ";
        		}
        		nb_shared=0;
        		prev_id=element;
        	}
        	nb_shared++;
        }
        outFile<<endl;

        read_id++;
    }


    outFile.close();
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

    int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
    if (getInput()->getStr(STR_URI_BANK_INPUT) == getInput()->getStr(STR_URI_QUERY_INPUT))
    	fingerprint_size=0;
    cout<<"fingerprint = "<<fingerprint_size<<endl;
    create_quasi_dictionary(fingerprint_size);
    fill_quasi_dictionary();

    int threshold = getInput()->getInt(STR_THRESHOLD);
    parse_query_sequences(threshold-1); //-1 avoids >=

    getInfo()->add (1, &LibraryInfo::getInfo());
    getInfo()->add (1, "input");
    getInfo()->add (2, "Reference bank:",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
    getInfo()->add (2, "Query bank:",  "%s",  getInput()->getStr(STR_URI_QUERY_INPUT).c_str());
    getInfo()->add (2, "Fingerprint size:",  "%d",  fingerprint_size);
    getInfo()->add (2, "Threshold size:",  "%d",  threshold);
    getInfo()->add (1, "output");
    getInfo()->add (2, "Results written in:",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());




}

