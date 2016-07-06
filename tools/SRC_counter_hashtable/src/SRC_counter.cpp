#include <SRC_counter.hpp>


using namespace std;


/********************************************************************************/


// We define some constant strings for names of command line parameters
static const char* STR_URI_BANK_INPUT = "-bank";
static const char* STR_URI_QUERY_INPUT = "-query";
static const char* STR_FINGERPRINT = "-fingerprint_size";
static const char* STR_GAMMA = "-gamma";
static const char* STR_THRESHOLD = "-kmer_threshold";
static const char* STR_OUT_FILE = "-out";
static const char* STR_CORE = "-core";


SRC_counter::SRC_counter ()  : Tool ("SRC_counter"){
	// We add some custom arguments for command line interface
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
	//getParser()->push_back (new OptionOneParam (STR_VERBOSE,   "verbosity (0:no display, 1: display kmers, 2: display distrib",  false, "0"));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
	getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT, "query input",    true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal number of shared kmers for considering 2 reads as similar",    false, "10"));
	getParser()->push_back (new OptionOneParam (STR_GAMMA, "gamma value - useless here",    false, "2"));
	getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "8"));
	getParser()->push_back (new OptionOneParam (STR_CORE, "Number of thread",    false, "1"));
}


// We define a functor that will be cloned by the dispatcher
struct FunctorIndexer
{
	ISynchronizer* synchro;
//	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char > &quasiDico;
	std::unordered_map<u_int64_t, unsigned char> &hashDico;
	int kmer_size;
    
	FunctorIndexer(ISynchronizer* synchro, std::unordered_map<u_int64_t, unsigned char>& hashDico , int kmer_size)  :  synchro(synchro), hashDico(hashDico), kmer_size(kmer_size) {
	}

	void operator() (Kmer<>::Count & itKmer){
        synchro->lock();
       
		hashDico[itKmer.value.getVal()] = itKmer.abundance>0xFF?0xFF:(unsigned char)itKmer.abundance;//        hashDico[itKmer.value.getVal()] = 42;
//        cout<<"psdoc"<<endl;
        
        synchro->unlock();

        //        cout<<"indexed = "<<(int)hashDico[itKmer.value.getVal()]<<"value = "<<itKmer.abundance<<endl;
//        cout<<hashDico.size()<<" "<<i<<endl;
	}
};



void SRC_counter::create_and_fill_quasi_dictionary (int fingerprint_size, const int nbCores){
	const int display = getInput()->getInt (STR_VERBOSE);
	// We get a handle on the HDF5 storage object.
	// Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
	auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (getInput()->getStr(STR_URI_GRAPH)));
	// We get the group for dsk
	Group& dskGroup = storage->getGroup("dsk");
	kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());
	// We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");
	nbSolidKmers = solidKmers.getNbItems();
    hashDico.reserve(nbSolidKmers);
    cout<<"Reserve "<<nbSolidKmers<<" in the hash table"<<endl;
	if(nbSolidKmers==0){
		cout<<"No solid kmers in bank -- exit"<<endl;
		exit(0);
	}
	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());
	ProgressIterator<Kmer<>::Count> itKmers (solidKmers.iterator(), "Indexing solid kmers", nbSolidKmers);
	Dispatcher dispatcher (nbCores, 10000);
    
	ISynchronizer* synchro = System::thread().newSynchronizer();
	dispatcher.iterate (itKmers, FunctorIndexer(synchro, hashDico, kmer_size));
    delete synchro;
    cout<<hashDico.size()<<" kmers are indexed"<<endl;
    cout<<"Filled hash table memory usage (MB) = "<<System::info().getMemorySelfUsed()/1024<<endl;

}

//iterator->item().abundance
//		iterator->item().value.getVal()



// We define a functor that will be cloned by the dispatcher
class FunctorQuery
{
public:
	ISynchronizer* synchro;
	FILE* outFile;
	int kmer_size;
//	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char>* quasiDico;
	std::unordered_map<u_int64_t, unsigned char> *hashDico;
	int threshold;
	vector<u_int32_t> associated_read_ids;
	Kmer<KMER_SPAN(1)>::ModelCanonical model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator* itKmer;

	FunctorQuery(const FunctorQuery& lol)
	{
		synchro=lol.synchro;
		outFile=lol.outFile;
		kmer_size=lol.kmer_size;
		hashDico=lol.hashDico;
		threshold=lol.threshold;
		associated_read_ids=lol.associated_read_ids;
		model=lol.model;
		itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}


	FunctorQuery (ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  std::unordered_map<u_int64_t, unsigned char> *hashDico, const int threshold)
	: synchro(synchro), outFile(outFile), kmer_size(kmer_size), hashDico(hashDico), threshold(threshold) {
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
		// itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}


	~FunctorQuery () {
	}
	//
	//	int median(vector<int> &v)
	//	{
	//		if(v.size()==0) return 0;
	//	    size_t n = v.size() / 2;
	//	    nth_element(v.begin(), v.begin()+n, v.end());
	//	    return v[n];
	//	}
	//
	//	float mean(vector<int> &v){
	//		float mean=0;
	//		for(auto element: v) mean+=element;
	//		return mean/(float)v.size();
	//	}

	bool mean_median_min_max(vector<int> &v, float & mean, int & median, int & min, int & max){
		if(v.size()==0) return false;
		mean=0;
		min=v[0];
		max=v[0];
		for(auto element: v) {
			mean+=element;
			if(element<min) min=element;
			if(element>max) max=element;
		}
		mean= mean/(float)v.size();
		size_t n = v.size() / 2;
		nth_element(v.begin(), v.begin()+n, v.end());
		median=v[n];
		return true;


	}

	void operator() (Sequence& seq){
		if(not valid_sequence(seq, kmer_size)){return;}

		bool exists;
		unsigned char count;
		itKmer->setData (seq.getData());
		vector<int> values;
        
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
			std::unordered_map<u_int64_t, unsigned char>::iterator got = hashDico->find((*itKmer)->value().getVal());
            if (got == hashDico->end()) {
                values.push_back(0);
                continue;}
//            cout<<"in"<<endl; //DEB
            
			values.push_back(got->second);
		}
        
		float mean;
		int median, min, max;
		if(mean_median_min_max(values, mean, median, min, max)){

			string toPrint (to_string(seq.getIndex())+" "+to_string(mean)+" "+to_string(median)+" "+to_string(min)+" "+to_string(max)+"\n");
			synchro->lock();
			fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
			synchro->unlock ();
		}

		else{
			string toPrint (to_string(seq.getIndex())+" none\n");
			synchro->lock();
			fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
			synchro->unlock ();
		}


	}
};


void SRC_counter::parse_query_sequences (int threshold, const int nbCores){
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_QUERY_INPUT));
    cout<<"Total nb indexed elements = "<<hashDico.size()<<endl;
	cout<<"Query "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_QUERY_INPUT)<<endl;
	FILE * pFile;
	pFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
	string message("#query_read_id mean median min max number of shared "+to_string(kmer_size)+"mers with banq read set\n");
	fwrite((message).c_str(), sizeof(char), message.size(), pFile);

	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	ISynchronizer* synchro = System::thread().newSynchronizer();
	Dispatcher dispatcher (nbCores, 10000);
	dispatcher.iterate (itSeq, FunctorQuery(synchro,pFile, kmer_size,&hashDico, threshold));
	fclose (pFile);
	delete synchro;
}


void SRC_counter::execute (){
	int nbCores = getInput()->getInt(STR_CORE);
	int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
	// IMPORTANT NOTE:
	// Actually, during the filling of the dictionary values, one may fall on non solid non indexed kmers
	// that are quasi dictionary false positives (ven with a non null fingerprint. This means that one nevers knows in advance how much
	// values are gonna be stored for all kmers. This is why I currently us a vector<u_int32_t> for storing read ids associated to a kmer.

	// We need a non null finger print because of non solid non indexed kmers
	//	if (getInput()->getStr(STR_URI_BANK_INPUT).compare(getInput()->getStr(STR_URI_QUERY_INPUT))==0)
	//		fingerprint_size=0;
	cout<<"fingerprint = "<<fingerprint_size<<endl;
	create_and_fill_quasi_dictionary(fingerprint_size, nbCores);

	int threshold = getInput()->getInt(STR_THRESHOLD);
	parse_query_sequences(threshold-1, nbCores); //-1 avoids >=

	getInfo()->add (1, &LibraryInfo::getInfo());
	getInfo()->add (1, "input");
	getInfo()->add (2, "Reference bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
	getInfo()->add (2, "Query bank",  "%s",  getInput()->getStr(STR_URI_QUERY_INPUT).c_str());
	getInfo()->add (2, "Fingerprint size",  "%d",  fingerprint_size);
	getInfo()->add (2, "Threshold size",  "%d",  threshold);
	getInfo()->add (1, "output");
	getInfo()->add (2, "Results written in",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());
}
