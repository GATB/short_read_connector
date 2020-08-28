#include <SRC_counter.hpp>


using namespace std;


/********************************************************************************/


// We define some constant strings for names of command line parameters
static const char* STR_MAKE_INDEX                   = "-make_index";
static const char* STR_DUMPED_QD_FILE_NAME          = "-dumped_quasi_dict";
static const char* STR_URI_QUERY_INPUT              = "-query";
static const char* STR_FINGERPRINT                  = "-fingerprint_size";
static const char* STR_GAMMA                        = "-gamma";
static const char* STR_KEEP_LOW_COMPLEXITY          = "-keep_low_complexity";
static const char* STR_OUT_FILE                     = "-out";
static const char* STR_CORE                         = "-core";
static const char* STR_THRESHOLD                    = "-coverage_threshold";
static const char* STR_WINDOWS_SIZE                 = "-windows_size";


SRC_counter::SRC_counter ()  : Tool ("SRC_counter"){
	// INDEX OR QUERY
	getParser()->push_back (new OptionNoParam (STR_MAKE_INDEX,                 "Create and dump the index",      false));

    // OPTIONS FOR CREATE INDEX AND PERFORM QUERIES
    getParser()->push_back (new OptionOneParam (STR_DUMPED_QD_FILE_NAME,        "Mandatory for creating the index (creating) and for performing queries (read only). Dumped index file name",      true));
    getParser()->push_back (new OptionOneParam (STR_CORE,                       "Indexing & querying: Number of thread(s)", false,  "1"));
    getParser()->push_back (new OptionNoParam  (STR_KEEP_LOW_COMPLEXITY,        "Indexing & querying: Conserve low complexity sequences during indexing and querying", false));

    // OPTIONS FOR CREATE INDEX
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH,          "graph input",   false));
	getParser()->push_back (new OptionOneParam (STR_GAMMA,              "gamma value",    false, "2"));
    getParser()->push_back (new OptionOneParam (STR_FINGERPRINT,        "fingerprint size",    false, "8"));


    // OPTIONS FOR QUERY
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE,           "output_file",    false));
    getParser()->push_back (new OptionOneParam (STR_THRESHOLD,          "Minimal percentage of shared kmer span for considering a query read as similar to a read set.  The kmer span is the number of bases from the read query covered by a kmer shared with the target bank read set. If a read of length 80 has a kmer-span of 60 with the bank, then the percentage of shared kmer span is 75%. If a least a windows (of size \"windows_size\" contains at least kmer_threshold percent of positionf covered by shared kmers, the read is output in the boolean vector).",    false, "50"));
    getParser()->push_back (new OptionOneParam (STR_WINDOWS_SIZE, 		"size of the window. If the windows size is zero (default value), then the full read is considered",    false, "0"));
	getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT,    "query input",    false));
}


// We define a functor that will be cloned by the dispatcher
struct FunctorIndexer
{
	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char > &quasiDico;
    int kmer_size;
    bool keep_low_complexity;
    
	FunctorIndexer(quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char >& quasiDico, int kmer_size, bool keep_low_complexity)  :  quasiDico(quasiDico), kmer_size(kmer_size), keep_low_complexity(keep_low_complexity) {}

	void operator() (Kmer<>::Count & itKmer){
//        if(not keep_low_complexity and not is_high_complexity(itKmer.getValue()),kmer_size)){return;}
		quasiDico.set_value(itKmer.value.getVal(), itKmer.abundance>0xFF?0xFF:(unsigned char)itKmer.abundance);
	}
};


void SRC_counter::create_and_fill_quasi_dictionary (){
//	const int display = getInput()->getInt (STR_VERBOSE);
	// We get a handle on the HDF5 storage object.
	// Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
	auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (getInput()->getStr(STR_URI_GRAPH)));
	// We get the group for dsk
	Group& dskGroup = storage->getGroup("dsk");
	kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());
	// We get the solid kmers collection 1) from the 'dsk' group  2) from the 'solid' collection
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");
	nbSolidKmers = solidKmers.getNbItems();
	if(nbSolidKmers==0){cout<<"No solid kmers in bank -- exit"<<endl;exit(0);}
	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());
	quasiDico = quasidictionaryKeyGeneric<IteratorKmerH5Wrapper, unsigned char> (nbSolidKmers, iteratorOnKmers, fingerprint_size, nbCores, gamma_value);
    
    cout<<"Empty quasi-ictionary memory usage (MB) = "<<System::info().getMemorySelfUsed()/1024<<endl;
    
    
	ProgressIterator<Kmer<>::Count> itKmers (solidKmers.iterator(), "Indexing solid kmers counts", nbSolidKmers);
	Dispatcher dispatcher (nbCores, 10000);
	dispatcher.iterate (itKmers, FunctorIndexer(quasiDico, kmer_size, keep_low_complexity));
    
    cout<<"Filled quasi-ictionary memory usage (MB) = "<<System::info().getMemorySelfUsed()/1024<<endl;
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
	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char>* quasiDico;
	vector<u_int32_t> associated_read_ids;
	Kmer<KMER_SPAN(1)>::ModelCanonical model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator* itKmer;
    bool keep_low_complexity;
    int threshold;
    int windows_size;
//    BooleanVector* bv;

	FunctorQuery(const FunctorQuery& lol)
	{
		synchro             = lol.synchro;
		outFile             = lol.outFile;
		kmer_size           = lol.kmer_size;
		quasiDico           = lol.quasiDico;
		associated_read_ids = lol.associated_read_ids;
		model               = lol.model;
		itKmer              = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
        keep_low_complexity = lol.keep_low_complexity;
        threshold           = lol.threshold;
        windows_size        = lol.windows_size;
//        bv                  = lol.bv;
	}


	FunctorQuery (ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char >* quasiDico,  const int keep_low_complexity, int threshold, int windows_size)//, BooleanVector * bv)
    : synchro(synchro), outFile(outFile), kmer_size(kmer_size), quasiDico(quasiDico), keep_low_complexity(keep_low_complexity), threshold(threshold), windows_size(windows_size){//, bv(bv) {
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
		// itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}


	~FunctorQuery () {
	}

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
		const size_t n = v.size() / 2;
		nth_element(v.begin(), v.begin()+n, v.end());
		median=v[n];
		return true;


	}
    
    int number_positions_covered_shared_kmer(vector<int> &v, const int size_seq){
        
        if(v.size()==0) return 0;
        int number_positions_covered    = 0;
        int next_uncovered_position     = -1;
        std::vector<int>::iterator it   = v.begin();
        for(int i=0; i<size_seq;i++){
            if (i==*it){
                next_uncovered_position=i+kmer_size;
                if(it+1 != v.end()) it++;
            }
            if (i<next_uncovered_position) number_positions_covered++;
        }
        return number_positions_covered;
    }

	void operator() (Sequence& seq){
        int used_windows_size=windows_size;
        if (windows_size==0 || windows_size>seq.getDataSize()){                                           // if windows size == 0 then we use the full read length as windows
            used_windows_size=seq.getDataSize();
        }
        if(not keep_low_complexity and not is_high_complexity(seq,kmer_size)){return;}
		bool exists;
		unsigned char count;
		itKmer->setData (seq.getData());
		vector<int> values;                                                                             // For each position: number of occurrences of the kmer starting at this position.
//        vector<int> covered_positions; // DEPRECATED [OLD WAY FOR COMPUTING SHARE KMER POSITIONS, FASTER BUT NON IMPLEMENTED WITH WINDOWS SIZE METHODS (max_populated_window)]
        vector<bool> position_shared =  vector<bool>(seq.getDataSize());                                // boolean vector. A position is true if it's covered by at least a shared kmer.
        for (int pos=0;pos<seq.getDataSize();pos++) position_shared[pos]=false;
        
        int position=0;
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
			quasiDico->get_value((*itKmer)->value().getVal(),exists,count);
			if(!exists) {
                count=0;
            }
            values.push_back(count);
//            if (count>0) covered_positions.push_back(position); // DEPRECATED [OLD WAY FOR COMPUTING SHARE KMER POSITIONS, FASTER BUT NON IMPLEMENTED WITH WINDOWS SIZE METHODS (max_populated_window)]
            if (count>0) { // TODO: OPTIMIZABLE.
                for (int pos=position;pos<position+kmer_size && pos<=seq.getDataSize();pos++) position_shared[pos]=true;
            }
            position++;
		}
        
//        float percentage_shared_positions = 100*number_positions_covered_shared_kmer(covered_positions, seq.getDataSize())/float(seq.getDataSize());  // DEPRECATED [OLD WAY FOR COMPUTING SHARE KMER POSITIONS, FASTER BUT NON IMPLEMENTED WITH WINDOWS SIZE METHODS (max_populated_window)]
        const int mpw = max_populated_window(position_shared,used_windows_size);
        const float percentage_span_kmer = 100*mpw/float(used_windows_size);
        
//        if (percentage_shared_positions !=percentage_span_kmer){cout<<percentage_shared_positions<< " == " <<percentage_span_kmer<<" ?"<<endl; exit(1);} // TO REMOVE
        
        
		float mean;
		int median, min, max;
		if(mean_median_min_max(values, mean, median, min, max)){
           
			string toPrint (to_string(seq.getIndex())+" "+to_string(mean)+" "+to_string(median)+" "+to_string(min)+" "+to_string(max)+" "+to_string(percentage_span_kmer));
            
            
            toPrint.append("\n");

            if (percentage_span_kmer>=threshold) {
                synchro->lock();
                //                bv->set(seq.getIndex());
                fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
                synchro->unlock ();

            }
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

unsigned long get_bank_nb_items(IBank* bank){
    if(bank->getNbItems() != -1) return bank->getNbItems();
    // else: non implemented
    Iterator<Sequence>* it = bank->iterator();
    LOCAL (it);
    // We loop over sequences.
    unsigned long size=0;
    for (it->first(); !it->isDone(); it->next()){++size;}
    return size;
}

void SRC_counter::parse_query_sequences (){
    std::string bank_filename = getInput()->getStr(STR_DUMPED_QD_FILE_NAME).substr(getInput()->getStr(STR_DUMPED_QD_FILE_NAME).find_last_of("/\\") + 1);

    if (getInput()->getStr(STR_URI_QUERY_INPUT) == ""){
        cerr<<STR_URI_QUERY_INPUT<<" missing "<<endl;
        exit(0);
    }
    BankAlbum banks (getInput()->getStr(STR_URI_QUERY_INPUT));

    const std::vector<IBank*>& banks_of_queries = banks.getBanks();
    const int number_of_read_sets = banks_of_queries.size();
    
	FILE * pFile;

    if (getInput()->getStr(STR_OUT_FILE) == ""){
        cerr<<STR_OUT_FILE<<" missing "<<endl;
        exit(0);
    }
	pFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
    
    
	cout<<"Query "<<kmer_size<<"-mers from "<<getInput()->getStr(STR_URI_QUERY_INPUT)<<endl;
    for( int bank_id=0;bank_id<number_of_read_sets;bank_id++){ // iterate each bank
        
        IBank* bank=banks_of_queries[bank_id];
        LOCAL (bank);

        
        string message("#query_read_id (from bank "+bank->getId()+") mean median min max percentage_shared_positions -- number of shared "+to_string(kmer_size)+"mers with banq "+getInput()->getStr(STR_DUMPED_QD_FILE_NAME)+"\n");
        fwrite((message).c_str(), sizeof(char), message.size(), pFile);
        string progressMessage("Querying read set "+bank->getId());
        ProgressIterator<Sequence> itSeq (*bank, progressMessage.c_str());
        ISynchronizer* synchro = System::thread().newSynchronizer();
        Dispatcher dispatcher (nbCores, 10000);
        dispatcher.iterate (itSeq, FunctorQuery(synchro,pFile, kmer_size,&quasiDico, keep_low_complexity, threshold, windows_size));//, &bv));
        delete synchro;
        std::string query_filename = bank->getId().substr(bank->getId().find_last_of("/\\") + 1);
//        cout<<bv.nb_one()<<" reads in out_"+query_filename+"_in_"+bank_filename+".bv"<<endl;
//        bv.print("out_"+query_filename+"_in_"+bank_filename+".bv");
    }
    fclose (pFile);
}


void SRC_counter::load_quasi_dictionary(std::string dumped_file_name){
    std::ifstream is = std::ifstream(dumped_file_name, std::ios_base::binary); 
    is.read(reinterpret_cast<char*>(&this->kmer_size), sizeof(this->kmer_size));
	quasiDico = quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, unsigned char >();
    quasiDico.load(is);
}



void SRC_counter::write_quasi_dictionary(std::string dumped_file_name){
    std::ofstream os = std::ofstream(dumped_file_name, std::ios_base::binary); 
    os.write(reinterpret_cast<char const*>(&this->kmer_size), sizeof(this->kmer_size));
    quasiDico.save(os);
}
void SRC_counter::execute (){


    std::string dumped_file_name    = getInput()->getStr(STR_DUMPED_QD_FILE_NAME);
    if (dumped_file_name == ""){
        cerr<<STR_DUMPED_QD_FILE_NAME<<" missing "<<endl;
        exit(0);
    }


	nbCores = getInput()->getInt(STR_CORE);
    keep_low_complexity = getInput()->get(STR_KEEP_LOW_COMPLEXITY)?true:false;

	// IMPORTANT NOTE:
	// Actually, during the filling of the dictionary values, one may fall on non solid non indexed kmers
	// that are quasi dictionary false positives (even with a non null fingerprint. This means that one nevers knows in advance how much
	// values are gonna be stored for all kmers. This is why I currently use a vector<u_int32_t> for storing read ids associated to a kmer.



    bool make_index                 	= getInput()->get(STR_MAKE_INDEX)?true:false;
    // MAKE INDEX:
    if (make_index){
	    fingerprint_size                = getInput()->getInt(STR_FINGERPRINT);
        gamma_value                     = getInput()->getInt(STR_GAMMA);


		create_and_fill_quasi_dictionary();
        write_quasi_dictionary(dumped_file_name); 

		getInfo()->add (1, &LibraryInfo::getInfo());
        getInfo()->add (1, "input");

        getInfo()->add (2, "Kmer size",  "%d",  kmer_size);
        getInfo()->add (2, "Fingerprint size",  "%d",  fingerprint_size);
        getInfo()->add (2, "gamma",  "%d",  gamma_value);
        if(keep_low_complexity)
            getInfo()->add (2, "Low complexity sequences were kept");
        else
            getInfo()->add (2, "Low complexity sequences were removed");
            
        getInfo()->add (1, "output");
        getInfo()->add (2, "Index dumped in",  "%s",  dumped_file_name.c_str());
    }
	else{
   		windows_size         = getInput()->getInt(STR_WINDOWS_SIZE);
		// output named 		  getIput()->getStr(STR_OUT_FILE) is generated in the "parse_query_sequences" method
		// input queries 		  getIput()->getStr(STR_URI_QUERY_INPUT) is generated in the "parse_query_sequences" method
    	threshold 						= getInput()->getInt(STR_THRESHOLD);

		load_quasi_dictionary(dumped_file_name); 
		parse_query_sequences();

		getInfo()->add (1, &LibraryInfo::getInfo());
		getInfo()->add (1, "input");
        getInfo()->add (2, "Used index",  "%s",  getInput()->getStr(STR_DUMPED_QD_FILE_NAME).c_str());
		getInfo()->add (2, "Query bank",  "%s",  getInput()->getStr(STR_URI_QUERY_INPUT).c_str());
		getInfo()->add (2, "Kmer size",  "%d",  kmer_size);
		getInfo()->add (2, "windows_size", "%d", windows_size);
		getInfo()->add (2, "Minimal kmer span percentage", "%d", threshold);
		if(keep_low_complexity)
			getInfo()->add (2, "Low complexity query sequences were kept");
		else
			getInfo()->add (2, "Low complexity query sequences were removed");
		
		getInfo()->add (1, "output");
		getInfo()->add (2, "Results written in",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());
	}

}
