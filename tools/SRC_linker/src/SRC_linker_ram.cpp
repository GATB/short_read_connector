#include <SRC_linker_ram.hpp>


using namespace std;


/********************************************************************************/


// We define some constant strings for names of command line parameters
static const char* STR_URI_BANK_INPUT               = "-bank";
static const char* STR_URI_QUERY_INPUT              = "-query";
static const char* STR_FINGERPRINT                  = "-fingerprint_size";
static const char* STR_GAMMA                        = "-gamma";
static const char* STR_WINDOWS_SIZE                 = "-windows_size";
static const char* STR_THRESHOLD                    = "-kmer_threshold";
static const char* STR_OUT_FILE                     = "-out";
static const char* STR_CORE                         = "-core";
static const char* STR_COMMET_LIKE                  = "-no_sharing_detail";
static const char* STR_KEEP_LOW_COMPLEXITY          = "-keep_low_complexity";
static const char* STR_ZERO_DENSITY_WINDOWS_SIZE    = "-zero_density_windows_size";
static const char* STR_ZERO_DENSITY_THRESHOLD       = "-zero_density_threshold";


SRC_linker_ram::SRC_linker_ram ()  : Tool ("SRC_linker_ram"){
	// We add some custom arguments for command line interface
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH,                  "graph input",      true));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT,             "bank input",       true));
	getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT,            "query input",      true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE,                   "output_file",      true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD,                  "Minimal percentage of shared kmer span for considering 2 reads as similar.  The kmer span is the number of bases from the read query covered by a kmer shared with the target read. If a read of length 80 has a kmer-span of 60 with another read from the bank (of unkonwn size), then the percentage of shared kmer span is 75%. If a least a windows (of size \"windows_size\" contains at least kmer_threshold percent of positionf covered by shared kmers, the read couple is conserved).",    false, "75"));
	getParser()->push_back (new OptionOneParam (STR_WINDOWS_SIZE, "size of the window. If the windows size is zero (default value), then the full read is considered",    false, "0"));
    getParser()->push_back (new OptionOneParam (STR_GAMMA,                      "gamma value",      false,  "2"));
	getParser()->push_back (new OptionOneParam (STR_FINGERPRINT,                "fingerprint size", false,  "8"));
    getParser()->push_back (new OptionOneParam (STR_CORE,                       "Number of thread(s)", false,  "1"));
    getParser()->push_back (new OptionNoParam  (STR_COMMET_LIKE,                "Output ids of reads from query input that are shared with at least one read from reference bank input. With this option no information with whom a read is shared is provided, one only knows that a read is shared.", false));
    
    getParser()->push_back (new OptionNoParam  (STR_KEEP_LOW_COMPLEXITY,        "Conserve low complexity sequences during indexing and querying", false));
    getParser()->push_back (new OptionOneParam (STR_ZERO_DENSITY_WINDOWS_SIZE,  "If defined (>0): two reads are linked if they DO NOT contain a window of this size, with a percentage of zero higher than \"-zero_density_threshold\". Note: this test is performed over the full read length, not limited to \"-windows_size\"", false,  "0"));
    getParser()->push_back (new OptionOneParam (STR_ZERO_DENSITY_THRESHOLD,     "See \"-zero_density_windows_size\"", false, "80"));
}


void SRC_linker_ram::create_quasi_dictionary (int fingerprint_size, int nbCores){
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
	if(nbSolidKmers==0){
		cout<<"No solid kmers in bank -- exit"<<endl;
		exit(0);
	}
	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());
	quasiDico = quasidictionaryVectorKeyGeneric<IteratorKmerH5Wrapper, u_int32_t> (nbSolidKmers, iteratorOnKmers, fingerprint_size, nbCores, gamma_value);
}


struct FunctorIndexer{
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t > &quasiDico;
	int                                                                 kmer_size;
    bool                                                                keep_low_complexity;
    Kmer<KMER_SPAN(1)>::ModelCanonical                                  model;

    FunctorIndexer(quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >& quasiDico, int kmer_size, bool keep_low_complexity)  :
        quasiDico(quasiDico),
        kmer_size(kmer_size),
        keep_low_complexity(keep_low_complexity) {
        model = Kmer<KMER_SPAN(1)>::ModelCanonical(kmer_size);
	}

	void operator() (Sequence& seq){
//        read_id++;          // we do not use the seq.getIndex() id as it is limited to a read file and not a read set.
		if(not keep_low_complexity and not is_high_complexity(seq,kmer_size)){return;}
        Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);
		itKmer.setData (seq.getData());
        
//        if(repeated_kmers(model, itKmer)){return;}
        u_int32_t read_id = static_cast<u_int32_t>(seq.getIndex());
//        if (cur_read_id>read_id) read_id = cur_read_id;
//        cout<<" indexing seq "<<read_id<<endl;
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			// Adding the read id to the list of ids associated to this kmer.note that the kmer may not exist in the dictionary if it was under the solidity threshold.in this case, nothing is done
			quasiDico.set_value((itKmer)->value().getVal(), read_id);
		}
	}
};


void SRC_linker_ram::fill_quasi_dictionary (const int nbCores){
	bool exists;
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_BANK_INPUT));
	cout<<"Index "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	Dispatcher dispatcher (nbCores, 10000);
	dispatcher.iterate (itSeq, FunctorIndexer(quasiDico, kmer_size, keep_low_complexity));
}



class FunctorQuerySpanKmers // FunctorQuery used after claires discussion: number of positions covered by a shared kmer.
{
public:
	ISynchronizer*                                                      synchro;
	FILE*                                                               outFile;
	int                                                                 kmer_size;
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t>* quasiDico;
    int                                                                 threshold;
    int                                                                 windows_size;
    bool                                                                commet_like;
	vector<u_int32_t>                                                   associated_read_ids;
    Kmer<KMER_SPAN(1)>::ModelCanonical                                  model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator*                       itKmer;
    int                                                                 zero_density_windows_size;
    int                                                                 zero_density_threshold;
    bool                                                                keep_low_complexity;
    
	FunctorQuerySpanKmers(const FunctorQuerySpanKmers& lol)
	{
		synchro                     =   lol.synchro;
		outFile                     =   lol.outFile;
		kmer_size                   =   lol.kmer_size;
        windows_size                =   lol.windows_size;
		quasiDico                   =   lol.quasiDico;
		threshold                   =   lol.threshold;
		associated_read_ids         =   lol.associated_read_ids;
		model                       =   lol.model;
        commet_like                 =   lol.commet_like;
		itKmer                      =   new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
        zero_density_windows_size   =   lol.zero_density_windows_size;
        zero_density_threshold      =   lol.zero_density_threshold;
        keep_low_complexity         =   lol.keep_low_complexity;
        
	}
    
	FunctorQuerySpanKmers (ISynchronizer* synchro,
                           FILE* outFile,
                           const int kmer_size,
                           quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper,
                           u_int32_t >* quasiDico,
                           const int threshold,
                           const int windows_size,
                           const bool commet_like,
                           const int zero_density_windows_size,
                           const int zero_density_threshold,
                           const bool keep_low_complexity)
	:
    synchro                     (synchro),
    outFile                     (outFile),
    kmer_size                   (kmer_size),
    quasiDico                   (quasiDico),
    threshold                   (threshold),
    windows_size                (windows_size),
    commet_like                 (commet_like),
    zero_density_windows_size   (zero_density_windows_size),
    zero_density_threshold      (zero_density_threshold),
    keep_low_complexity         (keep_low_complexity)
    {
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
		// itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}
    
	FunctorQuerySpanKmers () {
	}
    
	void operator() (Sequence& seq){
        bool optimization_dont_check_all=false;                         // TODO CREATE AN OPTION IF WE KEEP THIS FEATURE
        std::unordered_set<u_int32_t>                                             similar_read_ids ; // conserve read ids that share enought similarity
        
        
        int used_windows_size=windows_size;
        if (windows_size==0){                                           // if windows size == 0 then we use the full read length as windows
            used_windows_size=seq.getDataSize();
        }
        if(not keep_low_complexity and not is_high_complexity(seq, kmer_size)){return;}
		bool exists;
        std::unordered_map<u_int32_t, vector<bool>>                         similar_read_ids_position ; // each read id --> vector of positions covered by a shared kmer
 		similar_read_ids_position={};                                   // tmp list of couples <last used position, kmer spanning>
		itKmer->setData (seq.getData());
		
        u_int i=0; // position on the read
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
			quasiDico->get_value((*itKmer)->value().getVal(),exists,associated_read_ids);
			if(!exists) {++i;continue;}
			for(auto &read_id: associated_read_ids){
				std::unordered_map<u_int32_t, vector<bool>>::const_iterator element = similar_read_ids_position.find(read_id);
                vector<bool> position_shared;
				if(element == similar_read_ids_position.end()) {// not inserted yet create an empty vector
                    position_shared = vector<bool>(seq.getDataSize());
                    for (int pos=i;pos<seq.getDataSize();pos++) position_shared[pos]=false;
                    
                }
                else{
                    position_shared = element->second;
                }
                
                
                for (int pos=i;pos<i+kmer_size && pos<=seq.getDataSize();pos++){
                    position_shared[pos]=true;
                }
                similar_read_ids_position[read_id] = position_shared;
                if (optimization_dont_check_all){
                    int start = i+kmer_size-used_windows_size;
                    if(start<0) start=0;
                    // if this target read was not already found as similar and if the last reachable window has higher score than the threshold, then add this target read id in the list of similar reads.
                    if(similar_read_ids.count(read_id) == 0 && 100*max_populated_window(position_shared, used_windows_size, start>=0?start:0, i+kmer_size)/float(used_windows_size)>=threshold)
                    {
                        similar_read_ids.insert(read_id);
//                        cout<< "AJOUT de "<<read_id<<endl;
                    }
                }
                
            }
            ++i;
		}
        
        if (optimization_dont_check_all){
            print_read_similarities_dont_check_all (seq,  similar_read_ids);
        }
        else{
            if (commet_like)    print_read_similarities_commet_like (seq, used_windows_size, similar_read_ids_position);
            else                print_read_similarities             (seq, used_windows_size, similar_read_ids_position);
        }
    }
private:
    void print_read_similarities_dont_check_all (Sequence&  seq,  std::unordered_set<u_int32_t> similar_read_ids ){
        if (similar_read_ids.empty()) return;
        string toPrint;
        toPrint=to_string(seq.getIndex())+":";
        for (auto &target_read_id:similar_read_ids)
            toPrint+=to_string(target_read_id)+" ";
        toPrint+="\n";
        synchro->lock();
        fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
        synchro->unlock ();
        
    }
    
    void print_read_similarities_commet_like (Sequence& seq, const int used_windows_size, const std::unordered_map<u_int32_t, vector<bool>>  similar_read_ids_position){
        for (auto &matched_read:similar_read_ids_position){
            
            if (zero_density_windows_size > 0 && contains_high_zero_density_windows(matched_read.second)){ continue; }
            const int mpw = max_populated_window(matched_read.second,used_windows_size);
            const float percentage_span_kmer = 100*mpw/float(used_windows_size);
            if (percentage_span_kmer >= threshold) {
                string toPrint  = to_string(seq.getIndex())+"\n";
                synchro->lock();
                fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
                synchro->unlock ();
                return;
            }
        }
    }
    
    void print_read_similarities(Sequence& seq, const int used_windows_size, const std::unordered_map<u_int32_t, vector<bool>>  similar_read_ids_position){
        string toPrint;
        bool read_id_printed=false; // Print (and sync file) only if the read is similar to something.
        for (auto &matched_read:similar_read_ids_position){
            
            if (zero_density_windows_size > 0 && contains_high_zero_density_windows(matched_read.second)){ continue; }
            const int mpw = max_populated_window(matched_read.second,used_windows_size);
            const float percentage_span_kmer = 100*mpw/float(used_windows_size);
            
            if (percentage_span_kmer >= threshold) {
                if (not read_id_printed){
                    read_id_printed=true;
                    //					synchro->lock();
                    toPrint=to_string(seq.getIndex())+":";
                }
                toPrint+=to_string(matched_read.first)+"-"+to_string(mpw)+"-"+to_string(float(percentage_span_kmer))+" ";
            }
            
        }
        if(read_id_printed){
            toPrint+="\n";
            synchro->lock();
            fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
            synchro->unlock ();
        }

    }
    
    
    bool contains_high_zero_density_windows(const vector<bool> populated){
        const int max_number_of_zeros = zero_density_windows_size*zero_density_threshold/100;
        
        int number_zeros=0;
        for(int i=0;i<zero_density_windows_size;i++){
            if (populated[i]==0) {
                number_zeros++;
            }
        }
        
        if (number_zeros>max_number_of_zeros)   return true;
        
        const int size_vector = populated.size();
        const int last_excluded_starting_window_position = size_vector-zero_density_windows_size+1;
        
        
        for (int pos=1;pos<last_excluded_starting_window_position;pos+=1){
            if (populated[pos-1]==0){number_zeros--;}
            if (populated[pos+windows_size-1]==0){number_zeros++;}
            if (number_zeros>max_number_of_zeros)   return true;
        }
        
        return false;
    }

    
    int max_populated_window(const vector<bool> populated, const int windows_size, const int start=0, int stop=-1){
        
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
//        cout<<max<<" "<<windows_size<<endl; //DEBUG
        return max;
    }
    
};



void SRC_linker_ram::parse_query_sequences (int threshold, const int nbCores, const int windows_size, const bool commet_like){
    BankAlbum banks (getInput()->getStr(STR_URI_QUERY_INPUT));
    const std::vector<IBank*>& banks_of_queries = banks.getBanks();
    const int number_of_read_sets = banks_of_queries.size();
    
    
	cout<<"Query "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_QUERY_INPUT)<<endl;
	FILE * outFile;
	outFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
    if (commet_like){
        string message("#query_read_id \n#Target read set: "+getInput()->getStr(STR_URI_BANK_INPUT)+"\n");
        fwrite((message).c_str(), sizeof(char), message.size(), outFile);
    }
    else{
        string message("#query_read_id [target_read_id-kmer_span (k="+to_string(kmer_size)+")-kmer_span query percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)\n"+"#Target read set: "+getInput()->getStr(STR_URI_BANK_INPUT)+"\n");
        fwrite((message).c_str(), sizeof(char), message.size(), outFile);
    }
    
    
    for( int bank_id=0;bank_id<number_of_read_sets;bank_id++){ // iterate each bank
        
        
        IBank* bank=banks_of_queries[bank_id];
        LOCAL (bank);
        string message("#Query read set number "+bank->getId()+"\n");
        fwrite((message).c_str(), sizeof(char), message.size(), outFile);
        string progressMessage("Querying read set "+bank->getId());
        ProgressIterator<Sequence> itSeq (*bank, progressMessage.c_str());
        ISynchronizer* synchro = System::thread().newSynchronizer();
        Dispatcher dispatcher (nbCores, 1000);
        dispatcher.iterate (itSeq, FunctorQuerySpanKmers(synchro,outFile, kmer_size,&quasiDico, threshold, windows_size, commet_like,zero_density_windows_size,zero_density_threshold, keep_low_complexity));
        delete synchro;
    }
	fclose (outFile);
    
    
//	IBank* bank = Bank::open (getInput()->getStr(STR_URI_QUERY_INPUT));
//	cout<<"Query "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_QUERY_INPUT)<<endl;
//	FILE * pFile;
//	pFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
//	string message("#query_read_id [target_read_id-kmer_span (k="+to_string(kmer_size)+")-kmer_span query percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)\n");
//	fwrite((message).c_str(), sizeof(char), message.size(), pFile);
//	LOCAL (bank);
//	ProgressIterator<Sequence> itSeq (*bank);
//	ISynchronizer* synchro = System::thread().newSynchronizer();
//	Dispatcher dispatcher (nbCores, 10000);
//	dispatcher.iterate (itSeq, FunctorQuerySpanKmers(synchro,pFile, kmer_size,&quasiDico, threshold));
//	fclose (pFile);
//	delete synchro;
}


void SRC_linker_ram::execute (){
    
	int nbCores                     = getInput()->getInt(STR_CORE);
	int fingerprint_size            = getInput()->getInt(STR_FINGERPRINT);
    gamma_value                     = getInput()->getInt(STR_GAMMA);
    zero_density_windows_size       = getInput()->getInt(STR_ZERO_DENSITY_WINDOWS_SIZE);
    zero_density_threshold          = getInput()->getInt(STR_ZERO_DENSITY_THRESHOLD);
	// IMPORTANT NOTE:
	// Actually, during the filling of the dictionary values, one may fall on non solid non indexed kmers
	// that are quasi dictionary false positives (ven with a non null fingerprint. This means that one nevers knows in advance how much
	// values are gonna be stored for all kmers. This is why I currently us a vector<u_int32_t> for storing read ids associated to a kmer.

	// We need a non null finger print because of non solid non indexed kmers
	//	if (getInput()->getStr(STR_URI_BANK_INPUT).compare(getInput()->getStr(STR_URI_QUERY_INPUT))==0)
	//		fingerprint_size=0;

    int threshold            = getInput()->getInt(STR_THRESHOLD);
    int windows_size         = getInput()->getInt(STR_WINDOWS_SIZE);
    bool commet_like         = getInput()->get(STR_COMMET_LIKE)>0?true:false;
    keep_low_complexity = getInput()->get(STR_KEEP_LOW_COMPLEXITY)>0?true:false;
	create_quasi_dictionary(fingerprint_size, nbCores);
	fill_quasi_dictionary(nbCores);

	parse_query_sequences(threshold, nbCores, windows_size, commet_like);

	getInfo()->add (1, &LibraryInfo::getInfo());
	getInfo()->add (1, "input");
	getInfo()->add (2, "Reference bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
    getInfo()->add (2, "Query bank",  "%s",  getInput()->getStr(STR_URI_QUERY_INPUT).c_str());
    getInfo()->add (2, "windows_size size",  "%d",  windows_size);
    getInfo()->add (2, "Kmer size",  "%d",  kmer_size);
	getInfo()->add (2, "Fingerprint size",  "%d",  fingerprint_size);
    getInfo()->add (2, "gamma",  "%d",  gamma_value);
	getInfo()->add (2, "Minimal kmer span percentage",  "%d",  threshold);
    if(keep_low_complexity)
        getInfo()->add (2, "Low complexity sequences were kept");
    else
        getInfo()->add (2, "Low complexity sequences were removed");
        
	getInfo()->add (1, "output");
    if (commet_like)
        getInfo()->add (2, "Output only ids of read shared (no complete links)");
	getInfo()->add (2, "Results written in",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());
}
