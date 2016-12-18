#include <SRC_linker_fuzzy.hpp>


using namespace std;


/********************************************************************************/


// We define some constant strings for names of command line parameters
static const char* STR_URI_BANK_INPUT = "-bank";
static const char* STR_URI_QUERY_INPUT = "-query";
static const char* STR_FINGERPRINT = "-fingerprint_size";
static const char* STR_CONTEXT_SIZE = "-context_size";
static const char* STR_GAMMA = "-gamma";
static const char* STR_WINDOWS_SIZE = "-windows_size";
static const char* STR_THRESHOLD = "-kmer_threshold";
static const char* STR_MAX_EDIT_DISTANCE = "-max_edit_distance";
static const char* STR_OUT_FILE = "-out";
static const char* STR_CORE = "-core";



SRC_linker_fuzzy::SRC_linker_fuzzy ()  : Tool ("SRC_linker_fuzzy"){
	// We add some custom arguments for command line interface
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
	getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT, "query input",    true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal percentage of shared kmer span for considering 2 reads as similar.  The kmer span is the number of bases from the read query covered by a kmer shared with the target read. If a read of length 80 has a kmer-span of 60 with another read from the bank (of unkonwn size), then the percentage of shared kmer span is 75%. If a least a windows (of size \"windows_size\" contains at least kmer_threshold percent of positionf covered by shared kmers, the read couple is conserved.",    false, "75"));
	getParser()->push_back (new OptionOneParam (STR_WINDOWS_SIZE, "size of the window. If the windows size is zero (default value), then the full read is considered",    false, "0"));
    getParser()->push_back (new OptionOneParam (STR_GAMMA, "gamma value",    false, "2"));
    getParser()->push_back (new OptionOneParam (STR_CONTEXT_SIZE, "size of each context",    false, "8"));
    getParser()->push_back (new OptionOneParam (STR_MAX_EDIT_DISTANCE, "maximal edit distance (cumulated for both left and right contexts)",    false, "2"));
    getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "8"));
	getParser()->push_back (new OptionOneParam (STR_CORE, "Number of thread",    false, "1"));
}


void SRC_linker_fuzzy::create_quasi_dictionary (int fingerprint_size, int nbCores){
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
	quasiDico = quasidictionaryVectorKeyGeneric<IteratorKmerH5Wrapper,  pair<u_int32_t,string>> (nbSolidKmers, iteratorOnKmers, fingerprint_size, nbCores, gamma_value);
}


struct FunctorIndexer{
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper,  pair<u_int32_t,string> > &quasiDico;
	int kmer_size;
    int context_size;

	FunctorIndexer(quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper,  pair<u_int32_t,string> >& quasiDico, int kmer_size, int context_size)  :  quasiDico(quasiDico), kmer_size(kmer_size), context_size(context_size) {
	}

	void operator() (Sequence& seq){
		if(not valid_sequence(seq,kmer_size)){return;}
        const string string_seq = seq.toString();
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);
		itKmer.setData (seq.getData());
        
//        if(repeated_kmers(model, itKmer)){return;}
		u_int32_t read_id = static_cast<u_int32_t>(seq.getIndex()+1);
//        cout<<string_seq<<endl; //DEBUG
        int position=-1;
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			// Adding the read id to the list of ids associated to this kmer.note that the kmer may not exist in the dictionary if it was under the solidity threshold.in this case, nothing is done
            position++;
            if (position<context_size) continue; // avoid first positions with no context
            if (position>seq.getDataSize()-context_size-kmer_size) break; // avoid last positions with no context
            string context=string_seq.substr (position-context_size,context_size);
            context+=string_seq.substr (position+kmer_size,context_size);
//            cout<< "read "<<read_id<<" position "<<position<<" contexts are "<<context<<endl;
            pair<u_int32_t,string> this_pair (read_id,context);
			quasiDico.set_value((itKmer)->value().getVal(), this_pair);
		}
	}
};


void SRC_linker_fuzzy::fill_quasi_dictionary (const int nbCores){
	bool exists;
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_BANK_INPUT));
	cout<<"Index "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	Dispatcher dispatcher (nbCores, 10000);
	dispatcher.iterate (itSeq, FunctorIndexer(quasiDico, kmer_size, context_size));
}



class FunctorQuerySpanKmers // FunctorQuery used after claires discussion: number of positions covered by a shared kmer.
{
public:
	ISynchronizer* synchro;
    FILE* outFile;
    int kmer_size;
    int context_size;
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, pair<u_int32_t,string>>* quasiDico;
    int threshold;
    int max_edit_distance;
    int windows_size;
	vector<pair<u_int32_t,string>> associated_read_ids_and_contexts;
	std::unordered_map<u_int32_t, vector<bool>> similar_read_ids_position; // each read id --> vector of positions covered by a shared kmer
	Kmer<KMER_SPAN(1)>::ModelCanonical model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator* itKmer;
    
	FunctorQuerySpanKmers(const FunctorQuerySpanKmers& lol)
	{
		synchro=lol.synchro;
		outFile=lol.outFile;
		kmer_size=lol.kmer_size;
        max_edit_distance=lol.max_edit_distance;
        context_size=lol.context_size;
        windows_size = lol.windows_size;
		quasiDico=lol.quasiDico;
		threshold=lol.threshold;
		associated_read_ids_and_contexts=lol.associated_read_ids_and_contexts;
		similar_read_ids_position=lol.similar_read_ids_position;
		model=lol.model;
		itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}
    
	FunctorQuerySpanKmers (ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper,  pair<u_int32_t,string>>* quasiDico, const int threshold, const int windows_size, const int context_size, const int max_edit_distance)
	: synchro(synchro), outFile(outFile), kmer_size(kmer_size), quasiDico(quasiDico), threshold(threshold), windows_size(windows_size), context_size(context_size), max_edit_distance(max_edit_distance){
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
		// itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}
    
	FunctorQuerySpanKmers () {
	}
    
    void operator() (Sequence& seq){
        const string string_seq = seq.toString();

        int used_windows_size=windows_size;
        if (windows_size==0){ // if windows size == 0 then we use the full read length as windows
            used_windows_size=seq.getDataSize();
        }
		if(not valid_sequence(seq, kmer_size)){return;}
		bool exists;
		associated_read_ids_and_contexts={}; // list of the ids of reads from the bank where a kmer occurs
 		similar_read_ids_position={}; // tmp list of couples <last used position, kmer spanning>
		itKmer->setData (seq.getData());
        
//        if(repeated_kmers(model, *itKmer)){return;}
        u_int position_on_read=-1; // position on the read
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
            
            position_on_read++;
            if (position_on_read<context_size) continue; // avoid first positions with no context
            if (position_on_read>seq.getDataSize()-context_size-kmer_size) break; // avoid last positions with no context
            
            
			quasiDico->get_value((*itKmer)->value().getVal(),exists,associated_read_ids_and_contexts);
            if(!exists) {continue;}

            
			for(auto &read_id_and_contexts: associated_read_ids_and_contexts){
                const u_int32_t read_id = read_id_and_contexts.first;
                const string indexed_context_left = read_id_and_contexts.second.substr(0,context_size); //prefix of size context_size
                const string indexed_context_right = read_id_and_contexts.second.substr(context_size); // suffix of size context_size (starting position context_size)
                
                const string this_context_left = string_seq.substr (position_on_read-context_size,context_size);
                const string this_context_right = string_seq.substr (position_on_read+kmer_size,context_size);
                
                size_t current_edit_distance=edit_distance(indexed_context_left, this_context_left);
                if (current_edit_distance>max_edit_distance) continue;
                    
                current_edit_distance+=edit_distance(indexed_context_right, this_context_right);
                if(current_edit_distance>max_edit_distance) continue;
                
                
                
				std::unordered_map<u_int32_t, vector<bool>>::const_iterator element = similar_read_ids_position.find(read_id);
                vector<bool> position_shared;
				if(element == similar_read_ids_position.end()) {// not inserted yet create an empty vector
                    position_shared = vector<bool>(seq.getDataSize());
                    for (int pos=0;pos<seq.getDataSize();pos++) position_shared[pos]=false;
                    
                }
                else{
                    position_shared = element->second;
                }
                
                
//                cout<<"hohohoh "<<i<<" "<<seq.getDataSize()<<endl; //DEB
                for (int pos=position_on_read-context_size;pos<position_on_read+kmer_size+context_size;pos++){
                    if (pos>seq.getDataSize()) {
//                        cout<<"ok !!"<<endl; //DEB
                        break;
                    }
                    position_shared[pos]=true; 
                }
                
                
                
                similar_read_ids_position[read_id] = position_shared;
                
            }
		}
        
        
		string toPrint;
		bool read_id_printed=false; // Print (and sync file) only if the read is similar to something.
		for (auto &matched_read:similar_read_ids_position){
            const int mpw = max_populated_window(matched_read.second,used_windows_size);
            const float percentage_span_kmer = 100*mpw/float(used_windows_size);
            
			if (percentage_span_kmer >= threshold) {
				if (not read_id_printed){
					read_id_printed=true;
//					synchro->lock();
					toPrint=to_string(seq.getIndex()+1)+":";
				}
				toPrint+=to_string(matched_read.first)+"-"+to_string(mpw)+"-"+to_string(float(percentage_span_kmer))+" ";
			}
            
		}
		if(read_id_printed){
            synchro->lock();
            toPrint+="\n";
            fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
//			fwrite("\n", sizeof(char), 1, outFile);
			synchro->unlock ();
		}
	}
private:
    int max_populated_window(const vector<bool> populated, const int windows_size){
        const int size_vector = populated.size();
        const int last_excluded_starting_window_position = size_vector-windows_size+1;
        int max = 0;
        int number_populated=0;
        
        for(int i=0;i<windows_size;i++){
            if (populated[i]) {
                number_populated++;
            }
        }
        max=number_populated;

        for (int pos=1;pos<last_excluded_starting_window_position;pos+=1){
            if (populated[pos-1]){number_populated--;}
            if (populated[pos+windows_size-1]){number_populated++;}
            if (number_populated>max) max=number_populated;
        }
//        cout<<max<<" "<<windows_size<<endl; //DEBUG
        return max;
    }

    
    int max_start_or_end(const vector<bool> populated, const int windows_size){
        const int start = max_starting_populated_window(populated,windows_size);
        const int end = max_ending_populated_window(populated,windows_size);
        
        if (start>end) return start;
        return end;
        
    }
    
    int max_starting_populated_window(const vector<bool> populated, const int windows_size){
        int number_populated=0;
        
        for(int i=0;i<windows_size;i++){
            if (populated[i]) {
                number_populated++;
            }
        }
        return number_populated;
    }
    
    
    
    int max_ending_populated_window(const vector<bool> populated, const int windows_size){
        const int size_vector = populated.size();
        const int last_excluded_starting_window_position = size_vector-windows_size+1;
        int number_populated=0;
        
        for(int i=size_vector-windows_size;i<size_vector;i++){
            if (populated[i]) {
                number_populated++;
            }
        }
        return number_populated;
    }

};



void SRC_linker_fuzzy::parse_query_sequences (int threshold, const int nbCores, const int windows_size){
    BankAlbum banks (getInput()->getStr(STR_URI_QUERY_INPUT));
    const std::vector<IBank*>& banks_of_queries = banks.getBanks();
    const int number_of_read_sets = banks_of_queries.size();
    
    
	cout<<"Query "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_QUERY_INPUT)<<endl;
	FILE * outFile;
	outFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
    string message("#query_read_id [target_read_id-kmer_span (k="+to_string(kmer_size)+")-kmer_span query percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)\n"+"#Target read set: "+getInput()->getStr(STR_URI_BANK_INPUT)+"\n");
    fwrite((message).c_str(), sizeof(char), message.size(), outFile);
    
    
    for( int bank_id=0;bank_id<number_of_read_sets;bank_id++){ // iterate each bank
        
        
        IBank* bank=banks_of_queries[bank_id];
        LOCAL (bank);
        string message("#Query read set number "+bank->getId()+"\n");
        fwrite((message).c_str(), sizeof(char), message.size(), outFile);
        string progressMessage("Querying read set "+bank->getId());
        ProgressIterator<Sequence> itSeq (*bank, progressMessage.c_str());
        ISynchronizer* synchro = System::thread().newSynchronizer();
        Dispatcher dispatcher (nbCores, 1000);
        dispatcher.iterate (itSeq, FunctorQuerySpanKmers(synchro,outFile, kmer_size,&quasiDico, threshold, windows_size,context_size,max_edit_distance));
        delete synchro;
    }
	fclose (outFile);
    
    
}


void SRC_linker_fuzzy::execute (){
	int nbCores = getInput()->getInt(STR_CORE);
	int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
    max_edit_distance = getInput()->getInt(STR_MAX_EDIT_DISTANCE);
    gamma_value = getInput()->getInt(STR_GAMMA);
    context_size = getInput()->getInt(STR_CONTEXT_SIZE);
	// IMPORTANT NOTE:
	// Actually, during the filling of the dictionary values, one may fall on non solid non indexed kmers
	// that are quasi dictionary false positives (ven with a non null fingerprint. This means that one nevers knows in advance how much
	// values are gonna be stored for all kmers. This is why I currently us a vector<u_int32_t> for storing read ids associated to a kmer.

	// We need a non null finger print because of non solid non indexed kmers
	//	if (getInput()->getStr(STR_URI_BANK_INPUT).compare(getInput()->getStr(STR_URI_QUERY_INPUT))==0)
	//		fingerprint_size=0;
	cout<<"fingerprint = "<<fingerprint_size<<endl;
	create_quasi_dictionary(fingerprint_size, nbCores);
	fill_quasi_dictionary(nbCores);

	int threshold = getInput()->getInt(STR_THRESHOLD);
    int windows_size = getInput()->getInt(STR_WINDOWS_SIZE);
	parse_query_sequences(threshold, nbCores, windows_size);

	getInfo()->add (1, &LibraryInfo::getInfo());
	getInfo()->add (1, "input");
	getInfo()->add (2, "Reference bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
    getInfo()->add (2, "Query bank",  "%s",  getInput()->getStr(STR_URI_QUERY_INPUT).c_str());
    getInfo()->add (2, "windows_size size",  "%d",  windows_size);
    getInfo()->add (2, "Kmer size",  "%d",  kmer_size);
    getInfo()->add (2, "Context size",  "%d",  context_size);
	getInfo()->add (2, "Fingerprint size",  "%d",  fingerprint_size);
    getInfo()->add (2, "gamma",  "%d",  gamma_value);
	getInfo()->add (2, "Minimal kmer span percentage",  "%d",  threshold);
	getInfo()->add (1, "output");
	getInfo()->add (2, "Results written in",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());
}













