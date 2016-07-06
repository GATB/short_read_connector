#include <SRC_linker_ram.hpp>


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


SRC_linker_ram::SRC_linker_ram ()  : Tool ("SRC_linker_ram"){
	// We add some custom arguments for command line interface
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
	getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT, "query input",    true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal percentage of shared kmer span for considering 2 reads as similar. The kmer span is the number of bases from the read query covered by a kmer shared with the target read. If a read of length 80 has a kmer-span of 60 with another read from the bank (of unkonwn size), then the percentage of shared kmer span is 75%.",    false, "75"));
	getParser()->push_back (new OptionOneParam (STR_GAMMA, "gamma value",    false, "2"));
	getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "8"));
	getParser()->push_back (new OptionOneParam (STR_CORE, "Number of thread",    false, "1"));
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
	int kmer_size;

	FunctorIndexer(quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >& quasiDico, int kmer_size)  :  quasiDico(quasiDico), kmer_size(kmer_size) {
	}

	void operator() (Sequence& seq){
		if(not valid_sequence(seq,kmer_size)){return;}
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);
		itKmer.setData (seq.getData());
        
//        if(repeated_kmers(model, itKmer)){return;}
		u_int32_t read_id = static_cast<u_int32_t>(seq.getIndex()+1);
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
	dispatcher.iterate (itSeq, FunctorIndexer(quasiDico, kmer_size));
}



class FunctorQuerySpanKmers // FunctorQuery used after claires discussion: number of positions covered by a shared kmer.
{
public:
	ISynchronizer* synchro;
	FILE* outFile;
	int kmer_size;
	quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t>* quasiDico;
	int threshold;
	vector<u_int32_t> associated_read_ids;
	std::unordered_map<u_int32_t, std::pair <u_int,u_int>> similar_read_ids_position_count; // each bank read id --> couple<next viable position (without overlap), number of shared kmers>
	Kmer<KMER_SPAN(1)>::ModelCanonical model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator* itKmer;
    
	FunctorQuerySpanKmers(const FunctorQuerySpanKmers& lol)
	{
		synchro=lol.synchro;
		outFile=lol.outFile;
		kmer_size=lol.kmer_size;
		quasiDico=lol.quasiDico;
		threshold=lol.threshold;
		associated_read_ids=lol.associated_read_ids;
		similar_read_ids_position_count=lol.similar_read_ids_position_count;
		model=lol.model;
		itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}
    
	FunctorQuerySpanKmers (ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  quasidictionaryVectorKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >* quasiDico, const int threshold)
	: synchro(synchro), outFile(outFile), kmer_size(kmer_size), quasiDico(quasiDico), threshold(threshold) {
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
		// itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}
    
	FunctorQuerySpanKmers () {
	}
    
	void operator() (Sequence& seq){
		if(not valid_sequence(seq, kmer_size)){return;}
		bool exists;
		associated_read_ids={}; // list of the ids of reads from the bank where a kmer occurs
 		similar_read_ids_position_count={}; // tmp list of couples <last used position, kmer spanning>
		itKmer->setData (seq.getData());
		
//        if(repeated_kmers(model, *itKmer)){return;}
        u_int i=0; // position on the read
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
			quasiDico->get_value((*itKmer)->value().getVal(),exists,associated_read_ids);
			if(!exists) {++i;continue;}
			for(auto &read_id: associated_read_ids){
				std::unordered_map<u_int32_t, std::pair <u_int,u_int>>::const_iterator element = similar_read_ids_position_count.find(read_id);
				if(element == similar_read_ids_position_count.end()) {// not inserted yet:
					similar_read_ids_position_count[read_id]=std::make_pair(i, kmer_size);
				}else{  // a kmer is already shared with this read
					std::pair <int,int> lastpos_spankmer = (element->second);
                    // update spanning, up to a kmer size
                    if ((i-lastpos_spankmer.first)<kmer_size)   lastpos_spankmer.second += i-lastpos_spankmer.first;
                    else                                        lastpos_spankmer.second += kmer_size;
                    lastpos_spankmer.first=i;                                            // update last position of a shared kmer with this read
                    similar_read_ids_position_count[read_id] = lastpos_spankmer;
				}
			}
			++i;
		}
		string toPrint;
		bool read_id_printed=false; // Print (and sync file) only if the read is similar to something.
		for (auto &matched_read:similar_read_ids_position_count){
            float percentage_span_kmer = 100*std::get<1>(matched_read.second)/float(seq.getDataSize());
			if (percentage_span_kmer >= threshold) {
				if (not read_id_printed){
					read_id_printed=true;
//					synchro->lock();
					toPrint=to_string(seq.getIndex()+1)+":";
//					fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
				}
				toPrint+=to_string(matched_read.first)+"-"+to_string(std::get<1>(matched_read.second))+"-"+to_string(float(percentage_span_kmer))+" ";
//				fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
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
};



void SRC_linker_ram::parse_query_sequences (int threshold, const int nbCores){
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
        dispatcher.iterate (itSeq, FunctorQuerySpanKmers(synchro,outFile, kmer_size,&quasiDico, threshold));
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
	int nbCores = getInput()->getInt(STR_CORE);
	int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
    gamma_value = getInput()->getInt(STR_GAMMA);
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
	parse_query_sequences(threshold, nbCores);

	getInfo()->add (1, &LibraryInfo::getInfo());
	getInfo()->add (1, "input");
	getInfo()->add (2, "Reference bank",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
	getInfo()->add (2, "Query bank",  "%s",  getInput()->getStr(STR_URI_QUERY_INPUT).c_str());
    getInfo()->add (2, "Kmer size",  "%d",  kmer_size);
	getInfo()->add (2, "Fingerprint size",  "%d",  fingerprint_size);
    getInfo()->add (2, "gamma",  "%d",  gamma_value);
	getInfo()->add (2, "Minimal kmer span percentage",  "%d",  threshold);
	getInfo()->add (1, "output");
	getInfo()->add (2, "Results written in",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());
}
