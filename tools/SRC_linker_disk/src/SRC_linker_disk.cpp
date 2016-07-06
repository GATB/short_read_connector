#include <SRC_linker_disk.hpp>


using namespace std;


/********************************************************************************/


static const char* STR_URI_BANK_INPUT = "-bank";
static const char* STR_URI_QUERY_INPUT = "-query";
static const char* STR_FINGERPRINT = "-fingerprint_size";
static const char* STR_GAMMA = "-gamma";
static const char* STR_THRESHOLD = "-kmer_threshold";
static const char* STR_OUT_FILE = "-out";
static const char* STR_CORE = "-core";


static const uint buff(5);

SRC_linker_disk::SRC_linker_disk ()  : Tool ("SRC_linker_disk"){
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
	getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT, "query input",    true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal percentage of shared kmer span for considering 2 reads as similar. The kmer span is the number of bases from the read query covered by a kmer shared with the target read. If a read of length 80 has a kmer-span of 60 with another read from the bank (of unkonwn size), then the percentage of shared kmer span is 75%.",    false, "75"));
	getParser()->push_back (new OptionOneParam (STR_GAMMA, "gamma value",    false, "2"));
	getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "8"));
	getParser()->push_back (new OptionOneParam (STR_CORE, "Number of thread",    false, "1"));
}


struct FunctorCount{
	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t > &quasiDico;
	int kmer_size;

	FunctorCount(quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t >& quasiDico, int kmer_size)  :  quasiDico(quasiDico), kmer_size(kmer_size) {}

	void operator() (Sequence& seq){
		if(not valid_sequence(seq, kmer_size)){return;}
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);//TODO we can do better than create this each time
        
		itKmer.setData (seq.getData());
        
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			uint64_t count;
			bool exists;
			quasiDico.get_value((itKmer)->value().getVal(),exists,count);
			if(exists){
				bool lol(quasiDico.set_value((itKmer)->value().getVal(), count+1));
			}
		}
	}
};


struct FunctorWriter{
	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t > &quasiDico;
	int kmer_size;
	uint64_t position;
	FILE * pFile;
	ISynchronizer* synchro;

	FunctorWriter(quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t >& quasiDico, int kmer_size, FILE * pFile,ISynchronizer* synchro)
	:  quasiDico(quasiDico), kmer_size(kmer_size), pFile(pFile),synchro(synchro) {
		position=0;
	}

	void operator() (Kmer<>::Count & it){
		bool exists;
		uint64_t abundance;
		quasiDico.get_value(it.value.getVal(),exists,abundance);
		quasiDico.set_value(it.value.getVal(), position);
		string toPrint(8*(abundance+1),0);
		// cout<<abundance<<endl;
		// cin.get();
		fwrite(toPrint.c_str(), 1, 8*(abundance+1), pFile);
		// cout<<lol<<endl;
		// cin.get();
		position+=8*(abundance+1);
		// fflush(pFile);
		// cout<<position<<" "<<ftell(pFile)<<endl;
	}
};


void SRC_linker_disk::create_quasi_dictionary (int fingerprint_size, int nbCores){
	const int display = getInput()->getInt (STR_VERBOSE);
	auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (getInput()->getStr(STR_URI_GRAPH)));
	Group& dskGroup = storage->getGroup("dsk");
	kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");
	nbSolidKmers = solidKmers.getNbItems();
	if(nbSolidKmers==0){cout<<"No solid kmers in bank -- exit"<<endl;exit(0);}
	//we compute the quasidico
	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());
	quasiDico = quasidictionaryKeyGeneric<IteratorKmerH5Wrapper, uint64_t> (nbSolidKmers, iteratorOnKmers, fingerprint_size, gamma_value);
	//we count the occurence of kmer in the bank file (including false positive)
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_BANK_INPUT));
	ISynchronizer* synchro = System::thread().newSynchronizer();
	ProgressIterator<Sequence> itSeq (*bank);
	Dispatcher dispatcher (nbCores, 1000);
	dispatcher.iterate (itSeq, FunctorCount(quasiDico, kmer_size));
	//we create the blank file
	pFile = fopen ("Erase_Me", "w+");
	ProgressIterator<Kmer<>::Count> itKmers (solidKmers.iterator(), "Indexing solid kmers", nbSolidKmers);
	Dispatcher dispatcher2 (1, 1000);
	dispatcher2.iterate (itKmers, FunctorWriter(quasiDico, kmer_size,pFile,synchro));
	string toPrint((buff+1)*8,0);
	fwrite(toPrint.c_str(),1,(buff+1)*8, pFile);
	fflush(pFile);
	fclose(pFile);
}


struct FunctorIndexer{
	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t > *quasiDico;
	int kmer_size;
	FILE * pFile;
	ISynchronizer* synchro;

	FunctorIndexer(quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t >* quasiDico, int kmer_size,ISynchronizer* synchro)
	:  quasiDico(quasiDico), kmer_size(kmer_size),synchro(synchro) {
	}

	FunctorIndexer(const FunctorIndexer& lol)
	{
		quasiDico=lol.quasiDico;
		kmer_size=lol.kmer_size;
		pFile = fopen ("Erase_Me", "r+");
		synchro=lol.synchro;
	}

	void operator() (Sequence& seq){
		if(not valid_sequence(seq, kmer_size)){return;}
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);
		itKmer.setData (seq.getData());
		uint64_t read_id(static_cast<uint64_t>(seq.getIndex())+1);
		uint64_t position;
		uint64_t id[buff];
		bool exists;
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			// find the position in the file on write the ReadID there
			// find the initial pos with the quasi-dico and find the first unused (non 0)
			quasiDico->get_value(itKmer->value().getVal(),exists,position);
			if(not exists) {continue;}
			fseek (pFile , position , SEEK_SET);
			while(true){
				fread(id,1,8*buff,pFile);
				// cout<<"lol"<<endl;
				for(uint ii(0);ii<buff;++ii){
					if(id[ii]==0){
						fseek ( pFile , position+8*ii , SEEK_SET );
						fwrite(&read_id, 1, 8,  pFile);
						goto end;
					}
				}
				position+=8*buff;
			}
			end:
			;
			// cout<<"lol"<<endl;
			// cin.get();
			// synchro->unlock();
		}
	}
};


void SRC_linker_disk::fill_quasi_dictionary (const int nbCores){
	bool exists;
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_BANK_INPUT));
	cout<<"Index "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	ISynchronizer* synchro(System::thread().newSynchronizer());
	Dispatcher dispatcher (1, 1000);
	dispatcher.iterate (itSeq, FunctorIndexer(&quasiDico, kmer_size,synchro));
}


class FunctorQuery{
public:
	ISynchronizer* synchro;
	FILE* outFile;
	FILE* pFile;
	int kmer_size;
	quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t>* quasiDico;
	int threshold;
	uint64_t position;
	std::vector<uint64_t> associated_read_ids;
	std::unordered_map<uint64_t, std::pair <uint,uint>> similar_read_ids_position_count; // each bank read id --> couple<next viable position (without overlap), number of shared kmers>
	Kmer<KMER_SPAN(1)>::ModelCanonical model;
	Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator* itKmer;

	FunctorQuery(const FunctorQuery& lol)
	{
		synchro=lol.synchro;
		outFile=lol.outFile;
		kmer_size=lol.kmer_size;
		quasiDico=lol.quasiDico;
		threshold=lol.threshold;
		associated_read_ids=lol.associated_read_ids;
		similar_read_ids_position_count=lol.similar_read_ids_position_count;
		model=lol.model;
		pFile = fopen ("Erase_Me", "r");
		itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}


	FunctorQuery (ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  quasidictionaryKeyGeneric <IteratorKmerH5Wrapper, uint64_t >* quasiDico, const int threshold)
	: synchro(synchro), outFile(outFile), kmer_size(kmer_size) , quasiDico(quasiDico), threshold(threshold) {
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
	}


	~FunctorQuery () {}


	void operator() (Sequence& seq){
		if(not valid_sequence(seq, kmer_size)){return;}
		bool exists;
		associated_read_ids={};
 		similar_read_ids_position_count={};
		itKmer->setData (seq.getData());
		uint i(0);
		uint64_t id[buff]; // position on the read
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
			associated_read_ids={};
			quasiDico->get_value((*itKmer)->value().getVal(),exists,position);
			if(not exists) {++i;continue;}
			fseek ( pFile , position , SEEK_SET );
			while(true){
				fread (id,1,8*buff,pFile);
				for(uint ii(0);ii<buff;++ii){
					if(id[ii]!=0){
						associated_read_ids.push_back(id[ii]);
					}else{
						goto end;
					}
				}
			}
			end:
//			for(auto &read_id: associated_read_ids){
//				std::unordered_map<uint64_t, std::pair <uint,uint>>::const_iterator element = similar_read_ids_position_count.find(read_id);
//				if(element == similar_read_ids_position_count.end()) {// not inserted yet:
//					similar_read_ids_position_count[read_id]=std::make_pair(i+kmer_size, 1);
//				}else{  // a kmer is already shared with this read
//					std::pair <int,int> viablepos_nbshared = (element->second);
//					if(i>=viablepos_nbshared.first){ // the current position does not overlap the previous shared kmer
//						viablepos_nbshared.first = i+kmer_size; // next non overlapping position
//						viablepos_nbshared.second = viablepos_nbshared.second+1; // a new kmer shared.
//						similar_read_ids_position_count[read_id] = viablepos_nbshared;
//					}
//				}
//			}
            
            for(auto &read_id: associated_read_ids){
				std::unordered_map<uint64_t, std::pair <uint,uint>>::const_iterator element = similar_read_ids_position_count.find(read_id);
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
					toPrint=to_string(seq.getIndex()+1)+":";
				}
				toPrint+=to_string(matched_read.first)+"-"+to_string(std::get<1>(matched_read.second))+"-"+to_string(float(percentage_span_kmer))+" ";
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


void SRC_linker_disk::parse_query_sequences (int threshold, const int nbCores){
    
    BankAlbum banks (getInput()->getStr(STR_URI_QUERY_INPUT));
    const std::vector<IBank*>& banks_of_queries = banks.getBanks();
    const int number_of_read_sets = banks_of_queries.size();
    
    
	cout<<"Query "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_QUERY_INPUT)<<endl;
	FILE * outFile;
	outFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
    string message("#query_read_id [target_read_id-kmer_span (k="+to_string(kmer_size)+")-kmer_span query percentage]* or U (unvalid read, containing not only ACGT characters or low complexity read)\n");
    fwrite((message).c_str(), sizeof(char), message.size(), outFile);
    
    for( int bank_id=0;bank_id<number_of_read_sets;bank_id++){ // iterate each bank
        
        string message("#Query read set number "+to_string(bank_id)+"\n");
        fwrite((message).c_str(), sizeof(char), message.size(), outFile);
        
        IBank* bank=banks_of_queries[bank_id];
        LOCAL (bank);
        string progressMessage("Querying read set number "+to_string(bank_id));
        ProgressIterator<Sequence> itSeq (*bank, progressMessage.c_str());
        ISynchronizer* synchro = System::thread().newSynchronizer();
        Dispatcher dispatcher (nbCores, 1000);
        dispatcher.iterate (itSeq, FunctorQuery(synchro,outFile, kmer_size,&quasiDico, threshold));
        delete synchro;
    }
	fclose (pFile);
}


void SRC_linker_disk::execute (){
	int nbCores = getInput()->getInt(STR_CORE);
	int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
    gamma_value = getInput()->getInt(STR_GAMMA);
	cout<<1<<endl;
	create_quasi_dictionary(fingerprint_size,nbCores);
	cout<<2<<endl;
	fill_quasi_dictionary(nbCores);
	cout<<3<<endl;
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
