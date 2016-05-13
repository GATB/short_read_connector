#include <commet_linked_disk.hpp>


using namespace std;


/********************************************************************************/


static const char* STR_URI_BANK_INPUT = "-bank";
static const char* STR_URI_QUERY_INPUT = "-query";
static const char* STR_FINGERPRINT = "-fingerprint_size";
static const char* STR_THRESHOLD = "-kmer_threshold";
static const char* STR_OUT_FILE = "-out";
static const char* STR_CORE = "-core";


commet_linked_disk::commet_linked_disk ()  : Tool ("commet_linked_disk"){
	getParser()->push_back (new OptionOneParam (STR_URI_GRAPH, "graph input",   true));
	getParser()->push_back (new OptionOneParam (STR_URI_BANK_INPUT, "bank input",    true));
	getParser()->push_back (new OptionOneParam (STR_URI_QUERY_INPUT, "query input",    true));
	getParser()->push_back (new OptionOneParam (STR_OUT_FILE, "output_file",    true));
	getParser()->push_back (new OptionOneParam (STR_THRESHOLD, "Minimal number of shared kmers for considering 2 reads as similar",	false, "10"));
	getParser()->push_back (new OptionOneParam (STR_FINGERPRINT, "fingerprint size",    false, "8"));
	getParser()->push_back (new OptionOneParam (STR_CORE, "Number of thread",    false, "1"));
}


static int NT2int(char nt){return (nt>>1)&3;}


bool correct(Sequence& seq){
	const char* data = seq.getDataBuffer();
	int DUSTSCORE[64]={0}; // all tri-nucleotides

	size_t lenseq =seq.getDataSize();
	if (data[0]!='A' && data[0]!='C' && data[0]!='G' && data[0]!='T')  { return false; }
	if (data[1]!='A' && data[1]!='C' && data[1]!='G' && data[1]!='T')  { return false; }

	for (int j=2; j<lenseq; ++j){
		++DUSTSCORE[NT2int(data[j-2])*16 + NT2int(data[j-1])*4 + NT2int(data[j])];
		if (data[j]!='A' && data[j]!='C' && data[j]!='G' && data[j]!='T')  { return false; }
	}
	int m,s=0;

	for (int i=0; i<64; ++i)
	{
		m = DUSTSCORE[i];
		s  += (m*(m-1))/2;
	}

	return s<((lenseq-2)/4 * (lenseq-6)/4)/2;
}


struct FunctorCount{
	quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t > &quasiDico;
	int kmer_size;

	FunctorCount(quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >& quasiDico, int kmer_size)  :  quasiDico(quasiDico), kmer_size(kmer_size) {}

	void operator() (Sequence& seq){
		if(not correct(seq)){return;}
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);//TODO we can do better than create this each time
		itKmer.setData (seq.getData());
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			uint32_t count;
			bool exists;
			quasiDico.get_value((itKmer)->value().getVal(),exists,count);
			if(exists){
				bool lol(quasiDico.set_value((itKmer)->value().getVal(), count+1));
			}
		}
	}
};


struct FunctorWriter{
	quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t > &quasiDico;
	int kmer_size;
	uint32_t position;
	FILE * pFile;
	ISynchronizer* synchro;

	FunctorWriter(quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >& quasiDico, int kmer_size, FILE * pFile,ISynchronizer* synchro)
	:  quasiDico(quasiDico), kmer_size(kmer_size), pFile(pFile),synchro(synchro) {
		position=0;
	}

	void operator() (Kmer<>::Count & it){
		bool exists;
		uint32_t abundance;
		quasiDico.get_value(it.value.getVal(),exists,abundance);
		quasiDico.set_value(it.value.getVal(), position);
		string toPrint(4*(abundance+1),0);
		fwrite(toPrint.c_str(), 1, 4*(abundance+1), pFile);
		position+=4*(abundance+1);
	}
};


void commet_linked_disk::create_quasi_dictionary (int fingerprint_size, int nbCores){
	const int display = getInput()->getInt (STR_VERBOSE);
	auto_ptr<Storage> storage (StorageFactory(STORAGE_HDF5).load (getInput()->getStr(STR_URI_GRAPH)));
	Group& dskGroup = storage->getGroup("dsk");
	kmer_size = atoi(dskGroup.getProperty("kmer_size").c_str());
	Partition<Kmer<>::Count>& solidKmers = dskGroup.getPartition<Kmer<>::Count> ("solid");
	nbSolidKmers = solidKmers.getNbItems();
	if(nbSolidKmers==0){cout<<"No solid kmers in bank -- exit"<<endl;exit(0);}
	//we compute the quasidico
	IteratorKmerH5Wrapper iteratorOnKmers (solidKmers.iterator());
	int gamma(10);//TODO parameter gamma
	quasiDico = quasiDictionnaryKeyGeneric<IteratorKmerH5Wrapper, uint32_t> (nbSolidKmers, iteratorOnKmers, fingerprint_size, gamma);
	//we count the occurence of kmer in the bank file (including false positive)
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_BANK_INPUT));
	ISynchronizer* synchro = System::thread().newSynchronizer();
	ProgressIterator<Sequence> itSeq (*bank);
	Dispatcher dispatcher (nbCores, 1000);
	dispatcher.iterate (itSeq, FunctorCount(quasiDico, kmer_size));
	//we create the blank file
	pFile = fopen ("Erase_Me", "w+");
	uint position(0);
	uint i(0);
	uint abundance(0);
	bool exists;
	ProgressIterator<Kmer<>::Count> itKmers (solidKmers.iterator(), "Indexing solid kmers", nbSolidKmers);
	Dispatcher dispatcher2 (1, 10000);
	dispatcher2.iterate (itKmers, FunctorWriter(quasiDico, kmer_size,pFile,synchro));
	string toPrint(4,0);
	fwrite(toPrint.c_str(), 1, 4, pFile);
}


struct FunctorIndexer{
	quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t > &quasiDico;
	int kmer_size;
	FILE * pFile;
	ISynchronizer* synchro;

	FunctorIndexer(quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >& quasiDico, int kmer_size, FILE * pFile,ISynchronizer* synchro)
	:  quasiDico(quasiDico), kmer_size(kmer_size), pFile(pFile),synchro(synchro) {
	}

	void operator() (Sequence& seq){
		if(not correct(seq)){return;}
		Kmer<KMER_SPAN(1)>::ModelCanonical model (kmer_size);
		Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator itKmer (model);//TODO we can do better than create this each time
		itKmer.setData (seq.getData());
		uint32_t read_id(static_cast<u_int32_t>(seq.getIndex())+1),position;
		uint32_t id(0);
		bool exists;
		for (itKmer.first(); !itKmer.isDone(); itKmer.next()){
			// find the position in the file on write the ReadID there
			// find the initial pos with the quasi-dico and find the first unused (non 0)
			quasiDico.get_value(itKmer->value().getVal(),exists,position);
			if(not exists) {continue;}
			synchro->lock();
			fseek (pFile , position , SEEK_SET);
			while(true){
				int lol=fread(&id,1,4,pFile); //we can read much and then read on an array
				if(id==0){
					fseek ( pFile , -4 , SEEK_CUR );
					fwrite(&read_id, 1, 4,  pFile);
					break;
				}
			}
			synchro->unlock();
		}
	}
};


void commet_linked_disk::fill_quasi_dictionary (const int nbCores){
	bool exists;
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_BANK_INPUT));
	cout<<"Index "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_BANK_INPUT)<<endl;
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	ISynchronizer* synchro(System::thread().newSynchronizer());
	Dispatcher dispatcher (nbCores, 10000);
	dispatcher.iterate (itSeq, FunctorIndexer(quasiDico, kmer_size,pFile,synchro));
}


class FunctorQuery{
public:
	ISynchronizer* synchro;
	FILE* outFile;
	FILE* pFile;
	int kmer_size;
	quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t>* quasiDico;
	int threshold;
	uint32_t position;
	int count;
	std::vector<uint32_t> associated_read_ids;
	std::unordered_map<u_int32_t, std::pair <u_int,u_int>> similar_read_ids_position_count; // each bank read id --> couple<next viable position (without overlap), number of shared kmers>
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
		pFile=lol.pFile;
		itKmer = new Kmer<KMER_SPAN(1)>::ModelCanonical::Iterator (model);
	}


	FunctorQuery (ISynchronizer* synchro, FILE* outFile,  const int kmer_size,  quasiDictionnaryKeyGeneric <IteratorKmerH5Wrapper, u_int32_t >* quasiDico, const int threshold,  FILE * pFile)
	: synchro(synchro), outFile(outFile), kmer_size(kmer_size), pFile(pFile), quasiDico(quasiDico), threshold(threshold) {
		model=Kmer<KMER_SPAN(1)>::ModelCanonical (kmer_size);
	}


	~FunctorQuery () {}


	void operator() (Sequence& seq){
		if(not correct(seq)){return;}
		bool exists;
		associated_read_ids={};
 		similar_read_ids_position_count={};
		itKmer->setData (seq.getData());
		uint32_t i(0);
		uint32_t id(0); // position on the read
		for (itKmer->first(); !itKmer->isDone(); itKmer->next()){
			associated_read_ids={};
			quasiDico->get_value((*itKmer)->value().getVal(),exists,position);
			if(not exists) {++i;continue;}
			fseek ( pFile , position , SEEK_SET );
			while(true){
				fread (&id,1,4,pFile);//TODO we can read much and then look at an array
				if(id!=0){
					associated_read_ids.push_back(id);
				}else{
					break;
				}
			}
			for(auto &read_id: associated_read_ids){
				std::unordered_map<u_int32_t, std::pair <u_int,u_int>>::const_iterator element = similar_read_ids_position_count.find(read_id);
				if(element == similar_read_ids_position_count.end()) {// not inserted yet:
					similar_read_ids_position_count[read_id]=std::make_pair(i+kmer_size, 1);
				}else{  // a kmer is already shared with this read
					std::pair <int,int> viablepos_nbshared = (element->second);
					if(i>=viablepos_nbshared.first){ // the current position does not overlap the previous shared kmer
						viablepos_nbshared.first = i+kmer_size; // next non overlapping position
						viablepos_nbshared.second = viablepos_nbshared.second+1; // a new kmer shared.
						similar_read_ids_position_count[read_id] = viablepos_nbshared;
					}
				}
			}
			++i;
		}
		string toPrint;
		bool read_id_printed=false; // Print (and sync file) only if the read is similar to something.
		for (auto &matched_read:similar_read_ids_position_count){
			if (std::get<1>(matched_read.second) > threshold) {
				if (not read_id_printed){
					read_id_printed=true;
					synchro->lock();
					toPrint=to_string(seq.getIndex()+1)+":";
					fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
				}
				toPrint=to_string(matched_read.first+1)+"-"+to_string(std::get<1>(matched_read.second))+" ";
				fwrite(toPrint.c_str(), sizeof(char), toPrint.size(), outFile);
			}
		}
		if(read_id_printed){
			fwrite("\n", sizeof(char), 1, outFile);
			synchro->unlock ();
		}
	}
};


void commet_linked_disk::parse_query_sequences (int threshold, const int nbCores){
	IBank* bank = Bank::open (getInput()->getStr(STR_URI_QUERY_INPUT));
	cout<<"Query "<<kmer_size<<"-mers from bank "<<getInput()->getStr(STR_URI_QUERY_INPUT)<<endl;
	FILE * outFile;
	outFile = fopen (getInput()->getStr(STR_OUT_FILE).c_str(), "wb");
	string message("#query_read_id [target_read_id number_shared_"+to_string(kmer_size)+"mers]* or U (unvalid read, containing not only ACGT characters or low complexity read)\n");
	fwrite((message).c_str(), sizeof(char), message.size(), outFile);
	LOCAL (bank);
	ProgressIterator<Sequence> itSeq (*bank);
	ISynchronizer* synchro = System::thread().newSynchronizer();
	Dispatcher dispatcher (nbCores, 10000);
	dispatcher.iterate (itSeq, FunctorQuery(synchro,outFile, kmer_size,&quasiDico, threshold, pFile));
	fclose (pFile);
	delete synchro;
}


void commet_linked_disk::execute (){
	int nbCores = getInput()->getInt(STR_CORE);
	int fingerprint_size = getInput()->getInt(STR_FINGERPRINT);
	cout<<"fingerprint = "<<fingerprint_size<<endl;
	create_quasi_dictionary(fingerprint_size,nbCores);
	fill_quasi_dictionary(nbCores);
	int threshold = getInput()->getInt(STR_THRESHOLD);
	parse_query_sequences(threshold-1, nbCores); //-1 avoids >=

	getInfo()->add (1, &LibraryInfo::getInfo());
	getInfo()->add (1, "input");
	getInfo()->add (2, "Reference bank:",  "%s",  getInput()->getStr(STR_URI_BANK_INPUT).c_str());
	getInfo()->add (2, "Query bank:",  "%s",  getInput()->getStr(STR_URI_QUERY_INPUT).c_str());
	getInfo()->add (2, "Fingerprint size:",  "%d",  fingerprint_size);
	getInfo()->add (2, "Threshold size:",  "%d",  threshold);
	getInfo()->add (1, "output");
	getInfo()->add (2, "Results written in:",  "%s",  getInput()->getStr(STR_OUT_FILE).c_str());
}
