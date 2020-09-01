// We include what we need for the test
#include <gatb/gatb_core.hpp>
#include <iostream>
#include "boolean_vector.hpp"
// We use the required packages
using namespace std;
/********************************************************************************/
/*                         Bank conversion with some filtering                  */
/*                                                                              */
/* This snippet shows how to copy a bank to another one with some filtering     */
/* criteria. The criteria is described through a functor.                       */
/* Here, we can keep only every 'modulo' sequence that have no 'N' in data.     */
/*                                                                              */
/* Cmd-line: bank9 <fasta/q file> <outfile> <modulo>                            */
/*                                                                              */
/* Sample: bank9 gatb-core/gatb-core/test/db/reads1.fa reads1_filtered.fa 5     */
/*                                                                              */
/********************************************************************************/
// We define a functor for our sequence filtering.
struct FilterFunctor
{
    BooleanVector bv;
    const std::string & bv_file_name;
    FilterFunctor (const std::string & bv_file_name) : bv_file_name(bv_file_name) {
		// Read the boolean vector in bv file
		bv.read(bv_file_name);
    }
    bool operator ()  (Sequence& seq)
    {
        //DEBUG
    //     cout<<"seq "<<seq.getComment()<<" "<<seq.getIndex()<<endl;
    //    if (bv.is_set(seq.getIndex())){
    //        cout<<"OK "<<seq.getIndex()<<endl;
    //    }
    //    else{
    //        cout<<"KO "<<seq.getIndex()<<endl;
    //    }
        //END DEBUG
        return (bv.is_set(seq.getIndex()));
    }
};
/********************************************************************************/
int main (int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "you must provide an input bank name, a boolean vector and an output bank name:" << endl;
        cerr << "   1) input read file URI"  << endl;
        cerr << "   2) input boolean vector (.bv) file" << endl;
        cerr << "   3) output URI" << endl;
        return EXIT_FAILURE;
    }
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        // We declare an input Bank and use it locally
        IBank* inputBank = Bank::open (argv[1]);
        LOCAL (inputBank);
        // We declare an output Bank
        BankFasta outputBank (argv[3]);

        
        // We create a filter on the input bank.
        // Note the first argument of the FilterIterator is an iterator of the input bank and the
        // second argument is the filtering functor.
        FilterIterator<Sequence,FilterFunctor>* itFilter = new FilterIterator<Sequence,FilterFunctor>(
            inputBank->iterator(),
            FilterFunctor (argv[2])
        );
            
        // We create an iterator that will provide progress notification.
        // Note the estimated number of items takes into account the modulo
        ProgressIterator<Sequence> itSeq (itFilter, "Filtering bank", inputBank->estimateNbItems() );
        // We loop over sequences.
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            // We insert the current sequence into the output bank.
            outputBank.insert (*itSeq);
        }
        // We make sure that the output bank is flushed correctly.
        outputBank.flush ();
    }
    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
    }
}