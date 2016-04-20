/*
 * GzFileReader.h
 *
 *  Created on: 16 mars 2016
 *      Author: gbenoit
 */

#ifndef QUASI_DICTIONNARY_SRC_ITERATORKMERH5_H_
#define QUASI_DICTIONNARY_SRC_ITERATORKMERH5_H_


#include <gatb/gatb_core.hpp>


//template<typename Key>
class IteratorKmerH5 : public std::iterator<std::forward_iterator_tag, const Kmer<>::Count>//, const Key>
{
public:
	IteratorKmerH5()  : iterator(0), pos(0) {}

	IteratorKmerH5(Iterator<Kmer<>::Count>* iterator, const int kmer_size)  : iterator(iterator), pos(0)  {
		iterator->first();
		model = Kmer<>::ModelCanonical(kmer_size);
	}

	unsigned long long  const& operator*()  {
		Kmer<>::Count& count = iterator->item();
		unsigned long long kmer_int_value=oahash(count.value);
		return kmer_int_value;
		//cout << std::get<0>(iterator->item()) << endl;
		//return 0;
		//return _item;
		//KmerLca kmerLca = iterator->item();
		//cout << kmerLca.first << " " << kmerLca.second << endl;
		//return kmerLca;
	}

	IteratorKmerH5& operator++()
	{
        iterator->next();
        pos++;
        if (iterator->isDone())
        {
            iterator = nullptr;
            pos = 0;
        }
        return *this;
	}



	friend bool operator==(IteratorKmerH5 const& lhs, IteratorKmerH5 const& rhs)
	{
		if (!lhs.iterator || !rhs.iterator)  {  if (!lhs.iterator && !rhs.iterator) {  return true; } else {  return false;  } }
		return rhs.pos == lhs.pos;
	}

	friend bool operator!=(IteratorKmerH5 const& lhs, IteratorKmerH5 const& rhs)  {  return !(lhs == rhs);  }

private:
	Iterator<Kmer<>::Count>* iterator;
	unsigned long pos;
	Kmer<>::ModelCanonical model;
	//KmerLca _item;
	//u_int64_t _item;
};

////template<typename Key>
//class IteratorGzMPHFWrapperKey
//{
//public:
//	IteratorGzMPHFWrapperKey(){}
//	IteratorGzMPHFWrapperKey (Iterator<KmerLca>* iterator) : iterator(iterator) {}
//	IteratorGzMPHFWrapperKey(IteratorGzMPHFWrapperKey& copy){}
//
//	IteratorGzMPHFAdaptatorKey begin() const  {  return IteratorGzMPHFAdaptatorKey (iterator); }
//	IteratorGzMPHFAdaptatorKey end  () const  {  return IteratorGzMPHFAdaptatorKey ();         }
//    size_t        size () const  {  return 0;                        }
//
//private:
//    // noncopyble // FIXME: made it copyable because boophf needed it; need to see if it's correct
//    //iterator_wrapper(iterator_wrapper const&);
//    //iterator_wrapper& operator=(iterator_wrapper const&);
//    Iterator<KmerLca>* iterator;
//};


//template<typename Key>
class IteratorKmerH5Wrapper
{
public:
	IteratorKmerH5Wrapper(){}
	IteratorKmerH5Wrapper (Iterator<Kmer<>::Count>* iterator, const int kmer_size) : iterator(iterator), kmer_size(kmer_size) {}
	IteratorKmerH5Wrapper(IteratorKmerH5& copy){}

	IteratorKmerH5 begin() const  {  return IteratorKmerH5 (iterator, kmer_size); }
	IteratorKmerH5 end  () const  {  return IteratorKmerH5 ();         }
    size_t        size () const  {  return 0;                        }

private:
    // noncopyble // FIXME: made it copyable because boophf needed it; need to see if it's correct
    //iterator_wrapper(iterator_wrapper const&);
    //iterator_wrapper& operator=(iterator_wrapper const&);
    Iterator<Kmer<>::Count>* iterator;
	Kmer<>::ModelCanonical model;
	int kmer_size;
};








#endif /* QUASI_DICTIONNARY_SRC_ITERATORKMERH5_H_ */
