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

	IteratorKmerH5(Iterator<Kmer<>::Count>* iterator)  : iterator(iterator), pos(0)  { iterator->first();}

	u_int64_t  operator*(){
//		cout<<iterator<<endl;
//		u_int64_t res = *new u_int64_t(iterator->item().value.getVal());
//		return *new u_int64_t(iterator->item().value.getVal());
//		cout<<"return "<<(u_int64_t)iterator->item().value.getVal()<<" "<<res<<endl;
//		return res;
		return iterator->item().value.getVal();
//		Kmer<>::Count& count = iterator->item();
//		return count.value;
//		cout<<"count = "<<count.value<<" value = "<<oahash(count.value)<<endl;
//		u_int64_t kmer_int_value=oahash(count.value);
//		return *new u_int64_t(kmer_int_value); //QUESTION: pourquoi utiliser ce new ici. (sinon warning et valeur renvoyée = fausse, meme avec un return 42)
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
};



class IteratorKmerH5Wrapper
{
public:
	IteratorKmerH5Wrapper(){}
	IteratorKmerH5Wrapper (Iterator<Kmer<>::Count>* iterator) : iterator(iterator) {}
	IteratorKmerH5Wrapper(IteratorKmerH5& copy){}

	IteratorKmerH5 begin() const  {   return IteratorKmerH5 (iterator); }
	IteratorKmerH5 end  () const  {  return IteratorKmerH5 ();         }
    size_t        size () const  {  return 0;                        }

private:
    Iterator<Kmer<>::Count>* iterator;
};



//
////template<typename Key>
//class IteratorKmerH5Count : public std::iterator<std::forward_iterator_tag, const Kmer<>::Count>//, const Key>
//{
//public:
//	IteratorKmerH5Count()  : iterator(0), pos(0) {}
//
//	IteratorKmerH5Count(Iterator<Kmer<>::Count>* iterator)  : iterator(iterator), pos(0)  { iterator->first();}
//
//	CountNumber  const& operator*(){
////		cout<<iterator<<endl;
//		return iterator->item().abundance;
//	}
//
//	IteratorKmerH5Count& operator++()
//	{
//        iterator->next();
//        pos++;
//        if (iterator->isDone())
//        {
//        	iterator = nullptr;
//            pos = 0;
//        }
//        return *this;
//	}
//
//
//
//	friend bool operator==(IteratorKmerH5Count const& lhs, IteratorKmerH5Count const& rhs)
//	{
//		if (!lhs.iterator || !rhs.iterator)  {  if (!lhs.iterator && !rhs.iterator) {  return true; } else {  return false;  } }
//		return rhs.pos == lhs.pos;
//	}
//
//	friend bool operator!=(IteratorKmerH5Count const& lhs, IteratorKmerH5Count const& rhs)  {  return !(lhs == rhs);  }
//
//private:
//	Iterator<Kmer<>::Count>* iterator;
//	unsigned long pos;
//};
//
//
//
//class IteratorKmerH5CountWrapper
//{
//public:
//	IteratorKmerH5CountWrapper(){}
//	IteratorKmerH5CountWrapper (Iterator<Kmer<>::Count>* iterator) : iterator(iterator) {}
//	IteratorKmerH5CountWrapper(IteratorKmerH5Count& copy){}
//
//	IteratorKmerH5Count begin() const  {   return IteratorKmerH5Count (iterator); }
//	IteratorKmerH5Count end  () const  {  return IteratorKmerH5Count ();         }
//    size_t        size () const  {  return 0;                        }
//
//private:
//    Iterator<Kmer<>::Count>* iterator;
//};
//



#endif /* QUASI_DICTIONNARY_SRC_ITERATORKMERH5_H_ */
