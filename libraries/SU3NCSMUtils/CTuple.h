#ifndef CTuple_h
#define CTuple_h
#include <cstring>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
// Class CTuple<T, n> implements an array of N elements of T datatype.  As a
// consequence, sizeof(CTuple<T, N>) = N*sizeof(T)
template<class T, size_t N>
class CTuple {
    private:
	T m_array[N];
    public:
	typedef T Data_Type;
	enum {NELEMENTS = N};

	CTuple(const T* array) {memcpy(m_array, array, N*sizeof(T));}
	CTuple(T* array) {memcpy(m_array, array, N*sizeof(T));}
	CTuple(const CTuple<T, N>& Tuple) {memcpy(m_array, Tuple.Get(), N*sizeof(T));}
	CTuple(const T& Val) {for (size_t i = 0; i < N; m_array[i++] = Val){}}
	//	Warning .... this may not work for array of doubles!!!
	CTuple() {memset(m_array, 0, N*sizeof(T));};
	inline void Set(const T* array) {memcpy(m_array, array, N*sizeof(T));}
	inline void Set(const CTuple<T, N>& Tuple) {memcpy(m_array, Tuple.Get(), N*sizeof(T));}
	inline const T* Get() const {return m_array;}
	inline T* begin() {return m_array;}
	inline T* end() {return m_array+N;}

	inline const T* begin() const {return m_array;} 
	inline const T* end() const {return m_array+N;}

	T& operator[](size_t i) {return m_array[i];}
	const T& operator[](size_t i) const {return m_array[i];}
	inline bool operator!=(const CTuple<T, N>& Tuple) const {return (memcmp(m_array, Tuple.Get(), N*sizeof(T)) != 0); }
	inline bool operator==(const CTuple<T, N>& Tuple) const {return (memcmp(m_array, Tuple.Get(), N*sizeof(T)) == 0); }
	inline bool operator<(const CTuple<T, N>& Tuple) const  {return (memcmp(m_array, Tuple.Get(), N*sizeof(T)) < 0);}

	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		for (int i = 0; i < N; ++i)
		{
        	ar & m_array[i];
		}
    }
	
};

template<class T, size_t N>
inline std::istream& operator>>(std::istream& is, CTuple<T, N>& Indices)  
{
	for (size_t i = 0; i < N; ++i)
	{
		is >> Indices[i];
	}
	return is;
}
#endif

// Struct CCompareLabels implements various sorting criterions for CTuple and
// different set of containers that has class CTuple<n> as a template parameter.
// Example of usage: CSU3Master.cpp line 158 and CSU3Master.h line 99
/* 
template<class T>
struct CCompareLabels 
{
    bool operator() (const T& Labels1, const T& Labels2) { return (memcmp(&Labels1, &Labels2, sizeof(T)) < 0); }
    bool operator() (const std::pair<T, double>& val1, const T& val2) { return (memcmp(&val1.first, &val2, sizeof(T)) < 0); }
    bool operator() (const std::pair<T, std::vector<std::pair<CLabels<4>, double> > >& val1, const T& val2)
    {return (memcmp(&val1.first, &val2, sizeof(T)) < 0); } 
    bool operator() (
	    const std::pair<T, std::vector<std::pair<CLabels<4>, double> > >& val1, 
	    const std::pair<T, std::vector<std::pair<CLabels<4>, double> > >& val2) 
	    {return (memcmp(&val1.first, &val2.first, sizeof(T)) < 0); }
};
*/
