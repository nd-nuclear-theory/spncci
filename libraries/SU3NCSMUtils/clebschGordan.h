#ifndef clebschGordan_h
#define clebschGordan_h
#include <UNU3SU3/UNU3SU3Basics.h>
#include <SU3NCSMUtils/factorial.h>
#include <SU3NCSMUtils/CTuple.h>
#include <LookUpContainers/HashFixed.h>


#include <cstdlib>
#include <iostream>

using std::cout;
using std::endl;

double clebschGordan( int a2, int al2, int b2, int bt2, int c2, int gm2 );
double racah( int a2, int b2, int e2, int d2, int c2, int f2 );
double wigner3jm( int a2,  int b2,  int c2, 
                  int al2, int bt2, int gm2 );
double wigner6j( int a2, int b2, int c2, 
                 int d2, int e2, int f2 );
double wigner9j( int a2, int b2, int c2, 
                 int d2, int e2, int f2,
                 int g2, int h2, int j2 );
double racahU( int a2, int b2, int e2, int d2, int c2, int f2 );
double unitary9j( int a2, int b2, int c2, 
                  int d2, int e2, int f2,
                  int g2, int h2, int j2 );

inline bool isOdd(int a)
{
    return (a & 0x0001 == 1);
}

inline bool isOdd2(int a)
{
    return (a & 0x0002 == 2);
}

inline bool isAllowed(int a2, int b2, int c2)
{
    return (!isOdd(a2 + b2 + c2) && c2 >= abs(a2 - b2) && c2 <= a2 + b2 );
}

inline int multiplicity(int a2, int b2, int c2)
{
    if (isAllowed(a2, b2, c2))
        return (1);
    else 
        return (0);
}

inline double logDelta2(int a2, int b2, int c2)
{
    return (  logFact((a2 + b2 - c2)/2)
             +logFact((c2 + a2 - b2)/2)
             +logFact((b2 + c2 - a2)/2)
             -logFact((a2 + b2 + c2)/2 + 1) );
}

struct hashCTuple
{
  inline std::size_t operator() (const CTuple<SU2::LABEL, 9>& c) const
  {
	 return (std::size_t) c[0] | ((std::size_t)c[1]<<4)| ((std::size_t)c[2]<<8)| ((std::size_t)c[3]<<12)| ((std::size_t)c[4]<<16)| ((std::size_t)c[5]<<20)| ((std::size_t)c[6]<<24) |  ((std::size_t)c[7]<<30) |  ((std::size_t)c[8]<<34);
	 //	 return foo | c[0] | (c[1]<<4)| (c[2]<<8)| (c[3]<<12);
	 //	 return foo | c[0] | (c[1]<<8)| (c[2]<<16)| (c[3]<<24);
  }
};

struct DontDoAnything
{
	inline std::size_t operator() (uint64_t const & c) const
	{
		return c;
	}
};

template <typename DOUBLE = float, uint8_t maxParam = 65>
class CWig9jLookUpTable
{
	typedef HashFixed<uint64_t, DOUBLE, DontDoAnything> CACHE;
	typedef CACHE* WIGNER_9J_TABLE;

	public:
	CWig9jLookUpTable(): wig9jTable_((WIGNER_9J_TABLE)new char[sizeof(CACHE)*maxParam])
	{
		for (int32_t i = 0; i < maxParam; ++i)
		{
			new(wig9jTable_ + i) CACHE(4096);
		}
	}

	~CWig9jLookUpTable()
	{
		for (int32_t i = 0; i < maxParam; ++i)
		{
			(wig9jTable_ + i)->~CACHE();
		}
		delete[] (char*)(wig9jTable_);
	}


	DOUBLE GetWigner9j(SU2::LABEL a1, SU2::LABEL a2, SU2::LABEL a3, SU2::LABEL a4, SU2::LABEL a5, SU2::LABEL a6, SU2::LABEL a7, SU2::LABEL a8, SU2::LABEL a9)
	{
		assert(a1 <= maxParam);
		CACHE& wig9j_a1 = wig9jTable_[a1];
		typename CACHE::iterator wig9j = wig9j_a1.find(INT64BIT_KEY(a2, a3, a4, a5, a6, a7, a8, a9));

		if (wig9j == wig9j_a1.end())
		{
			DOUBLE value = wigner9j(a1, a2, a3, a4, a5, a6, a7, a8, a9);
			wig9j_a1.insert(INT64BIT_KEY(a2, a3, a4, a5, a6, a7, a8, a9), value);
			return value;
		}
		else
		{ 
			return *wig9j;
		}
	}
	private:
/** Construct integer key from the values of \f$L_{i}, S_{i}, L_{f}, S_{f}\f$. */
	inline uint64_t INT64BIT_KEY(SU2::LABEL a2, SO3::LABEL a3, SU2::LABEL a4, SU2::LABEL a5, SU2::LABEL a6, SU2::LABEL a7, SU2::LABEL a8, SU2::LABEL a9) const 
	{ 
		return ((uint64_t)a2 << 56) | ((uint64_t)a3 << 48) | ((uint64_t)a4 << 40) | ((uint64_t)a5 << 32) | ((uint64_t)a6 << 24) | ((uint64_t)a7 << 16) | ((uint64_t)a8 << 8) | (uint64_t)a9; 
	}
	WIGNER_9J_TABLE wig9jTable_;
};

#endif /* clebschGordan */
