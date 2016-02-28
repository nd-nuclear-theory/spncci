#ifndef UNU3SU3basics_h
#define UNU3SU3basics_h
#include <SU3NCSMUtils/CTuple.h>
#include <vector>
#include <map>
#include <fstream>
#include <bitset>
#include <cmath>
#include <limits>

#include <cassert>
#include <iostream>

//#define USER_FRIENDLY_OUTPUT		

const double PRECISION_LIMIT = std::numeric_limits<float>::epsilon();

// Class U3 generates and stores a list of HO single particle states for each
// active shell up to nmax (argument of constructor)
//
// Secondly this class defines types U3::LABELS and U3:SPS 
// whose usage makes code more readable
class U3 {
	public:
	typedef CTuple<unsigned char, 3> LABELS;
	typedef std::vector<U3::LABELS> SPS;
	enum {NZ, NX, NY};
};

// This class was created mainly because I want to use UN::LABELS in my code in
// order to make it more readable
class UN {
	public:
//////////////////////////////////////////////////////////////////////////////////////////	
//	New naming is better. One can use UN::SU3 and UN::SU3xSU2 which are different from ordinary SU3 and SU3xSU2
//	due to the outer multiplicity mult arising from reduction U(N)>SU(3)xSU(2)
	struct SU3 {
		typedef unsigned short mult_type; // and I reserve mult = 0xFFFF for vacuum
		typedef unsigned char mu_type;
		typedef unsigned char lm_type;
//		unsigned long mult; // unsigned long is too much, e.g. 64Ge Nmax=10 ==> max(mult) == 567 ==> unsigned short just fine
		mult_type mult; // 64Ge Nmax=14 ==> max(mult) == 1493 ==> unsigned short will be just fine
		lm_type lm;
		mu_type mu;
		SU3(const mult_type mult_, const lm_type lm_, const mu_type mu_): mult(mult_), lm(lm_), mu(mu_) {}
//		SU3(const unsigned char lm_, const unsigned char mu_): mult(1), lm(lm_), mu(mu_) {} // implicitly set mult=1
		SU3(): mult(+1), lm(0), mu(0) {} 
		unsigned long C2() const { return (lm*lm + lm*mu + mu*mu + 3*(lm + mu)); }
		bool operator==(const SU3& unsu3) const {return (unsu3.mult == mult && unsu3.lm == lm && unsu3.mu == mu);}
//		inline bool operator< (const SU3MULT& rhs) const { return ( (lm < rhs.lm) || ((lm == rhs.lm && mu < rhs.mu)) ); }
//		inline bool operator< (const SU3MULT& rhs) const { return (C2() < rhs.C2()); }
	};
	

	struct SU3xSU2: public SU3 {
/*		
 *	Intel compiler had some issues with IDENTITY_OPERATOR [undefined symbol
 *	when linking] and so I decided to comment definition of IDENTITY_OPERATOR
 *	and use instead SU3xSU2::LABELS(+1, 0, 0, 0) instead

		//	The following construction allows me to use UN::SU3xSU2::VACUUM as a parameter of UN::SU3xSU2 constructor
		//	Example: 
		//	if (A == 0) 
		//	{
		//		UN::SU3xSU2 Irrep(UN::SU3xSU2::VACUUM);
		//	}
		struct SetVacuum {
			const static UN::SU3::mult_type mult = +1;
			const static UN::SU3::lm_type lm = 0;
			const static UN::SU3::mu_type mu = 0;
			const static unsigned char S2 = 0;
		};
		static SetVacuum VACUUM;
		SU3xSU2(SetVacuum vacuum) {mult = vacuum.mult; lm = vacuum.lm; mu = vacuum.mu; S2 = vacuum.S2;} 
*/
		unsigned char S2;
		SU3xSU2(const SU3& SU3_, const unsigned char S2_): SU3(SU3_), S2(S2_) {}
		SU3xSU2(): SU3(), S2(0) {}
		SU3xSU2(const SU3::mult_type mult_, SU3::lm_type lm_, SU3::mu_type mu_, const unsigned char S2_): SU3(mult_, lm_, mu_), S2(S2_) {}

		inline bool HasSameLabelsAsVacuum() const {return ((mult == +1) && (lm == 0) && (mu == 0) && (S2 == 0));}

		bool operator==(const SU3xSU2& unsu3su2) const {return (unsu3su2.mult == mult && unsu3su2.lm == lm && unsu3su2.mu == mu && unsu3su2.S2 == S2);}
	};

	typedef std::vector<UN::SU3> SU3_VEC;
	typedef std::vector<UN::SU3xSU2> SU3xSU2_VEC;
//////////////////////////////////////////////////////////////////////////////////////////	

	public:


	typedef std::vector<char> LABELS; 
	typedef std::vector<unsigned char> BASIS_STATE_WEIGHT_VECTOR;
	typedef std::map<U3::LABELS, unsigned long>  U3MULT_LIST;

	typedef std::pair<std::bitset<32>, std::bitset<32> > BASIS_STATE_32BITS; 	// ho shells: 0, 1, 2, 3, 4, 5, 6
	typedef std::pair<std::bitset<64>, std::bitset<64> > BASIS_STATE_64BITS; 	// ho shells: 7, 8, 9

//	TODO: one may want to implement representation of states in terms of indices of occupied single-particle-states	
//	However it is not the highest priority right now 
	typedef std::pair<std::bitset<128>, std::bitset<128> > BASIS_STATE_128BITS; // ho shells: 10, 11, 12, 13, 14 
	typedef std::pair<std::bitset<256>, std::bitset<256> > BASIS_STATE_256BITS; // ho shells: 15, 16, 17, 18, 19, 20, 21 
	public:
	template <typename UN_BASIS_STATE>
	struct CMP_BASIS_STATES
	{
		public:
		bool operator() (const UN_BASIS_STATE& s1, const UN_BASIS_STATE& s2)
		{
			return (s1.first.to_ulong() < s2.first.to_ulong()  || (!(s2.first.to_ulong() < s1.first.to_ulong()) && s1.second.to_ulong() < s2.second.to_ulong()));
		}
	};

	template<typename UN_BASIS_STATE>
	struct BASIS_STATES_EQUAL
	{
		public:
		bool operator() (const UN_BASIS_STATE& s1, const UN_BASIS_STATE& s2)
		{
			return (s1.first.to_ulong() == s2.first.to_ulong() && s1.second.to_ulong() == s2.second.to_ulong());
		}
	};

	public:
	static unsigned long dim(const UN::LABELS& f); // Calculate dimension of U(N) irrep with [f]=[f_{1},\dots,f_{N}] Young shape 
};

template <>
struct UN::CMP_BASIS_STATES<UN::BASIS_STATE_128BITS>
{
	public:
	bool operator() (const BASIS_STATE_128BITS& s1, const BASIS_STATE_128BITS& s2)
	{
		return (s1.first.to_string<char,std::char_traits<char>,std::allocator<char> >() < s2.first.to_string<char,std::char_traits<char>,std::allocator<char> >()  || (!(s2.first.to_string<char,std::char_traits<char>,std::allocator<char> >() < s1.first.to_string<char,std::char_traits<char>,std::allocator<char> >()) && s1.second.to_string<char,std::char_traits<char>,std::allocator<char> >() < s2.second.to_string<char,std::char_traits<char>,std::allocator<char> >()));
	}
};

template<>
struct UN::BASIS_STATES_EQUAL<UN::BASIS_STATE_128BITS>
{
	public:
	bool operator() (const BASIS_STATE_128BITS& s1, const BASIS_STATE_128BITS& s2)
	{
		return (s1.first.to_string<char,std::char_traits<char>,std::allocator<char> >() == s2.first.to_string<char,std::char_traits<char>,std::allocator<char> >() && s1.second.to_string<char,std::char_traits<char>,std::allocator<char> >() == s2.second.to_string<char,std::char_traits<char>,std::allocator<char> >());
	}
};

//	TODO: this will be awfully slow ... one may consider a better implementation
template <>
struct UN::CMP_BASIS_STATES<UN::BASIS_STATE_256BITS>
{
	public:
	bool operator() (const BASIS_STATE_256BITS& s1, const BASIS_STATE_256BITS& s2)
	{
		return (s1.first.to_string<char,std::char_traits<char>,std::allocator<char> >() < s2.first.to_string<char,std::char_traits<char>,std::allocator<char> >()  || (!(s2.first.to_string<char,std::char_traits<char>,std::allocator<char> >() < s1.first.to_string<char,std::char_traits<char>,std::allocator<char> >()) && s1.second.to_string<char,std::char_traits<char>,std::allocator<char> >() < s2.second.to_string<char,std::char_traits<char>,std::allocator<char> >()));
	}
};

template<>
struct UN::BASIS_STATES_EQUAL<UN::BASIS_STATE_256BITS>
{
	public:
	bool operator() (const BASIS_STATE_256BITS& s1, const BASIS_STATE_256BITS& s2)
	{
		return (s1.first.to_string<char,std::char_traits<char>,std::allocator<char> >() == s2.first.to_string<char,std::char_traits<char>,std::allocator<char> >() && s1.second.to_string<char,std::char_traits<char>,std::allocator<char> >() == s2.second.to_string<char,std::char_traits<char>,std::allocator<char> >());
	}
};

inline std::ostream& operator << (std::ostream& os, const UN::SU3& labels)
{
	return (os << (int)labels.mult << "(" << (int)labels.lm << " " << (int)labels.mu << ")");
}
inline std::ostream& operator << (std::ostream& os, const UN::SU3xSU2& labels)
{
	return (os << (const UN::SU3&)labels << (int)labels.S2);
}

inline std::istream& operator >> (std::istream& is, UN::SU3& labels)
{
		int tmp;
		is >> tmp;
		labels.mult = tmp;
		is >> tmp;
		labels.lm = tmp;
		is >> tmp;
		labels.mu = tmp;
		return (is);
}

inline std::istream& operator >> (std::istream& is, UN::SU3xSU2& labels)
{
		int tmp;
		is >> tmp;
		labels.mult = tmp;
		is >> tmp;
		labels.lm = tmp;
		is >> tmp;
		labels.mu = tmp;
		is >> tmp;
		labels.S2 = tmp;
		return (is);
}


/*----------------------------------------------------*
 * 				FORWARD DECLARATIONS
 *----------------------------------------------------*/ 				
namespace SU2 {
	typedef unsigned char LABEL;
	typedef char CANONICAL;
}
typedef std::vector<SU2::LABEL> SU2_VEC;

namespace U1 {
	typedef char LABEL;
};

namespace U1 {
	struct LabelGenerator {
		private:
			SU2::LABEL m_S2;
		public:
			LabelGenerator(const SU2::LABEL& S2): m_S2(S2+2) {}
			U1::LABEL operator() () {return (m_S2 -= 2);}
	};
};

namespace SU2 {
	class IRREP 
	{
	private:
		SU2::LABEL m_S2;
	public:
		IRREP(const SU2::LABEL& S2): m_S2(S2) {}
		operator SU2::LABEL() {return m_S2;}
		inline size_t dim() {return (m_S2 + 1);}
		size_t GetU1Labels(std::vector<U1::LABEL>& Sigma) 
		{
			Sigma.resize(0);
			Sigma.reserve(dim());
			std::generate_n(back_inserter(Sigma), dim(), U1::LabelGenerator(m_S2));
			return Sigma.size();
		}
	};	

	inline size_t dim(const SU2::LABEL& S2) {return (S2 + 1);}
	// WARNING: always make sure that S1, S2, and S3 holds 2*angular momenta.
	inline size_t mult(const SU2::LABEL& S1, const SU2::LABEL& S2, const SU2::LABEL& S3) {return (S3 >= abs(S1 - S2) && S3 <= S1 + S2) && !((S1 + S2 + S3) & 0x0001);}
	inline size_t GetU1Labels(const SU2::LABEL& S2, std::vector<U1::LABEL>& Sigma) { return (SU2::IRREP(S2)).GetU1Labels(Sigma); }
	inline void Couple(const SU2::LABEL J12, const SU2::LABEL J22, std::vector<SU2::LABEL>& vJ32) 
	{
		for (size_t j32 = abs(J12-J22); j32 <= J12+J22; j32 += 2)
		{
			vJ32.push_back(j32);
		}
	}
/*	
	inline size_t Couple(const SU2::LABEL J12, const SU2::LABEL J22, std::vector<SU2::LABEL>& vJ32) 
	{
		size_t minJJ = abs(J12 - J22);
		vJ32.reserve((J12+J22 - minJJ + 2)<<2);
		vJ32.push_back(minJJ);
		for (; vJ32.back() < J12+J22; vJ32.push_back(vJ32.back() + 2));
	}
*/
};

namespace SO3 = SU2;
namespace _SO3 = SO3; // to be able to use use type SO3::LABEL in namespace SU3::SO3 by using _SO3::LABEL

namespace SU3 {
	typedef double WIGNER;
	typedef char rho_type;
	typedef char RHO_MAX_TYPE;
	typedef unsigned char mu_type;
	typedef unsigned char lm_type;

	struct LABELS 
	{

		LABELS(const SU3::LABELS& R): rho(R.rho), lm(R.lm), mu(R.mu) {} 
		LABELS(int rhot, int lmt, int mut):rho(rhot), lm(lmt), mu(mut) {}
		LABELS(int lmt, int mut):rho(+1), lm(lmt), mu(mut) {} // if no rho in input ==> implicitly set rhomax/rho = 1
		LABELS(): rho(+1), lm(0), mu(0) {}

		inline int C2() const { return (lm*lm + lm*mu + mu*mu + 3*(lm + mu)); }

		inline bool operator< (const SU3::LABELS& rhs) const 
		{
			return (lm < rhs.lm) || (lm == rhs.lm && mu < rhs.mu);
		}

		inline bool operator== (const SU3::LABELS& rhs) const {return ((lm == rhs.lm) && (mu == rhs.mu) && (rho == rhs.rho));}
		inline bool operator!= (const SU3::LABELS& rhs) const {return ((lm != rhs.lm) || (mu != rhs.mu) || (rho != rhs.rho));}
		inline size_t dim() {return (size_t)((lm + 1.0) * (mu + 1.0) * ((lm + mu)/2.0 + 1.0));}
		rho_type rho; // this variable should rahter be  called rhomax 
		lm_type lm;
		mu_type mu;
	};


	namespace SO3 
	{
		struct LABELS {
			LABELS(const char& k, const _SO3::LABEL& l2): K(k), L2(l2) {}
			LABELS(): K(0), L2(0) {}
			inline bool operator< (const SU3::SO3::LABELS& rhs) const {return (K < rhs.K) || ((K == rhs.K) && (L2 < rhs.L2));}
			inline bool operator== (const SU3::SO3::LABELS& rhs) const {return (K == rhs.K && L2 ==  rhs.L2);}
			char K;			//	kappa multiplicity in SU(3)>SO(3) reduction
			_SO3::LABEL L2;	//	SO(3)
		};
	};

	namespace SU2U1 {
		struct  LABELS {
			LABELS(const U1::LABEL& eps_, const SU2::LABEL& LM_): eps(eps_), LM(LM_) {}
			LABELS(): eps(0), LM() {}
			inline bool operator< (const SU3::SU2U1::LABELS& rhs) const {return (eps < rhs.eps) || ((eps == rhs.eps) && LM < rhs.LM);}
			inline bool operator==(const SU3::SU2U1::LABELS& rhs) const {return ((eps == rhs.eps) && (LM == rhs.LM));}

			U1::LABEL eps;		// epsilon	...		U(1) 
			SU2::LABEL LM;		// 2*LAMBDA	...		SU(2)
		};

		struct LabelGenerator 
		{
		private:
			const int mu, epstop;
			int p, q;
		public:
			LabelGenerator(const SU3::LABELS& lmmu): mu(lmmu.mu), epstop(2*lmmu.lm+lmmu.mu) , p(0), q(-1){}
			inline SU3::SU2U1::LABELS operator() ()
			{
				if (++q > mu) {
					p++;
					q = 0;
				}
				return SU3::SU2U1::LABELS(epstop - 3*(p + q), (mu + p - q));
			} 
		};
	};

	struct CANONICAL: SU3::SU2U1::LABELS { // full SU(3)>SU(2)xU(1)>U(1) labeling: epsilon, 2LM, 2*M_LM
		CANONICAL(const SU3::SU2U1::LABELS& SU2U1Labels, const U1::LABEL& m): SU3::SU2U1::LABELS(SU2U1Labels) {M_LM = m;}
		CANONICAL(const U1::LABEL& eps_, const SU2::LABEL& LM_, const U1::LABEL& M_LM_): SU3::SU2U1::LABELS(eps_, LM_) {M_LM = M_LM_;}
		CANONICAL() {M_LM = 0;}
		inline bool operator< (const SU3::CANONICAL& rhs) const {return (SU2U1::LABELS::operator<(rhs) || ((SU2U1::LABELS::operator==(rhs)) && M_LM < rhs.M_LM));}
		inline bool operator==(const SU3::CANONICAL& rhs) const {return (SU2U1::LABELS::operator==(rhs) && M_LM == rhs.M_LM);}
		U1::LABEL M_LM; 	//	2*M_LAMDBA		
	};

	inline int LM(const U3::LABELS& U3Labels) {return (U3Labels[U3::NZ] - U3Labels[U3::NX]);}
	inline int MU(const U3::LABELS& U3Labels) {return (U3Labels[U3::NX] - U3Labels[U3::NY]);}

	inline void HWS(const SU3::LABELS& SU3Labels, SU3::CANONICAL& HwsLabels) 
	{
		HwsLabels.eps = 2*SU3Labels.lm + SU3Labels.mu;
		HwsLabels.M_LM = HwsLabels.LM  = SU3Labels.mu;
	}

	inline void HWS(const  U3::LABELS&  U3Labels, SU3::CANONICAL& HwsLabels)
	{
		HwsLabels.eps = 2*SU3::LM(U3Labels) + SU3::MU(U3Labels);
		HwsLabels.M_LM = HwsLabels.LM  = SU3::MU(U3Labels);
	}


	inline int kmax(const SU3::LABELS& RhoLmMu, const int L) 
	{return(std::max(0,((short)RhoLmMu.lm+(short)RhoLmMu.mu+2-L)/2)-std::max(0,((short)RhoLmMu.lm+1-L)/2)-std::max(0,((short)RhoLmMu.mu+1-L)/2));}

	
	template<typename T>
	inline size_t dim(const T& RhoLmMu) {return (size_t)((RhoLmMu.lm + 1.0) * (RhoLmMu.mu + 1.0) * ((RhoLmMu.lm + RhoLmMu.mu)/2.0 + 1.0));}
//	int mult (const SU3::LABELS& RhoLmMu1, const SU3::LABELS& RhoLmMu2, const SU3::LABELS& RhoLmMu3);

	template<class T1, class T2, class T3>
	int mult(const T1& RhoLmMu1, const T2& RhoLmMu2, const T3& RhoLmMu3)
	{
		int lm1(RhoLmMu1.lm), mu1(RhoLmMu1.mu), lm2(RhoLmMu2.lm), mu2(RhoLmMu2.mu), lm3(RhoLmMu3.lm), mu3(RhoLmMu3.mu);
	    unsigned result = 0;

    	int phase, redphase, invl1, invl2;

	    phase = lm1 + lm2 - lm3 - mu1 - mu2 + mu3;

    	if (phase % 3 == 0)
	    {
    	   int l1, l2, l3, m1, m2, m3;

        	redphase = phase / 3;
	        if (phase >= 0){
    	        l1 = lm1; l2 = lm2; l3 = lm3;
        	    m1 = mu1; m2 = mu2; m3 = mu3;
	        } else {
    	        l1 = mu1; l2 = mu2; l3 = mu3;
        	    m1 = lm1; m2 = lm2; m3 = lm3;
            	redphase = -redphase;
	        }
    	    phase = redphase + m1 + m2 - m3;

	        invl1 = std::min(l1 - redphase, m2 );
    	    if (invl1 < 0) return (result);

	        invl2 = std::min(l2 - redphase, m1 );
    	    if (invl2 < 0) return (result);

        	result = std::max(std::min(phase,invl2) - std::max(phase-invl1,0) + 1, 0);
	    }
    	return (result);
	}
/*
	    if (phase % 3 == 0)
    	{
	        int invl1, invl2, redphase = phase / 3;
	        if (phase >= 0){
	        	invl1 = std::min(lm1 - redphase, mu2);
    	    	if (invl1 < 0) 
				{
					return 0;
				}
		        invl2 = std::min(lm2 - redphase, mu1);
	    	    if (invl2 < 0) 
				{
					return 0;
				}
    	    	phase = mu1 + mu2 - mu3 + redphase;
	        } else {
	        	invl1 = std::min(mu1 + redphase, lm2);
    		    if (invl1 < 0) 
				{
					return 0;
				}
		        invl2 = std::min(mu2 + redphase, lm1);
    		    if (invl2 < 0) 
				{
					return 0;
				}
    	    	phase = lm1 + lm2 - lm3 - redphase;
	        }
	        return std::max(std::min(phase, invl2) - std::max(phase - invl1, 0) + 1, 0);
    	}
	    return 0;
*/
	size_t GetSO3Labels(const SU3::LABELS& RhoLmMu, std::vector<SU3::SO3::LABELS>& KL);
	size_t GetAngularMomenta(const SU3::LABELS& RhoLmMu, std::vector<SU2::LABEL>& Ls);

	inline size_t GetSU2U1Labels(const SU3::LABELS& RhoLmMu, std::vector<SU3::SU2U1::LABELS>& epsLam)	
	{
		epsLam.resize(0);
		generate_n(back_inserter(epsLam), (RhoLmMu.lm + 1)*(RhoLmMu.mu + 1), SU3::SU2U1::LabelGenerator(RhoLmMu));
		return epsLam.size();
	}
	inline static void toNxNyNz(const SU3::CANONICAL& sps, int& nx, int& ny, int& nz) {nx = (sps.LM  + sps.M_LM)/2; ny = (sps.LM  - sps.M_LM)/2; nz = (sps.eps + sps.LM)/2;}
	size_t GetCanonicalBasis(const SU3::LABELS RhoLmMu, std::vector<SU3::CANONICAL>& SU3CanonicalBasis, bool bSortedOrder = true);

	inline std::ostream& operator << (std::ostream& os, const SU3::LABELS& labels)
	{
#ifdef USER_FRIENDLY_OUTPUT		
		if (labels.rho < 0) {
			return (os << "(" << (int)labels.lm << " " << (int)labels.mu << ")");
		} else {
			return (os << (int)labels.rho << "(" << (int)labels.lm << " " << (int)labels.mu << ")");
		}
#else
		return (os << (int)labels.rho << " " << (int)labels.lm << " " << (int)labels.mu << " ");
#endif		
	}

	inline std::istream& operator >> (std::istream& is, SU3::LABELS& labels)
	{
		int tmp;
		is >> tmp;
		labels.rho = tmp;
		is >> tmp;
		labels.lm = tmp;
		is >> tmp;
		labels.mu = tmp;
		return (is);
	}


	inline std::ostream& operator << (std::ostream& os, const SU3::CANONICAL& labels)
	{
		return (os << (int)labels.eps << "  " << (int)labels.LM << "  " << (int)labels.M_LM);
	}
};

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, SU3::LABELS& irrep, const unsigned int version)
{
   ar & irrep.rho;
   ar & irrep.lm;
   ar & irrep.mu;
}

template<class Archive>
void serialize(Archive & ar, UN::SU3xSU2& irrep, const unsigned int version)
{
   ar & irrep.mult;
   ar & irrep.lm;
   ar & irrep.mu;
   ar & irrep.S2;
}

}	// namespace serialization
}	// namespace boost

typedef std::vector<SU3::LABELS> SU3_VEC;
namespace SU3
{
	template<class T1, class T2>
	size_t Couple(const T1& ir1, const T2& ir2, SU3_VEC& ResultingIrreps)
	{
		SU3::LABELS ir3;
		size_t max = (ir1.lm + ir2.lm + ir1.mu + ir2.mu);
		ResultingIrreps.reserve(max);

	    for (ir3.lm = 0; ir3.lm <= max; ++ir3.lm) 
    	{
       		for (ir3.mu = 0; ir3.mu <= max; ++ir3.mu) 
			{
		    	ir3.rho  = SU3::mult(ir1, ir2, ir3);
			    if (ir3.rho != 0) 
		    	{
		    		ResultingIrreps.push_back(ir3);
				}
			}
    	}	
		return ResultingIrreps.size();
	}
}


namespace SU3xSU2 
{
	struct PHYSICALJ
	{
		PHYSICALJ(char K, char LL, char JJ, char MM): k(K), ll(LL), jj(JJ), mm(MM) {}
		PHYSICALJ(): k(0), ll(0), jj(0), mm(0) {}
		inline bool operator<(const SU3xSU2::PHYSICALJ& rhs) const 
		{
			return (k < rhs.k || (((k == rhs.k) && ll < rhs.ll) || ((k == rhs.k && ll == rhs.ll && jj < rhs.jj) || (k == rhs.k && ll == rhs.ll && jj == rhs.jj && mm < rhs.mm))));
		}
		char k, ll, jj, mm;
	};

	struct LABELS: SU3::LABELS 
	{
/*		
	Intel compiler had some issues with IDENTITY_OPERATOR [undefined symbol
	when linking] and so I decided to comment definition of IDENTITY_OPERATOR
	and use instead SU3xSU2::LABELS(+1, 0, 0, 0) instead

		//	The following construction allows me to use SU3xSU2::LABELS::IDENTITY_OPERATOR as a parameter of SU3xSU2::LABELS constructor
		//	and also IR0 = SU3xSU2::LABELS::IDENTITY_OPERATOR
		//	Example: 
		//	if (bra == ket && dA == 0) 
		//	{
		//		Tensor.push_back(SU3xSU2::LABELS::IDENTITY_OPERATOR);
		//	}
		struct IdentityOperatorLabels {
			const static SU3::rho_type rho = +1;
			const static SU3::lm_type lm = 0;
			const static SU3::mu_type mu = 0;
			const static unsigned char S2 = 0;
		};
		static IdentityOperatorLabels IDENTITY_OPERATOR;
		LABELS(IdentityOperatorLabels identity): SU3::LABELS(identity.rho, identity.lm, identity.mu), S2(identity.S2) {} 
*/

////////////////////////////////////////////////////////
//	adan > 0 ====> creation operator on HO shell n = adan - 1 			==>	1(lm=(adan-1) mu=0) 1/2
//	adan < 0 ====> annihilation operator on HO shell n = -(adan + 1)	==>	1(lm=0 mu=-(adan+1)) 1/2
		LABELS(const char adan):SU3::LABELS(1, (adan > 0) ? (adan - 1):0, (adan > 0) ? 0:-1*(adan + 1)), S2(1) {};

		LABELS(const SU3::LABELS& Su3Labels, const SU2::LABEL& Su2Label): SU3::LABELS(Su3Labels), S2(Su2Label) {}
		LABELS(const SU3::rho_type rho, const SU3::lm_type lm, const SU3::mu_type mu, const SU2::LABEL s2): SU3::LABELS(rho, lm, mu), S2(s2) {}
		LABELS(const SU3::lm_type lm, const SU3::mu_type mu, const SU2::LABEL s2): SU3::LABELS(lm, mu), S2(s2) {}
		LABELS(const UN::SU3xSU2& unsu3su2): SU3::LABELS(unsu3su2.lm, unsu3su2.mu), S2(unsu3su2.S2) {}
		LABELS(): SU3::LABELS() { S2 = 0;}

		inline bool operator<(const SU3xSU2::LABELS& rhs) const {return ((lm < rhs.lm || (lm == rhs.lm && mu < rhs.mu)) || ((lm == rhs.lm && mu == rhs.mu) && S2 < rhs.S2));}
		inline bool operator==(const SU3xSU2::LABELS& rhs) const {return (SU3::LABELS::operator==(rhs) && S2 == rhs.S2);}
		inline bool operator!=(const SU3xSU2::LABELS& rhs) const {return (S2 != rhs.S2 || SU3::LABELS::operator!=(rhs));}

		operator SU2::LABEL() const {return S2;}

		inline size_t dim() {return SU3::LABELS::dim()*SU2::dim(S2);}

		SU2::LABEL S2;
	};

	
	struct CANONICAL: SU3::CANONICAL { // Full SU(3)xSU(2)>SU(2)xU(1)xSU(2)>U(1)xU(1): eps 2LM 2L_LM S Sigma
		CANONICAL(const SU3::CANONICAL& epsLMM, const U1::LABEL& ms): SU3::CANONICAL(epsLMM), Sigma(ms) {}
		CANONICAL(const U1::LABEL& eps_, const SU2::LABEL& LM_, const U1::LABEL& M_LM_, const U1::LABEL& Sigma_): SU3::CANONICAL(eps_, LM_, M_LM_), Sigma(Sigma_) {}
		CANONICAL(): SU3::CANONICAL() {Sigma = 0;}
		inline bool operator<(const SU3xSU2::CANONICAL& rhs) const {return (SU3::CANONICAL::operator<(rhs) || (SU3::CANONICAL::operator==(rhs) && Sigma < rhs.Sigma));}
		U1::LABEL Sigma;
	};
//	Here is the definition of the SU(3)xSU(2) highest 	
//	weight state which is used for the generation of rme.
	inline void HWS(const SU3xSU2::LABELS Su3xSU2Labels, SU3xSU2::CANONICAL& Hws)
	{
		SU3::HWS(Su3xSU2Labels, Hws);
		Hws.Sigma = Su3xSU2Labels.S2;
	}

	size_t GetPhysicalJBasis(const SU3xSU2::LABELS& IR, std::vector<SU3xSU2::PHYSICALJ>& SU3xSU2PhysicalJBasis);
	size_t GetCanonicalBasis(const SU3xSU2::LABELS& RhoLmMuS2, std::vector<SU3xSU2::CANONICAL>& SU3xSU2CanonicalBasis, bool bSortedOrder = true);
	inline size_t dim(const SU3xSU2::LABELS& Ir1) { return SU3::dim((SU3::LABELS)Ir1)*SU2::dim(Ir1.S2);}
	inline std::ostream& operator << (std::ostream& os, const SU3xSU2::CANONICAL& labels)
	{
		return (os << (int)labels.eps << "  " << (int)labels.LM << "  " << (int)labels.M_LM << " " << (int)labels.Sigma);
	}

	inline std::ostream& operator << (std::ostream& os, const SU3xSU2::PHYSICALJ& labels)
	{
		return (os << (int)labels.k << " " << (int)labels.ll << "  " << (int)labels.jj << " " << (int)labels.mm);
	}


	inline std::ostream& operator << (std::ostream& os, const SU3xSU2::LABELS& labels)
	{
#ifdef USER_FRIENDLY_OUTPUT		
		return (os << (const SU3::LABELS&)labels << (int)labels.S2);
#else
		return (os << (const SU3::LABELS&)labels << " " << (int)labels.S2);
#endif		
	}

	inline std::istream& operator >> (std::istream& is, SU3xSU2::LABELS& labels)
	{
		int tmp;
		is >> tmp;
		labels.rho = tmp;
		is >> tmp;
		labels.lm = tmp;
		is >> tmp;
		labels.mu = tmp;
		is >> tmp;
		labels.S2 = tmp;
		return (is);
	}
	size_t Couple(const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, std::vector<SU3xSU2::LABELS>& ResultingIrreps);
};

namespace boost {
namespace serialization {

template<class Archive>
void serialize(Archive & ar, SU3xSU2::LABELS& irrep, const unsigned int version)
{
   ar & irrep.rho;
   ar & irrep.lm;
   ar & irrep.mu;
   ar & irrep.S2;
}

}	// namespace serialization
}	// namespace boost



namespace SU3{
	inline SU3xSU2::LABELS Conjugate(const SU3xSU2::LABELS& ir) {return SU3xSU2::LABELS(ir.rho, ir.mu, ir.lm, ir.S2);}
	inline SU3::LABELS Conjugate(const SU3::LABELS& ir) {return SU3::LABELS(ir.rho, ir.mu, ir.lm);}
};
typedef std::vector<SU3xSU2::LABELS> SU3xSU2_VEC;

size_t Couple(const SU3::LABELS& ir1, const SU3::LABELS& ir2, std::vector<SU3::LABELS>& ResultingIrreps);
#endif
