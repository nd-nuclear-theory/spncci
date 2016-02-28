#ifndef CSU3Master_h
#define CSU3Master_h
#include <UNU3SU3/UNU3SU3Basics.h>
#include <cstdlib>
#include <utility>

const size_t MAX_K = 9;

// Subroutines of original Fortran SU(3) library
extern "C" { 
#ifndef AIX
//		LM1   MU1   LM2   MU2   LM3   MU3   L1    L2    L3   K0MAX K1MAX K2MAX K3MAX
extern void wu3r3w_(int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, double[MAX_K][MAX_K][MAX_K][MAX_K]);
extern void wru3optimized_(int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&, double[], int&);
extern void wzu3optimized_(int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&, double[], int&);
extern void xewu3_(int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int[], int[], int[], int&, int&, int&);
extern void xwu3_(int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int&, int&, double[], int[], int[], int[], int[], int&, int&, int[], int&, int&, int&, int&);
extern void wu39lm_(int&, int& , int&, int&, int& , int& , int& , int&, int&, int&, int&, int&, int& , int& , int& , int&, int&, int&, double[], int&);
extern void blocks_(void);
#else
extern void wu3r3w(int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, 
		double[MAX_K][MAX_K][MAX_K][MAX_K]);
extern void wru3optimized(int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&, double[], int&);
extern void wzu3optimized(int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&,int&, double[], int&);
extern void xewu3(int&, int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int[], int[], int[], int&, int&, int&);
extern void xwu3(int&, int&, int&, int&, int&, int&, int&, int&, int&, double[], int&, int&, double[], int[], int[], int[], int[], int&, int&, int[], int&, int&, int&, int&);
extern void wu39lm(int&, int& , int&, int&, int& , int& , int& , int&, int&, int&, int&, int&, int& , int& , int& , int&, int&, int&, double[], int&);
extern void blocks(void);
#endif
}

// Before the first use of any subroutine from Fortran SU(3) library,
// subroutine blocks_ has to be called. This class implements logic of this
// functionality so that a user of class CSU3Master does not have to take care
// of that.
class CBlocks {
private:
	static bool m_fBlocksInitiated;
public:
	CBlocks() 
	{ 
	    if (!m_fBlocksInitiated) 
	    {
#ifndef AIX		
		blocks_();
#else
		blocks();
#endif		
		m_fBlocksInitiated = true;
	    }
	}
};

// Note: even tough methods LM?(), eps?(), or l?() should return SU2::LABEL, U1::LABEL, SO3::LABEL, they return int so
// that one can use cout >> LM1() 

struct EPS1LM1EPS2LM2 
{
	EPS1LM1EPS2LM2(const SU3::SU2U1::LABELS& eps1lm1, const SU3::SU2U1::LABELS& eps2lm2): epslm1(eps1lm1), epslm2(eps2lm2) {};
	EPS1LM1EPS2LM2(const U1::LABEL& e1, const SU2::LABEL& lm1, const U1::LABEL& e2, const SU2::LABEL& lm2): epslm1(e1, lm1), epslm2(e2, lm2) {};
	EPS1LM1EPS2LM2(): epslm1(), epslm2() {};

	bool operator< (const EPS1LM1EPS2LM2& rhs) const {return (epslm1 < rhs.epslm1) || (!(rhs.epslm1 < epslm1) && epslm2 < rhs.epslm2);};

	inline int eps1() const {return epslm1.eps;};
	inline int eps2() const {return epslm2.eps;};
	inline int LM1() const {return epslm1.LM;};
	inline int LM2() const {return epslm2.eps;};
	inline void eps1(const U1::LABEL& e) {epslm1.eps = e;};
	inline void eps2(const U1::LABEL& e) {epslm2.eps = e;};
	inline void LM1(const SU2::LABEL& LM) {epslm1.LM = LM;};
	inline void LM2(const SU2::LABEL& LM) {epslm2.LM = LM;};
private:
	SU3::SU2U1::LABELS epslm1, epslm2;
};
struct EPS1LM1EPS2LM2EPS3LM3 
{
	EPS1LM1EPS2LM2EPS3LM3 (const SU3::SU2U1::LABELS& eps1lm1, const SU3::SU2U1::LABELS& eps2lm2, const SU3::SU2U1::LABELS& eps3lm3): epslm1(eps1lm1), epslm2(eps2lm2), epslm3(eps3lm3) {};
	EPS1LM1EPS2LM2EPS3LM3 (const U1::LABEL e1, const SU2::LABEL& lm1, const U1::LABEL& e2, const SU2::LABEL& lm2, const U1::LABEL& e3, const SU2::LABEL& lm3):
	epslm1(e1, lm1), epslm2(e2, lm2), epslm3(e3, lm3) {};
	EPS1LM1EPS2LM2EPS3LM3 (): epslm1(), epslm2(), epslm3() {};

	bool operator== (const EPS1LM1EPS2LM2EPS3LM3& rhs) const 
	{
		return epslm1 == rhs.epslm1 && epslm2 == rhs.epslm2 && epslm3 == rhs.epslm3;
	}	


	bool operator< (const EPS1LM1EPS2LM2EPS3LM3& rhs) const 
	{
/*		
		std::cout << (int)epslm1.eps << " " << (int)epslm1.LM << "  " << (int)epslm2.eps << " " << (int)epslm2.LM << "  " << (int)epslm3.eps << " " << (int)epslm3.LM;
		std::cout << " < " << (int)rhs.eps1() << " " << (int)rhs.LM1() << "  " << (int)rhs.eps2() << " " << (int)rhs.LM2() << "  " << (int)rhs.eps3() << " " << (int)rhs.LM3();
		std::cout << " ---> returning " << ((epslm1 < rhs.epslm1) || (!(rhs.epslm1 < epslm1) && epslm2 < rhs.epslm2) || (!(rhs.epslm1 < epslm1) && !(rhs.epslm2 < epslm2) && epslm3 < rhs.epslm3)) << std::endl;
		std::cout <<  (!(rhs.epslm1 < epslm1) && epslm2 < rhs.epslm2) << std::endl;
		std::cout << (!(rhs.epslm1 < epslm1) && !(rhs.epslm2 < epslm2) && epslm3 < rhs.epslm3) << std::endl;
		std::cout << ((!(rhs.epslm1 < epslm1) && epslm2 < rhs.epslm2) || (!(rhs.epslm1 < epslm1) && !(rhs.epslm2 < epslm2) && epslm3 < rhs.epslm3)) << std::endl;
		std::cout << (epslm1 < rhs.epslm1) << std::endl;
		std::cout << ((epslm1 < rhs.epslm1) || ((!(rhs.epslm1 < epslm1) && epslm2 < rhs.epslm2) || (!(rhs.epslm1 < epslm1) && !(rhs.epslm2 < epslm2) && epslm3 < rhs.epslm3))) << std::endl;
*/
		return (epslm1 < rhs.epslm1) 
				|| (
				(epslm1 == rhs.epslm1 && epslm2 < rhs.epslm2) 
				|| (
				(epslm1 == rhs.epslm1 && epslm2 == rhs.epslm2 && epslm3 < rhs.epslm3)));
	}

	inline int eps1() const {return epslm1.eps;};
	inline int eps2() const {return epslm2.eps;};
	inline int eps3() const {return epslm3.eps;};
	inline void eps1(const U1::LABEL& e) {epslm1.eps = e;};
	inline void eps2(const U1::LABEL& e) {epslm2.eps = e;};
	inline void eps3(const U1::LABEL& e) {epslm3.eps = e;};

	inline int LM1() const {return epslm1.LM;};
	inline int LM2() const {return epslm2.LM;};
	inline int LM3() const {return epslm3.LM;};
	inline void LM1(const SU2::LABEL& LM) {epslm1.LM = LM;};
	inline void LM2(const SU2::LABEL& LM) {epslm2.LM = LM;};
	inline void LM3(const SU2::LABEL& LM) {epslm3.LM = LM;};
private:
	SU3::SU2U1::LABELS epslm1, epslm2, epslm3;
};

struct K1L1K2L2 
{
	K1L1K2L2(const char& k1, const SO3::LABEL& l1, const char& k2, const SO3::LABEL& l2): k1l1(k1, l1), k2l2(k2, l2) {};
	K1L1K2L2(const SU3::SO3::LABELS& so3labels1, const SU3::SO3::LABELS& so3labels2): k1l1(so3labels1), k2l2(so3labels2) {};
	K1L1K2L2(): k1l1(), k2l2() {};

	bool operator< (const K1L1K2L2& rhs) const {return k1l1 < rhs.k1l1 || (!(rhs.k1l1 < k1l1) && k2l2 < rhs.k2l2);}

	inline int k1() const {return k1l1.K;};
	inline int k2() const {return k2l2.K;};
	inline void k1(int k) {k1l1.K = k;};
	inline void k2(int k) {k2l2.K = k;};
		
	inline int  l1() const {return k1l1.L2;};
	inline int  l2() const {return k2l2.L2;};
	inline void l1(const SO3::LABEL& l2) {k1l1.L2 = l2;};
	inline void l2(const SO3::LABEL& l2) {k2l2.L2 = l2;};
private:
	SU3::SO3::LABELS k1l1;
	SU3::SO3::LABELS k2l2;
};
	
struct K1L1K2L2K3 
{
	K1L1K2L2K3(const SU3::SO3::LABELS& so3labels1, const SU3::SO3::LABELS& so3labels2, const char k):k1l1(so3labels1), k2l2(so3labels2), K3(k) {}
	K1L1K2L2K3(): k1l1(), k2l2(), K3(0) {};

	bool operator< (const K1L1K2L2K3& rhs) const {
		return ((k1l1 < rhs.k1l1) || (!(rhs.k1l1 < k1l1) && k2l2 < rhs.k2l2) || (!(rhs.k1l1 < k1l1) && !(rhs.k2l2 < k2l2) && K3 < rhs.K3));
	}

	inline int  k1() const {return k1l1.K;};
	inline int  k2() const {return k2l2.K;};
	inline int  k3() const {return K3;};
	inline void k1(const char k) {k1l1.K = k;};
	inline void k2(const char k) {k2l2.K = k;};
	inline void k3(const char k) {K3 = k;};

	inline int  l1() const {return k1l1.L2;};
	inline int  l2() const {return k2l2.L2;};
	inline void l1(const SO3::LABEL& l) {k1l1.L2 = l;};
	inline void l2(const SO3::LABEL& l) {k2l2.L2 = l;};
private:
	SU3::SO3::LABELS k1l1, k2l2;
	char K3;
};

struct K1L1K2L2K3L3 
{
	K1L1K2L2K3L3(const char& k1, const U1::LABEL& l1, const char& k2, const U1::LABEL& l2, const char& k3, const U1::LABEL& l3): 
	k1l1(k1, l1), k2l2(k2, l2), k3l3(k3, l3) {};
	K1L1K2L2K3L3(const SU3::SO3::LABELS& so3labels1, const SU3::SO3::LABELS& so3labels2, const SU3::SO3::LABELS& so3labels3): k1l1(so3labels1), k2l2(so3labels2), k3l3(so3labels3) {};
	K1L1K2L2K3L3(): k1l1(), k2l2(), k3l3() {};

	bool operator< (const K1L1K2L2K3L3& rhs) const {
		return ((k1l1 < rhs.k1l1) || (!(rhs.k1l1 < k1l1) && k2l2 < rhs.k2l2) || (!(rhs.k1l1 < k1l1) && !(rhs.k2l2 < k2l2) && k3l3 < rhs.k3l3));
	}
	bool operator== (const K1L1K2L2K3L3& rhs) const {return (k1l1 == rhs.k1l1 && k2l2 == rhs.k2l2 && k3l3 == rhs.k3l3);}

	inline int k1() const {return k1l1.K;};
	inline int k2() const {return k2l2.K;};
	inline int k3() const {return k3l3.K;};
	inline void k1(int k) {k1l1.K = k;};
	inline void k2(int k) {k2l2.K = k;};
	inline void k3(int k) {k3l3.K = k;};
	
	inline int l1() const {return k1l1.L2;};
	inline int l2() const {return k2l2.L2;};
	inline int l3() const {return k3l3.L2;};
	inline void l1(const SO3::LABEL& l) {k1l1.L2 = l;};
	inline void l2(const SO3::LABEL& l) {k2l2.L2 = l;};
	inline void l3(const SO3::LABEL& l) {k3l3.L2 = l;};
private:
	SU3::SO3::LABELS k1l1, k2l2, k3l3;
};

// class CSU3CGMaster is derived from CBlocks and hence the constructor of
// CBlocks is called automatically (required for the proper functioning of
// fortran SU(3) subroutines).
//
// None of the SU3CGMaster::Get methods invokes clear() on data structure which
// is provided for the output results. It is responsibility of the outer subroutine
// to make sure the structure is empty. Otherwise std::vector::push_back() is applied
// and hence the size will be ever incresing.
class CSU3CGMaster: public CBlocks {
	static double m_Zero;
	double m_dCG[MAX_K][MAX_K][MAX_K][MAX_K]; // this array is required by Fortran void wru3optimized_
    public: 
	CSU3CGMaster() {memset(m_dCG, sizeof(m_dCG), 0);}
	public:
	enum UZ6lmType {U6LM, Z6LM};
///////////////////////////// SU(3) > SO(3) //////////////////////////////////
//
// Returns <(lm1 mu1) ... ; (lm2 mu2) ... || (lm3 mu3) k3 l3>_{rho} Resulting
// data structure: vector whose elements are {[k1,l1,k2,l2], double}	NOTE:
// it is guaranteed that only non-zero SU(3) Wigner coefficients are returned
	size_t GetSO3(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir3, const SU3::SO3::LABELS& K3L3, std::vector< std::pair<K1L1K2L2, double> >& vCGs);

// Returns <(lm1 mu1) ... ; (lm2 mu2) ... || (lm3 mu3) ... l3>_{...} Resulting
// data structure: vector whose elements are {[k1,l1,k2,l2,k3],
// {rho_{0},rho_{1},...rho_{max}}}	NOTE: some of the resulting SU(3) Wigner
// coefficients may be equal to zero for a certain value of \rho_{i} and hence
// it is wrong to assume that resulting coefficients are non-zero.
	size_t GetSO3(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir3, const SO3::LABEL& L2, std::vector<std::pair<K1L1K2L2K3, std::vector<double> > >& vCGs);
	
// Returns a complete set of SU(3) Wigner coefficients <(lm1 mu1) ..., (lm2
// mu2) ... || (lm3 mu3) ...>_{...}	NOTE: some of the resulting SU(3) Wigner
// coefficients (stored in vector<double>) may be equal to zero for a certain
// value of \rho_{i} and hence it is wrong to assume that resulting
// coefficients are non-zero.
	 size_t GetSO3(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir3, std::vector<std::pair<K1L1K2L2K3L3, std::vector<double> > >& vCGs);
/////////////////////////// SU(3) > SU(2)xU(1) ///////////////////////////////
size_t GetSU2U1(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir3, std::vector<std::pair<EPS1LM1EPS2LM2EPS3LM3, std::vector<double> > >& vCGs);


// 	Get SU(3) C-G coefficients that are needed for r.m.e calculation. 
//	< ir1 eps1_{hw} LM1_{hw}; ir2 eps2 LM2 || ir3 eps3_{hw} LM3_{hw}>, where 
//		eps1_{hw} = 2lm1 + mu1
//		LM1_{hw} = mu1/2
//		eps2 = eps3_{hw} - eps1_{hw}
//		LM2 ... provided in vector<SU2::LABEL> LM2s
//		eps3_{hw} = 2lm3 + mu3
//		LM3_{hw} = mu3/2
template <typename T1, typename T2, typename T3>
size_t RMESU2U1(const T1& ir1, const T2& ir2, const T3& ir3, std::vector<SU2::LABEL>& LM2s, std::vector<std::vector<double> >& vCGs);

void Get6lm(	const UZ6lmType Type,
				const SU3::LABELS& ir1,
				const SU3::LABELS& ir2,
				const SU3::LABELS& ir,
				const SU3::LABELS& ir3,
				const SU3::LABELS& ir12,
				const SU3::LABELS& ir23, 
				std::vector<double>& U6lm);

void Get9lm(const SU3::LABELS& ir1,
			const SU3::LABELS& ir2,
			const SU3::LABELS& ir12,
			const SU3::LABELS& ir3,
			const SU3::LABELS& ir4,
			const SU3::LABELS& ir34, 
			const SU3::LABELS& ir13, 
			const SU3::LABELS& ir24, 
			const SU3::LABELS& ir, 
			std::vector<double>& su39lm);

void Get9lm(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir12,
			const SU3::LABELS& ir3, const SU3::LABELS& ir4, const SU3::LABELS& ir34, 
			const SU3::LABELS& ir13, const SU3::LABELS& ir24, const SU3::LABELS& ir, const int ntotal, double* su39lm);
void Get9lm(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir12,
			const SU3::LABELS& ir3, const SU3::LABELS& ir4, const SU3::LABELS& ir34, 
			const SU3::LABELS& ir13, const SU3::LABELS& ir24, const SU3::LABELS& ir, const int ntotal, float* su39lm);


template <typename T1, typename T2, typename T3>
void GetPhi(const T1& ir1, const T2& ir2, const T3& ir3, std::vector<double>& Phi);
};

template <typename T1, typename T2, typename T3>
void CSU3CGMaster::GetPhi(const T1& ir1, const T2& ir2, const T3& ir3, std::vector<double>& Phi)
{
	int ntotal = Phi.size();
	if (ntotal == 0) 
	{
		std::cerr << "CSU3CGMaster::GetPhi \t" << "size of Phi matrix is zero!? " << std::endl;
		exit(EXIT_FAILURE);
	}
	int zero = 0;
	int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu);

	int rho12max   = 1;
	int rho23max   = 1;
	int rho12_3max = SU3::mult(ir1, ir2, ir3);
	int rho1_23max = SU3::mult(ir1, ir2, ir3);

	double dZ6lm[ntotal];
	memset(dZ6lm, ntotal*sizeof(double), 0.0);
#ifndef AIX	
	wzu3optimized_(lm1, mu1, zero, zero, lm3, mu3, lm2, mu2, lm1, mu1, lm2, mu2, rho12max, rho12_3max, rho23max, rho1_23max, dZ6lm, ntotal);
	if (Phi.size() != rho12max*rho12_3max*rho23max*rho1_23max)
	{
		std::cerr << "CSU3CGMaster::GetPhi ... something wrong with multiplicities and size of the output vector " << std::endl;
		exit(EXIT_FAILURE);
	}
#else 
	wzu3optimized(lm1, mu1, zero, zero, lm3, mu3, lm2, mu2, lm1, mu1, lm2, mu2, rho12max, rho12_3max, rho23max, rho1_23max, dZ6lm, ntotal);
	if (Phi.size() != rho12max*rho12_3max*rho23max*rho1_23max)
	{
		std::cerr << "CSU3CGMaster::GetPhi ... something wrong with multiplicities and size of the output vector " << std::endl;
		exit(EXIT_FAILURE);
	}

#endif
	Phi.insert(Phi.begin(), dZ6lm, dZ6lm + ntotal);
}

// 	Get SU(3) C-G coefficients that are needed for r.m.e calculation. 
//	< ir1 eps1_{hw} LM1_{hw}; ir2 eps2 LM2 || ir3 eps3_{hw} LM3_{hw}>, where 
//		eps1_{hw} = 2lm1 + mu1
//		LM1_{hw} = mu1/2
//		eps2 = eps3_{hw} - eps1_{hw}
//		LM2 ... all possible values of LM2 for a given eps2 (stored in LM2s)
//		eps3_{hw} = 2lm3 + mu3
//		LM3_{hw} = mu3/2
//		
//		OUTPUT:
//		LM2s	contains values of LM2
//		vCGs is two-dimensional vector
//		Example:
//		<ir1 eps1_{hw} LM1_{hw}; ir2 eps2 LM2 || ir3 eps3_{hw} LM3_{hw}>_{rho}
//		is stored as
//		vCGs[rho][irow], where LM2 = LM2[irow]

template <typename T1, typename T2, typename T3>
size_t CSU3CGMaster::RMESU2U1(const T1& ir1, const T2& ir2, const T3& ir3, std::vector<SU2::LABEL>& LM2s, std::vector<std::vector<double> >& vCGs)
{
	int NCE1 = 9; 
	int NCE2 = 13244;
	int KIMAX1 = 3*NCE2;

	int NCW1 = 9;
	int NCW2 = 42;
	int NCW3 = 9030;
	int KIMAX2 = 3*NCW3;
	const size_t NCW22 = NCW2*NCW2;


	int I3  = 1;
	int NEC(0);
	int max_rho(0);

	int indmax(0);

	int j1ta[NCE2];
	int j2ta[NCE2];
	int iea[NCE2];
	int j1smax[NCW22];
	int j1tmax[NCW22];
	int j2smax[NCW2];
	int j2tmax[NCW2];
	int indmat[NCW22];
	
	double dewu3[KIMAX1];
	double dwu3[KIMAX2];

	int nCGs, eps2max;

    int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu);
	int eps1 = 2*lm1 + mu1;  // eps1_{hws}
	int eps3 = 2*lm3 + mu3;  // eps3_{hws}
	int eps2 = eps3 - eps1;  // since eps1 + eps2 = eps3
	int LM1 = mu1;
	int LM3 = mu3;
	int eps2_, LM1_;


	int irho, ind, minro, iesj2s, ies, j2s, j1s;
#ifndef AIX			
	xewu3_(lm1, mu1, lm2, mu2, lm3, mu3, I3, NEC, max_rho, indmax, dewu3, j1ta, j2ta, iea, NCE1, NCE2, KIMAX1);
#else
	xewu3 (lm1, mu1, lm2, mu2, lm3, mu3, I3, NEC, max_rho, indmax, dewu3, j1ta, j2ta, iea, NCE1, NCE2, KIMAX1);
#endif
	LM2s.resize(0);
//	xewu3 returns maximal multiplicity 	(ir1 x ir2) --> ir3
	vCGs.resize(max_rho);
// Order of the following loops is very important. Method uncoupling two
// tensors [constructor of class Tensor] iterates over e3 lm3 and hence we went
// them in sorted order.
#ifndef AIX			
	xwu3_(lm1, mu1, lm2, mu2, lm3, mu3, eps3, LM3, NEC, dewu3, max_rho, indmax, dwu3, j1smax, j1tmax, j2smax, j2tmax, nCGs, eps2max, indmat, NCW1, NCW2, NCW3, KIMAX2);
#else 
	xwu3 (lm1, mu1, lm2, mu2, lm3, mu3, eps3, LM3, NEC, dewu3, max_rho, indmax, dwu3, j1smax, j1tmax, j2smax, j2tmax, nCGs, eps2max, indmat, NCW1, NCW2, NCW3, KIMAX2);
#endif
// ATTENTION! I have realized that eps2 = eps1 - eps3 corresponds to ies = 0. I
// tested it for multiple cases and it seems to work. Just be aware that with
// the change of inner logic of xwu3 this may no longer be true.
	ies = 0;
	for (j2s = 0; j2s < j2smax[ies]; ++j2s)
	{
		LM2s.push_back(j2tmax[ies] - 2*j2s);
		iesj2s = ies + NCW2*j2s;
		for (j1s = 0; j1s < j1smax[iesj2s]; ++j1s)
		{
			LM1_ = j1tmax[iesj2s] - 2*j1s; 
			if (LM1 != LM1_) {
				continue;
			}
			ind = (indmat[iesj2s] - LM1)/2;
			minro = max_rho * (ind - 1);
			for (irho = 0; irho < max_rho; ++irho)
			{
				vCGs[irho].push_back(dwu3[minro + irho]); // store SU(3) wigner coefficient
			}
		}
	}
	return max_rho;
}

template <typename KEY>
struct less_first 
{
	bool operator()(const std::pair<KEY, std::vector<double> >& lhs, const std::pair<KEY, std::vector<double> >& rhs)
	{
		return lhs.first < rhs.first;
	}
};

template<typename KEY>
struct key_less
{
	bool operator()( const std::pair<KEY,std::vector<double> >& lhs, const KEY& rhs ) const { return lhs.first < rhs; };
   	bool operator()( const KEY& lhs, const std::pair<KEY,std::vector<double> >& rhs ) const { return lhs < rhs.first; };
};
#endif
