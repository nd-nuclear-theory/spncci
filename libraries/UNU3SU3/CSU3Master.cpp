#include <UNU3SU3/CSU3Master.h>
#include <SU3NCSMUtils/CTuple.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <cassert>

bool CBlocks::m_fBlocksInitiated = false;
double CSU3CGMaster::m_Zero = std::numeric_limits<float>::epsilon();

int NCE1 = 9; 
int NCE2 = 13244;
int KIMAX1 = 3*NCE2;

int NCW1 = 9;
int NCW2 = 42;
int NCW3 = 9030;
int KIMAX2 = 3*NCW3;
const size_t NCW22 = NCW2*NCW2;


// Note that the resulting structure of vCGs is sorted, that is, the order of
// [l1,l2,l3,k3,k2,k1] tuples is identical to the sorted order method
// sort(vCGs.begin(), vCGs.end(), CSU3CGMaster::CompareKeys); would produce.
size_t CSU3CGMaster::GetSO3(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir3, std::vector<std::pair<K1L1K2L2K3L3, std::vector<double> > >& vCGs)
{
	// This order of labels (coresponds to order of loops) is required to yield
	// elements of vCGs vector in sorted order.  This has the obvious advantage
	// that one need not to invoke sort(vCGs.begin(), vCGs.end()),
	// CSU3CGMaster::CompareKeys).
    int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu);
	int l1, l2, l3;
	int k1, k2, k3, max_k1, max_k2, max_k3, max_rho, irho;
	std::vector<SO3::LABEL> vl1, vl2, vl3;
	K1L1K2L2K3L3 Key;
	std::vector<double> CGrho(SU3::mult(ir1, ir2, ir3));

	int nl1 = SU3::GetAngularMomenta(ir1, vl1);   // number of L's in (lm1 mu1)
	int nl2 = SU3::GetAngularMomenta(ir2, vl2);   // number of L's in (lm2 mu2)
	int nl3 = SU3::GetAngularMomenta(ir3, vl3);   // number of L's in (lm3 mu3)
	size_t i1, i2, i3; // indices of vl1, vl2, and vl3 vectors 

	for (i1 = 0; i1 < nl1; i1++) // iterate over all l1 in (lm1 mu1)
	{
		l1 = vl1[i1];
		Key.l1(2*l1);
		for (i2 = 0; i2 < nl2; i2++) // iterate over all l2 in (lm2 mu2)
		{
			l2 =  vl2[i2];
			Key.l2(2*l2);
			for (i3 = 0; i3 < nl3; i3++)
			{
				l3 = vl3[i3];
				if (l3 < abs(l1-l2) || l3 > (l1+l2)) {
					continue;
				}
				Key.l3(2*l3);
#ifndef AIX			
				wu3r3w_(lm1, mu1, lm2, mu2, lm3, mu3, l1, l2, l3, max_rho, max_k1, max_k2, max_k3, m_dCG);
#else
				wu3r3w (lm1, mu1, lm2, mu2, lm3, mu3, l1, l2, l3, max_rho, max_k1, max_k2, max_k3, m_dCG);
#endif
			   	for (k3 = 0; k3 < max_k3; k3++) // for all values of k3
				{ 
					Key.k3(k3);
					for (k2 = 0; k2 < max_k2; k2++) // for all values of k2
			   		{
						Key.k2(k2);
						for (k1 = 0; k1 < max_k1; k1++)  // for all values of k1
			   		    {
							Key.k1(k1);
				    		for (irho = 0; irho < max_rho; irho++) // for all possible multiplicities
							{
								CGrho[irho] = m_dCG[k3][k2][k1][irho]; // store SU(3) wigner coefficient
		    					m_dCG[k3][k2][k1][irho] = 0.0; // Keep m_dCG clean TODO: test if I can remove it
	    				    }
							vCGs.push_back(std::make_pair(Key, CGrho));
			   			}
    	    		}
	    		}
			}
		}
	}
	return vCGs.size();
}

// Returns <(lm1 mu1) ... ; (lm2 mu2) ... || (lm3 mu3) ... l3>_{...}
// Resulting data structure: vector whose elements are {[l1,l2,k3,k2,k1], {rho_{0},rho_{1},...rho_{max}}}	
size_t CSU3CGMaster::GetSO3(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir3, const SO3::LABEL& L2, std::vector<std::pair<K1L1K2L2K3, std::vector<double> > >& vCGs)
{
    int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu);
	int l1, l2, l3(L2/2), k1, k2, k3;
	int max_k1, max_k2, max_k3, max_rho;
	int irho;
	std::vector<SO3::LABEL> vl1, vl2;
	K1L1K2L2K3 Key;
	std::vector<double> CGrho(SU3::mult(ir1, ir2, ir3));

//	std::cout << ir1 << " x " << ir2 << "\t" << ir3 << std::endl;
//	std::cout << "2*L0 = " << (int)L2 << std::endl;

	int nl1 = SU3::GetAngularMomenta(ir1, vl1);   // number of L's in (lm1 mu1)
	int nl2 = SU3::GetAngularMomenta(ir2, vl2);   // number of L's in (lm2 mu2)
	size_t i1, i2; // indices of vl1 and vl2 vectors 

	for (i1 = 0; i1 < nl1; i1++) // iterate over all l1 in (lm1 mu1)
	{
		l1 = vl1[i1];
		Key.l1(2*l1);
		for (i2 = 0; i2 < nl2; i2++) // iterate over all l2 in (lm2 mu2)
		{
			l2 =  vl2[i2];
			Key.l2(2*l2);
			if (l3 < abs(l1-l2) || l3 > (l1+l2)) {
				continue;
			}
#ifndef AIX			
			wu3r3w_(lm1, mu1, lm2, mu2, lm3, mu3, l1, l2, l3, max_rho, max_k1, max_k2, max_k3, m_dCG);
#else
			wu3r3w (lm1, mu1, lm2, mu2, lm3, mu3, l1, l2, l3, max_rho, max_k1, max_k2, max_k3, m_dCG);
#endif
		   	for (k3 = 0; k3 < max_k3; k3++) // for all values of k3
			{ 
				Key.k3(k3);
    			for (k2 = 0; k2 < max_k2; k2++)  // for all values of k1
		   		{
					Key.k2(k2);
    			    for (k1 = 0; k1 < max_k1; k1++) // for all values of k2
		   		    {
						Key.k1(k1);
			    		for (irho = 0; irho < max_rho; irho++) // for all possible multiplicities
						{
							CGrho[irho] = m_dCG[k3][k2][k1][irho]; // store SU(3) wigner coefficient
		    				m_dCG[k3][k2][k1][irho] = 0.0; // Keep m_dCG clean TODO: test if I can remove it
    				    }
						vCGs.push_back(std::make_pair(Key, CGrho));
		   			}
    	    	}
	    	}
		}
	}
//	std::cout << "Returning " << vCGs.size() << " CGS" << std::endl;
	return vCGs.size();
}


size_t CSU3CGMaster::GetSO3(
						const SU3::LABELS& ir1, 
						const SU3::LABELS& ir2, 
						const SU3::LABELS& ir3, const SU3::SO3::LABELS& K3L3, std::vector< std::pair<K1L1K2L2, double> >& vCGs)
{
	int Rho = ir3.rho;
	if (Rho == -1) {
		std::cerr << "You need to specify rho_{3}!" << std::endl;
		exit(EXIT_FAILURE);
	}

    int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu);
	int l1, l2, l3(K3L3.L2/2), k1, k2, k3;
	int max_k1, max_k2, max_k3, max_rho;
	int irho;
	int K3(K3L3.K);
	std::vector<SO3::LABEL> vl1, vl2;
	K1L1K2L2 Key;
	double dCGCoefficient;

	int nl1(SU3::GetAngularMomenta(ir1, vl1));   // number of L's in (lm1 mu1)
	int nl2(SU3::GetAngularMomenta(ir2, vl2));   // number of L's in (lm2 mu2)
	size_t i1, i2; // indices of vl1 and vl2 vectors 

	for (i1 = 0; i1 < nl1; i1++)
	{
		l1 = vl1[i1];
		Key.l1(2*l1);
		for (i2 = 0; i2 < nl2; i2++)
		{
			l2 = vl2[i2];
			Key.l2(2*l2);
#ifndef AIX			
			wu3r3w_(lm1, mu1, lm2, mu2, lm3, mu3, l1, l2, l3, max_rho, max_k1, max_k2, max_k3, m_dCG);
#else
			wu3r3w (lm1, mu1, lm2, mu2, lm3, mu3, l1, l2, l3, max_rho, max_k1, max_k2, max_k3, m_dCG);
#endif

	   		for (k3 = 0; k3 < max_k3; k3++) // I'm not sure if I need to clean m_dCG structure and
											// hence I iterate through all k3 even though I know
											// that I need just k3 == K 
			{									   
			    for (k2 = 0; k2 < max_k2; k2++)
				{
					Key.k2(k2);
    				for (k1 = 0; k1 < max_k1; k1++)
					{
						Key.k1(k1);
			    		for (irho = 0; irho < max_rho; irho++) // I'm not sure whether it is safe not to clean m_dCG 
															   // structure and hence I iterate through all irho even 
															   // though I know that I need just irho == Rho.
						{
    						dCGCoefficient = m_dCG[k3][k2][k1][irho]; // This is crucial!
		    				if (fabs(dCGCoefficient) > m_Zero)
    						{
								if ((Rho == irho) && (K3 == k3))
								{
									vCGs.push_back(std::make_pair(Key, dCGCoefficient));
								}
		    				    m_dCG[k3][k2][k1][irho] = 0.0; // Keep m_dCG clean TODO: test if I can remove it
		    				}
    				    }
		   			}
    	    	}
	    	}
		}
	}
	return vCGs.size();
}
// Returns <(lm1 mu1) ... ; (lm2 mu2) ... || (lm3 mu3) k3 l3>_{rho} Resulting
// data structure: vector whose elements are {[k1,l1,k2,l2], double}	NOTE:
// it is guaranteed that only non-zero SU(3) Wigner coefficients are returned
size_t CSU3CGMaster::GetSU2U1(	const SU3::LABELS& ir1, 
								const SU3::LABELS& ir2, 
								const SU3::LABELS& ir3, 
								std::vector<std::pair<EPS1LM1EPS2LM2EPS3LM3, std::vector<double> > >& vCGs)
{
	int I3  = 1;
	int NEC(0);
	int max_rho(0), indmax(0);

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
	EPS1LM1EPS2LM2EPS3LM3 Key;

	std::vector<SU3::SU2U1::LABELS> eps3LM3;
	size_t neps3twoLM2 = SU3::GetSU2U1Labels(ir3, eps3LM3);
	int eps3, LM3, eps2, eps1, LM2, LM1;
	int irho, ind, minro, iesj2s, ies, j2s, j1s;
	size_t iposition = 0;
#ifndef AIX			
	xewu3_(lm1, mu1, lm2, mu2, lm3, mu3, I3, NEC, max_rho, indmax, dewu3, j1ta, j2ta, iea, NCE1, NCE2, KIMAX1);
#else
	xewu3 (lm1, mu1, lm2, mu2, lm3, mu3, I3, NEC, max_rho, indmax, dewu3, j1ta, j2ta, iea, NCE1, NCE2, KIMAX1);
#endif
// Order of the following loops is very important. Method uncoupling two
// tensors [constructor of class Tensor] iterates over e3 lm3 and hence we went
// them in sorted order.
	for (size_t i = 0; i < neps3twoLM2; ++i)
	{
		eps3 	= eps3LM3[i].eps;
		Key.eps3(eps3);
		LM3 	= eps3LM3[i].LM;
		Key.LM3(LM3);
#ifndef AIX			
		xwu3_(lm1, mu1, lm2, mu2, lm3, mu3, eps3, LM3, NEC, dewu3, max_rho, indmax, dwu3, j1smax, j1tmax, j2smax, j2tmax, nCGs, eps2max, indmat, NCW1, NCW2, NCW3, KIMAX2);
#else 
		xwu3 (lm1, mu1, lm2, mu2, lm3, mu3, eps3, LM3, NEC, dewu3, max_rho, indmax, dwu3, j1smax, j1tmax, j2smax, j2tmax, nCGs, eps2max, indmat, NCW1, NCW2, NCW3, KIMAX2);
#endif
        for (ies = 0; ies < nCGs; ++ies)
		{
			eps2 = eps2max - 3*(nCGs - (ies + 1));
			Key.eps2 (eps2);
			eps1 = eps3 - eps2;
			Key.eps1(eps1);
			for (j2s = 0; j2s < j2smax[ies]; ++j2s)
			{
                  LM2 = j2tmax[ies] - 2*j2s;
				  Key.LM2(LM2);
                  iesj2s = ies + NCW2*j2s;
				  for (j1s = 0; j1s < j1smax[iesj2s]; ++j1s)
				  {
                      LM1 = j1tmax[iesj2s] - 2*j1s; 
					  Key.LM1(LM1);
                      ind = (indmat[iesj2s] - LM1)/2;
					  minro = max_rho * (ind - 1);

					  std::vector<double> CGrho(SU3::mult(ir1, ir2, ir3), 0.0);
					  for (irho = 0; irho < max_rho; ++irho)
					  {
						  CGrho[irho] = dwu3[minro + irho]; // store SU(3) wigner coefficient
					  }
					  vCGs.push_back(std::make_pair(Key, CGrho));
				  }
			}
		}
	}
	return vCGs.size();
}

//	This function calculates SU(3) U6 symbols.
//	NOTE: U6lm.size() == rho12*rho23*rho12_3*rho1_23
void CSU3CGMaster::Get6lm(	const CSU3CGMaster::UZ6lmType Type,
							const SU3::LABELS& ir1, 
							const SU3::LABELS& ir2, 
							const SU3::LABELS& ir, 
							const SU3::LABELS& ir3, 
							const SU3::LABELS& ir12, 
							const SU3::LABELS& ir23, 
							std::vector<double>& dUZ6lm)
{
	int rho12max   = SU3::mult(ir1, ir2, ir12);
	int rho23max   = SU3::mult(ir2, ir3, ir23);
	int rho12_3max = SU3::mult(ir12, ir3, ir);
	int rho1_23max = SU3::mult(ir1, ir23, ir);
	int ntotal = rho12max*rho23max*rho12_3max*rho1_23max;

	int lm1(ir1.lm); 
	int mu1(ir1.mu); 

	int lm2(ir2.lm); 
	int mu2(ir2.mu); 
	
	int lm3(ir3.lm); 
	int mu3(ir3.mu); 
	
	int lm(ir.lm); 
	int mu(ir.mu); 
	
	int lm12(ir12.lm); 
	int mu12(ir12.mu); 
	
	int lm23(ir23.lm); 
	int mu23(ir23.mu);

	dUZ6lm.resize(ntotal);
#ifndef AIX	
	if (Type == U6LM)
	{
		wru3optimized_(lm1, mu1, lm2, mu2, lm, mu, lm3, mu3, lm12, mu12, lm23, mu23, rho12max, rho12_3max, rho23max, rho1_23max, &dUZ6lm[0], ntotal);
	} 
	else 
	{
		wzu3optimized_(lm1, mu1, lm2, mu2, lm, mu, lm3, mu3, lm12, mu12, lm23, mu23, rho12max, rho12_3max, rho23max, rho1_23max, &dUZ6lm[0], ntotal);
	}
#else 
	if (Type == U6LM)
	{
		wru3optimized(lm1, mu1, lm2, mu2, lm, mu, lm3, mu3, lm12, mu12, lm23, mu23, rho12max, rho12_3max, rho23max, rho1_23max, &dUZ6lm[0], ntotal);
	}
	else
	{
		wzu3optimized(lm1, mu1, lm2, mu2, lm, mu, lm3, mu3, lm12, mu12, lm23, mu23, rho12max, rho12_3max, rho23max, rho1_23max, &dUZ6lm[0], ntotal);
	}
#endif
}

void CSU3CGMaster::Get9lm(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir12,
			const SU3::LABELS& ir3, const SU3::LABELS& ir4, const SU3::LABELS& ir34, 
			const SU3::LABELS& ir13, const SU3::LABELS& ir24, const SU3::LABELS& ir, 
			std::vector<double>& su39lm)
{
	int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu), lm4(ir4.lm), mu4(ir4.mu); 
	int lm12(ir12.lm), mu12(ir12.mu), lm34(ir34.lm), mu34(ir34.mu), lm13(ir13.lm), mu13(ir13.mu), lm24(ir24.lm), mu24(ir24.mu), lm(ir.lm), mu(ir.mu); 
	
	int ntotal = su39lm.size();
	if (ntotal == 0) 
	{
		std::cerr << "CSU3CGMaster::Get9lm \t" << "size of su39lm is zero!? " << std::endl;
		exit(EXIT_FAILURE);
	}
	double dsu39lm[ntotal];
#ifndef AIX	
	wu39lm_(lm1, mu1, lm2, mu2, lm12, mu12, 
			lm3, mu3, lm4, mu4, lm34, mu34, 
			lm13, mu13, lm24, mu24, lm, mu, 
			dsu39lm, ntotal);
#else 
	wu39lm(lm1, mu1, lm2, mu2, lm12, mu12, 
			lm3, mu3, lm4, mu4, lm34, mu34, 
			lm13, mu13, lm24, mu24, lm, mu, 
			dsu39lm, ntotal);
#endif
//	if su3lm.size() = n	
//	su39lm.insert(su39lm.begin(), dsu39lm, dsu39lm + ntotal);
//	su3lm.size() = ntotal + 1 after calling ... weird ...
	for (size_t i = 0; i < ntotal; ++i)
	{
		su39lm[i] = dsu39lm[i];
	}
}

void CSU3CGMaster::Get9lm(const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir12,
			const SU3::LABELS& ir3, const SU3::LABELS& ir4, const SU3::LABELS& ir34, 
			const SU3::LABELS& ir13, const SU3::LABELS& ir24, const SU3::LABELS& ir, int ntotal, double* su39lm)
{
	assert(ntotal);
//TODO: to speed up (1) do not declare and initialize local variables (2) declare method as inline
	int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu), lm4(ir4.lm), mu4(ir4.mu); 
	int lm12(ir12.lm), mu12(ir12.mu), lm34(ir34.lm), mu34(ir34.mu), lm13(ir13.lm), mu13(ir13.mu), lm24(ir24.lm), mu24(ir24.mu), lm(ir.lm), mu(ir.mu); 
#ifndef AIX	
	wu39lm_(lm1, mu1, lm2, mu2, lm12, mu12, 
			lm3, mu3, lm4, mu4, lm34, mu34, 
			lm13, mu13, lm24, mu24, lm, mu, 
			su39lm, ntotal);
#else 
	wu39lm(lm1, mu1, lm2, mu2, lm12, mu12, 
			lm3, mu3, lm4, mu4, lm34, mu34, 
			lm13, mu13, lm24, mu24, lm, mu, 
			su39lm, ntotal);
#endif
}

void CSU3CGMaster::Get9lm(	const SU3::LABELS& ir1, const SU3::LABELS& ir2, const SU3::LABELS& ir12,
							const SU3::LABELS& ir3, const SU3::LABELS& ir4, const SU3::LABELS& ir34, 
							const SU3::LABELS& ir13, const SU3::LABELS& ir24, const SU3::LABELS& ir, int ntotal, float* su39lm)
{
	assert(ntotal);
	double *tmp = new double[ntotal];
//TODO: to speed up (1) do not declare and initialize local variables (2) declare method as inline
	int lm1(ir1.lm), mu1(ir1.mu), lm2(ir2.lm), mu2(ir2.mu), lm3(ir3.lm), mu3(ir3.mu), lm4(ir4.lm), mu4(ir4.mu); 
	int lm12(ir12.lm), mu12(ir12.mu), lm34(ir34.lm), mu34(ir34.mu), lm13(ir13.lm), mu13(ir13.mu), lm24(ir24.lm), mu24(ir24.mu), lm(ir.lm), mu(ir.mu); 
#ifndef AIX	
	wu39lm_(lm1, mu1, lm2, mu2, lm12, mu12, lm3, mu3, lm4, mu4, lm34, mu34, lm13, mu13, lm24, mu24, lm, mu, tmp, ntotal);
#else 
	wu39lm(lm1, mu1, lm2, mu2, lm12, mu12, lm3, mu3, lm4, mu4, lm34, mu34, lm13, mu13, lm24, mu24, lm, mu, tmp, ntotal);
#endif
	std::copy(tmp, tmp + ntotal, su39lm);
	delete []tmp;
}
