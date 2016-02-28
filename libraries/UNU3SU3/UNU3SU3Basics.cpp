#include <UNU3SU3/UNU3SU3Basics.h>
size_t SU3xSU2::Couple(const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, std::vector<SU3xSU2::LABELS>& ResultingIrreps)
{
	std::vector<SU3::LABELS> vIR3;
	std::vector<SU2::LABEL> vS;

	SU3::Couple(ir1, ir2, vIR3);
	SU2::Couple(ir1.S2, ir2.S2, vS); 

	for (size_t iS = 0; iS < vS.size(); ++iS)
	{
		for (size_t i = 0; i < vIR3.size(); ++i)
		{
			ResultingIrreps.push_back(SU3xSU2::LABELS(SU3::LABELS(vIR3[i]), SU2::LABEL(vS[iS])));
		}
	}
	return ResultingIrreps.size();
}

size_t SU3::GetSO3Labels(const SU3::LABELS& RhoLmMu, std::vector<SU3::SO3::LABELS>& KL)
{
	int lmax = (RhoLmMu.lm + RhoLmMu.mu);
	int kmax;
	int k;

	for (int l = 0; l <= lmax; l++)
	{
		kmax = SU3::kmax(RhoLmMu, l);
		if (kmax) {
			for (k = 0; k < kmax; ++k)
			{	
				KL.push_back(SU3::SO3::LABELS(k, 2*l));
			}
		}
	}
	return KL.size();
}

size_t SU3::GetAngularMomenta(const SU3::LABELS& RhoLmMu, std::vector<_SO3::LABEL>& Ls)
{
	int lmax = (RhoLmMu.lm + RhoLmMu.mu);
	int kmax;
	for (int l = 0; l <= lmax; ++l)
	{
		if (SU3::kmax(RhoLmMu, l))
		{
			Ls.push_back(_SO3::LABEL(l));
		}
	}
	return Ls.size();
}

/*
int SU3::mult (const SU3::LABELS& RhoLmMu1, const SU3::LABELS& RhoLmMu2, const SU3::LABELS& RhoLmMu3)
{
	int lm1(RhoLmMu1.lm), mu1(RhoLmMu1.mu), lm2(RhoLmMu2.lm), mu2(RhoLmMu2.mu), lm3(RhoLmMu3.lm), mu3(RhoLmMu3.mu);
    int result = 0;

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
    return result;
}
*/

size_t SU3::GetCanonicalBasis(const SU3::LABELS RhoLmMu, std::vector<SU3::CANONICAL>& SU3CanonicalBasis, bool bSortedOrder)
{
	std::vector<SU3::SU2U1::LABELS> epsLam;
	std::vector<U1::LABEL> M_LMs;

	size_t nEpsLM = SU3::GetSU2U1Labels(RhoLmMu, epsLam);
	for (size_t i = 0; i < nEpsLM; ++i)
	{
		size_t nM_LMs = SU2::GetU1Labels(epsLam[i].LM, M_LMs);
		for (size_t j = 0; j < nM_LMs; ++j)
		{
			SU3CanonicalBasis.push_back(SU3::CANONICAL(epsLam[i], M_LMs[j]));
		}
	}
	if (bSortedOrder) {
		std::sort(SU3CanonicalBasis.begin(), SU3CanonicalBasis.end(), std::less<SU3::CANONICAL>());
	}
	return SU3CanonicalBasis.size();
}

size_t	SU3xSU2::GetPhysicalJBasis(const SU3xSU2::LABELS& IR, std::vector<SU3xSU2::PHYSICALJ>& SU3xSU2PhysicalJBasis)
{
	SU3xSU2PhysicalJBasis.reserve(SU3::dim(IR)*SU2::dim(IR.S2)); 
	int Lmax = IR.lm + IR.mu;
	for (int L = 0; L <= Lmax; ++L)
	{
		int kmax = SU3::kmax(IR, L);
		if (!kmax)
		{
			continue;
		}
		int LL = 2*L;
		for (int k = 0; k < kmax; ++k)
		{
			for (int JJ = abs(LL - IR.S2); JJ <= LL + IR.S2; JJ += 2)
			{
				for (int MM = -JJ; MM <= JJ; MM += 2)
				{
					SU3xSU2PhysicalJBasis.push_back(SU3xSU2::PHYSICALJ(k, LL, JJ, MM));
				}
			}
		}
	}
	return SU3xSU2PhysicalJBasis.size();
}

size_t SU3xSU2::GetCanonicalBasis(const SU3xSU2::LABELS& RhoLmMuS2, std::vector<SU3xSU2::CANONICAL>& SU3xSU2CanonicalBasis, bool bSortedOrder)
{
	std::vector<SU3::CANONICAL> SU3Basis;
	std::vector<SU2::CANONICAL> SU2Basis;

	size_t nEpsLMM_LM = SU3::GetCanonicalBasis(RhoLmMuS2, SU3Basis, bSortedOrder);
	size_t nSigma = SU2::GetU1Labels(RhoLmMuS2.S2, SU2Basis);

	SU3xSU2CanonicalBasis.resize(0);
	SU3xSU2CanonicalBasis.reserve(nEpsLMM_LM*nSigma);

	for (size_t i = 0; i < nEpsLMM_LM; ++i)
	{
		for (size_t j = 0; j < nSigma; ++j)
		{
			SU3xSU2CanonicalBasis.push_back(SU3xSU2::CANONICAL(SU3Basis[i], SU2Basis[j]));
		}
	}

	if (bSortedOrder) {
		std::sort(SU3xSU2CanonicalBasis.begin(), SU3xSU2CanonicalBasis.end(), std::less<SU3xSU2::CANONICAL>());
	}
	return SU3xSU2CanonicalBasis.size();
}
