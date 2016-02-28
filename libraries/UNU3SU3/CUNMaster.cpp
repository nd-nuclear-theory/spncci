#include <UNU3SU3/CUNMaster.h>
#include <cmath>
#include <iostream>
/* 
 * This function returns dimension of a U(N) irrep defined by Young shape
 * [f_{1},f_{2},...,f_{N}]=f[0...N-1];
 *
 * WARNING: This function has problems with rounding errors
 * For instance it used to report dimension that was less by one 
 * than the one calculated by summation of allowed SU(3) irreps.
 * That is why I added floor(exp(...)+0.5) -> bare in mind that this is just a fix!
 */
unsigned long UN::dim(const UN::LABELS& f)
{
    const size_t N = f.size();
    long double ldNominator   = 1.0;
    long double ldDenominator = 1.0;
    
    for (unsigned long l = 1; l < N; l++)
    {
		for (unsigned long k = 0; k < l; k++)
		{
	    	ldNominator   += log((long double)(f[k]-f[l]+l-k));
	    	ldDenominator += log((long double)(l-k));
		}
    }
    return floor(exp((ldNominator - ldDenominator))+0.5000000);
}

/*
U3::U3(const unsigned nmax)
{
    U3::LABELS U3Labels; // [0] -> nz ; [1] -> nx ; [2] -> ny

	for (unsigned n = 0 ; n <= nmax; n++)
	{
		U3::SPS ShellSPS; // set of sps (Nz, Nx, Ny) for harmonic oscillator shell n

		size_t N = (n+1)*(n+2)/2;
		std::vector<long double> CzxCoeff(N);
		std::vector<long double> CzyCoeff(N);
		std::vector<long double> CxyCoeff(N);
		size_t i = 0;
	    for (long k = 0; k <= n; k++)
    	{
			U3Labels[NZ] = n - k;
			for (long nx = k; nx >= 0; nx--)
			{
	    		U3Labels[NX] = nx;
		    	U3Labels[NY] = (n - U3Labels[NZ] - nx);
			    ShellSPS.push_back(U3Labels);

				CzxCoeff[i] = sqrt((U3Labels[NZ]+1.0)*nx);
				CzyCoeff[i] = sqrt((U3Labels[NZ]+1.0)*U3Labels[NY]);
				CxyCoeff[i] = sqrt((nx+1.0)*U3Labels[NY]);
				i++;
			}
		}
		m_AllShellSPS.push_back(ShellSPS);
		m_AllShells_CzxCoeff.push_back(CzxCoeff);
		m_AllShells_CzyCoeff.push_back(CzyCoeff);
		m_AllShells_CxyCoeff.push_back(CxyCoeff);
	}
}
*/

/*
 * INPUT:
 * (1) vWeights: 
 * A weight vector of a basis state of U(N) 
 * (2) U3:SPS 
 * a list of single particle states of HO in speedometer order for a given shell 
 * A distribution of HO quanta in (z, x, y) directions associated with each of
 * N levels of U(N)  
 *
 * TASK: 
 * calculate labels of U(3) irrep using equation (4) of CPC.
 *
 * OUTPUT: 
 * U(3) labels [N_{z}, N_{x}, N_{y}]
 */
void CUNMaster::Weight2U3Label(
	const UN::BASIS_STATE_WEIGHT_VECTOR& vWeights,
	const U3::SPS& ShellSPS, 
	U3::LABELS& vU3Labels)
{
    size_t N = ShellSPS.size();
    for (size_t i = 0; i < N; i++)
    {
		vU3Labels[U3::NZ] += vWeights[i]*ShellSPS[i][U3::NZ];
		vU3Labels[U3::NX] += vWeights[i]*ShellSPS[i][U3::NX];
		vU3Labels[U3::NY] += vWeights[i]*ShellSPS[i][U3::NY];
    }
}

/*
 * INPUT:
 * (1) vGelfandParentRow:
 * Young shape [f] = [f_{1},f_{2},...,f_{N}] that labels an U(N) irrep
 * (2) vvShellStructure:
 * a distribution of HO quanta in (z, x, y) directions associated with each of
 * N levels of U(N)  
 * (3) vWeights:
 * A weight vector of a basis state of U(N) irrep. This parameter is used in a 
 * recursive execution. An empty N-dimensional vector must be provided on input.
 *
 * TASK: 
 * evaluate all allowed U(3) patterns [N_{z}, N_{x}, N_{y}] of U(N) basis
 * states and their multiplicities. Algorithm described in Computer Physics
 * Communications 56 (1989) 279-290
 *
 * OUTPUT: 
 * Map of [N_{z}, N_{x}, N_{y}] labels associated with a corresponding multiplicity.
 * This map can be used to calculate allowed SU(3)xSU(2) irreps
 */
void CUNMaster::GenerateU3Labels(
	const UN::LABELS& vGelfandParentRow, 
	const U3::SPS& ShellSPS, 
	UN::BASIS_STATE_WEIGHT_VECTOR& vWeights,
    UN::U3MULT_LIST& mU3LabelsOccurance)
{
    size_t N = vGelfandParentRow.size()-1;
    
    std::vector<UN::LABELS> vvAllowedLabels;
    UN::LABELS vGelfandRow(N);
    std::vector<unsigned> vElemsPerChange(N, 1); //Fill with 1 cause vElemsPerChange[0] = 1;
    std::vector<size_t> vNLabels(N);

    // uSumGelfandParentRow is equal to the sum of a parent Gelfand row 
    unsigned uSumGelfandParentRow = 0;
    for (size_t i = 0; i < N+1; uSumGelfandParentRow += vGelfandParentRow[i++]);

    // evaluate all allowed Gelfand patterns based on a parent Gelfand row
    // (vGelfandParentRow) and store them in vvAllowedLabels
    // iNAllowedCombinations is equal to the number of allowed Gelfand patterns
    unsigned iNAllowedCombinations = 1;
    for (size_t i = 0; i < N; i++)
    {
		UN::LABELS vLabels;
		unsigned uLabelMin = std::min(vGelfandParentRow[i+1], vGelfandParentRow[i]);
		unsigned uLabelMax = std::max(vGelfandParentRow[i+1], vGelfandParentRow[i]);
		for (unsigned uLabel = uLabelMin; uLabel <= uLabelMax; vLabels.push_back(uLabel++));
		vvAllowedLabels.push_back(vLabels);
		vNLabels[i] = vLabels.size();
		iNAllowedCombinations *= vNLabels[i];
    }

    for (size_t i = 1; i < N; i++)
    {
		vElemsPerChange[i] = vElemsPerChange[i-1]*vNLabels[i-1]; //if i == 0 => vElemsPerChange[0]=1 since constructor
    }

    for (unsigned index = 0; index < iNAllowedCombinations; index++)
    {
		unsigned uSumGelfandRow = 0;
		for (size_t i = 0; i < N; i++)
		{
		    size_t iElement = (index / vElemsPerChange[i]) % vNLabels[i];
		    vGelfandRow[i]  = vvAllowedLabels[i][iElement];
		    uSumGelfandRow += vGelfandRow[i];
		}

	 	vWeights[N] = uSumGelfandParentRow - uSumGelfandRow;
		if (N > 1) {// this condition is due to n = 0 case when N = 0 and hence one wants to keep vWeights[0] = uSumGelfandRow;
			GenerateU3Labels(vGelfandRow, ShellSPS, vWeights, mU3LabelsOccurance);
		} else {
			if (N == 1) {
				vWeights[0] = uSumGelfandRow;
			}
			U3::LABELS vU3Labels((const unsigned char*)"\0\0\0");
			Weight2U3Label(vWeights, ShellSPS, vU3Labels);
			mU3LabelsOccurance[vU3Labels] += 1;
		}
    }
}

/*
 * INPUT:
 * (1) vGelfandParentRow:
 * Young shape [f] = [f_{1},f_{2},...,f_{N}] that labels an U(N) irrep
 * (2) ShellSPS:
 * a distribution of HO quanta in (z, x, y) directions associated with each of
 * N levels of U(N)  
 * (3) vWeights:
 * A weight vector of a basis state of U(N) irrep. This parameter is used in a 
 * recursive execution. An empty N-dimensional vector must be provided on input.
 *
 * TASK: 
 * Evaluate all basis state of U(N) irrep |[f] N_{x}; N_{y}; N_{z}> and calculate how many times they occur (degeneracy)
 * If N_{z}>=N_{x}>=N_{y} store bitset representation of |[f] N_{x}; N_{y}; N_{z}> in mHWSbitsets
 *
 * OUTPUT: 
 * 1) Map of [N_{z}, N_{x}, N_{y}] labels associated with a corresponding degeneracy. 
 * This map is used to calculate allowed SU(3)xSU(2) irreps.
 * 2) Map of [N_{z}, N_{x}, N_{y}] labels such that N_{z} >= N_{x} >= N_{y} and
 * associated vector of U(N) basis states represented as a pair of bitsets.
 * These data structures are used to calculate (1) multiplicities \alpha of SU(3) irreps in U(N)) and (2) HWS of SU(3) irreps
 */
void CUNMaster::GenerateU3LabelsUNWeights(
	const UN::LABELS& vGelfandParentRow,
	const U3::SPS& ShellSPS,
	UN::BASIS_STATE_WEIGHT_VECTOR& vWeights,
	UN::U3MULT_LIST& mU3LabelsOccurance,
	std::map<U3::LABELS, std::set<UN::BASIS_STATE_WEIGHT_VECTOR> >&  mUNweights)
{
    size_t N = vGelfandParentRow.size()-1;
    
    std::vector<UN::LABELS> vvAllowedLabels;
    UN::LABELS vGelfandRow(N);
    std::vector<unsigned> vElemsPerChange(N, 1); //Fill with 1 cause vElemsPerChange[0] = 1;
    std::vector<size_t> vNLabels(N);

    // uSumGelfandParentRow is equal to the sum of a parent Gelfand row 
    unsigned uSumGelfandParentRow = 0;
    for (size_t i = 0; i < N+1; uSumGelfandParentRow += vGelfandParentRow[i++]);

    // evaluate all allowed Gelfand patterns based on a parent Gelfand row
    // (vGelfandParentRow) and store them in vvAllowedLabels
    // iNAllowedCombinations is equal to the number of allowed Gelfand patterns
    unsigned iNAllowedCombinations = 1;
    for (size_t i = 0; i < N; i++)
    {
		UN::LABELS vLabels;
		unsigned uLabelMin = std::min(vGelfandParentRow[i+1], vGelfandParentRow[i]);
		unsigned uLabelMax = std::max(vGelfandParentRow[i+1], vGelfandParentRow[i]);
		for (unsigned uLabel = uLabelMin; uLabel <= uLabelMax; vLabels.push_back(uLabel++));
		vvAllowedLabels.push_back(vLabels);
		vNLabels[i] = vLabels.size();
		iNAllowedCombinations *= vNLabels[i];
    }

    for (size_t i = 1; i < N; i++)
    {
		vElemsPerChange[i] = vElemsPerChange[i-1]*vNLabels[i-1]; //if i == 0 => vElemsPerChange[0]=1 since constructor
    }

    for (unsigned index = 0; index < iNAllowedCombinations; index++)
    {
		unsigned uSumGelfandRow = 0;
		for (size_t i = 0; i < N; i++)
		{
		    size_t iElement = (index / vElemsPerChange[i]) % vNLabels[i];
		    vGelfandRow[i]  = vvAllowedLabels[i][iElement];
		    uSumGelfandRow += vGelfandRow[i];
		}

	 	vWeights[N] = uSumGelfandParentRow - uSumGelfandRow;
		if (N > 1) {
			GenerateU3LabelsUNWeights(vGelfandRow, ShellSPS, vWeights, mU3LabelsOccurance, mUNweights);
		} else { // Another U(N) basis state has been generated. It is stored in vWeights.
			if (N == 1) { // this condition is due to n = 0 case when N = 0 and hence one wants to keep vWeights[0] = uSumGelfandRow;
				vWeights[0] = uSumGelfandRow;
			}
			U3::LABELS vU3Labels((const unsigned char*)"\0\0\0");
			Weight2U3Label(vWeights, ShellSPS, vU3Labels); // Calculate U(3) labels from vWeights
			mU3LabelsOccurance[vU3Labels] += 1; // Increment the number of U(N) basis states that carry given U(3) labels.
			std::map<U3::LABELS, std::set<UN::BASIS_STATE_WEIGHT_VECTOR> >::iterator it = mUNweights.find(vU3Labels);
			if (it == mUNweights.end())  // Is given U(3) label already stored map ?
			{
				std::set<UN::BASIS_STATE_WEIGHT_VECTOR> Basis;	//	create a new set of U(N) irrep basis states weights
				Basis.insert(vWeights);					//	and store the current weight vector in it.
				mUNweights.insert(make_pair(vU3Labels, Basis));
			} else {
				// a single weight vector can appear multiple times and
				// hence we first have to see whether a given weight is
				// already stored
				it->second.insert(vWeights); // it->second == set<UN::BASIS_STATE_WEIGHT_VECTOR>
			}
		}
    }
}

unsigned CUNMaster::GetMultiplicity(const U3::LABELS vU3Labels, const UN::U3MULT_LIST& mU3IRs)
{
    UN::U3MULT_LIST::const_iterator cit;
    unsigned f1 = vU3Labels[0], f2 = vU3Labels[1], f3 = vU3Labels[2];
    unsigned usResult;
    U3::LABELS vU3LabelsTmp;

    cit = mU3IRs.find(vU3Labels);
    if (cit == mU3IRs.end())
    {
		std::cout << "Error: U3 irrep [" << f1 << ", " << f2 << ", " << f3 << "] is not present in list of allowed U3 irreps";
		std::cout << std::endl;
		exit(0);
    }
    usResult = cit->second;

    vU3LabelsTmp[0] = f1 + 1; vU3LabelsTmp[1] = f2 + 1; vU3LabelsTmp[2] = f3 - 2;
    cit = mU3IRs.find(vU3LabelsTmp);
    usResult += (cit == mU3IRs.end()) ? 0 : cit->second;

    vU3LabelsTmp[0] = f1 + 2; vU3LabelsTmp[1] = f2 - 1; vU3LabelsTmp[2] = f3 - 1;
    cit = mU3IRs.find(vU3LabelsTmp);
    usResult += (cit == mU3IRs.end()) ? 0 : cit->second;

    vU3LabelsTmp[0] = f1 + 2; vU3LabelsTmp[1] = f2; vU3LabelsTmp[2] = f3 - 2;
    cit = mU3IRs.find(vU3LabelsTmp);
    usResult -= (cit == mU3IRs.end()) ? 0 : cit->second;

    vU3LabelsTmp[0] = f1 + 1; vU3LabelsTmp[1] = f2 - 1; vU3LabelsTmp[2] = f3;
    cit = mU3IRs.find(vU3LabelsTmp);
    usResult -= (cit == mU3IRs.end()) ? 0 : cit->second;

    vU3LabelsTmp[0] = f1; vU3LabelsTmp[1] = f2 + 1; vU3LabelsTmp[2] = f3 - 1;
    cit = mU3IRs.find(vU3LabelsTmp);
    usResult -= (cit == mU3IRs.end()) ? 0 : cit->second;

    return usResult;
}

/*
 * INPUT:
 * (1) n:
 * harmonic oscillator shell
 * (2) A:
 * number of fermions
 *
 * TASK:
 * Use algorithm described in Computer Physics Communications 56 (1989) 279-290
 * to obtain a set of allowed SU(3)xSU(2) irreps for a system of A fermions in n-th HO shell
 *
 * OUTPUT:
 * std::vector of SU(3)xSU(2) irrep labels
 */
void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS, std::vector<UN::SU3_VEC>& SpinSU3MULTLabels)
{
	std::cout << "OBSOLETE! I would rather see you to use void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<std::pair<unsigned char, UN::SU3_VEC> >& SpinSU3MULTLabels)" << std::endl;
    const unsigned N = (n+1)*(n+2)/2;
    if (2*N < A) {
		std::cout << "Error: only " << 2*N << " states available for " << A << " particles!" << std::endl;
		exit(-1);
    }
    unsigned S2min = !(A%2) ? 0 : 1; // if A = 2k => S2min = 0; S2min = 1 otherwise
    unsigned S2max = A;
	unsigned usMult;

    for (unsigned S2 = S2min; S2 <= S2max; S2+=2)
    {
		UN::SU3_VEC SU3MULTList; // vector of {mult, lm, mu} ... for a given spin S2

		UN::LABELS U2Labels(2);
		U2Labels[0] = (A + S2)/2; // Note that this is also equal to the number of fermions with spin up
		U2Labels[1] = (A - S2)/2;

		if (U2Labels[0] > N) {
			continue;
		}

		UN::LABELS UNLabels(N, 0);
		for (size_t i = 0; i < U2Labels[0]; i++)
		{
		    UNLabels[i] = (i < U2Labels[1]) ? 2 : 1;
		}
//	Iterate over U(3) labels and calculate multiplicity Mult (corresponds to
//	alpha in [f] alpha (lm mu)S).  if Mult > 0 => add a given SU(3)xSU(2) irrep into SU3xSU2Labels
    	UN::U3MULT_LIST mU3_mult;
		
		GenerateU3Labels(UNLabels, ShellSPS, mU3_mult);

		UN::U3MULT_LIST::const_iterator citU3 = mU3_mult.begin();
		UN::U3MULT_LIST::const_iterator LAST_U3 = mU3_mult.end();
//	iterate over all [N_{z}, N_{x}, N_{y}] labels spanning U(N) irrep 
       	for(; citU3 != LAST_U3; citU3++)
    	{
    	    U3::LABELS U3Labels(citU3->first);
   	    	if (U3Labels[U3::NZ] >= U3Labels[U3::NX] && U3Labels[U3::NX] >= U3Labels[U3::NY]) // this could be U(3) irrep ... let's check it
    	    {
    			usMult = GetMultiplicity(U3Labels, mU3_mult);
//	Calculate multiplicity alpha in irrep of U(N)
	    		if (usMult) 
				{
				    SU3MULTList.push_back(UN::SU3(usMult, SU3::LM(U3Labels), SU3::MU(U3Labels)));
    			}
    	    }
    	}
		SpinSU3MULTLabels.push_back(SU3MULTList);
    }
}

void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<UN::SU3_VEC>& SpinSU3Labels)
{
	std::cout << "OBSOLETE! I would rather see you to use void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<std::pair<unsigned char, UN::SU3_VEC> >& SpinSU3MULTLabels)" << std::endl;
	U3::SPS ShellSPS; // set of sps (Nz, Nx, Ny) for harmonic oscillator shell n
	size_t N = (n + 1) * (n + 2) / 2;
    U3::LABELS U3Labels; // [0] -> nz ; [1] -> nx ; [2] -> ny

    for (long k = 0; k <= n; k++)
   	{
		U3Labels[U3::NZ] = n - k;
		for (long nx = k; nx >= 0; nx--)
		{
    		U3Labels[U3::NX] = nx;
	    	U3Labels[U3::NY] = (n - U3Labels[U3::NZ] - nx);
		    ShellSPS.push_back(U3Labels);
		}
	}
	GetAllowedSU3xSU2Irreps(n, A, ShellSPS, SpinSU3Labels);
}

void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS,  std::vector<std::pair<SU2::LABEL, UN::SU3_VEC> >& SpinSU3MULTLabels)
{
    const unsigned N = (n+1)*(n+2)/2;
    if (2*N < A) {
		std::cout << "Error: only " << 2*N << " states available for " << A << " particles!" << std::endl;
		exit(-1);
    }
    SU2::LABEL S2min = !(A%2) ? 0 : 1; // if A = 2k => S2min = 0; S2min = 1 otherwise
    SU2::LABEL S2max = A;
	unsigned usMult;

    for (SU2::LABEL S2 = S2min; S2 <= S2max; S2+=2)
    {
		UN::SU3_VEC SU3MULTList; // vector of {mult, lm, mu} ... for a given spin S2

		UN::LABELS U2Labels(2);
		U2Labels[0] = (A + S2)/2; // Note that this is also equal to the number of fermions with spin up
		U2Labels[1] = (A - S2)/2;

		if (U2Labels[0] > N) {
			continue;
		}

		UN::LABELS UNLabels(N, 0);
		for (size_t i = 0; i < U2Labels[0]; i++)
		{
		    UNLabels[i] = (i < U2Labels[1]) ? 2 : 1;
		}
//	Iterate over U(3) labels and calculate multiplicity Mult (corresponds to
//	alpha in [f] alpha (lm mu)S).  if Mult > 0 => add a given SU(3)xSU(2) irrep into SU3xSU2Labels
    	UN::U3MULT_LIST mU3_mult;
		
		GenerateU3Labels(UNLabels, ShellSPS, mU3_mult);

		UN::U3MULT_LIST::const_iterator citU3 = mU3_mult.begin();
		UN::U3MULT_LIST::const_iterator LAST_U3 = mU3_mult.end();
//	iterate over all [N_{z}, N_{x}, N_{y}] labels spanning U(N) irrep 
       	for(; citU3 != LAST_U3; citU3++)
    	{
    	    U3::LABELS U3Labels(citU3->first);
   	    	if (U3Labels[U3::NZ] >= U3Labels[U3::NX] && U3Labels[U3::NX] >= U3Labels[U3::NY]) // this could be U(3) irrep ... let's check it
    	    {
    			usMult = GetMultiplicity(U3Labels, mU3_mult);
//	Calculate multiplicity alpha in irrep of U(N)
	    		if (usMult) 
				{
				    SU3MULTList.push_back(UN::SU3(usMult, SU3::LM(U3Labels), SU3::MU(U3Labels)));
    			}
    	    }
    	}
		if (!SU3MULTList.empty())
		{
			SpinSU3MULTLabels.push_back(make_pair(S2, SU3MULTList));
		}
    }
}

void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<std::pair<SU2::LABEL, UN::SU3_VEC> >& SpinSU3MULTLabels)
{
	U3::SPS ShellSPS; // set of sps (Nz, Nx, Ny) for harmonic oscillator shell n
	size_t N = (n + 1) * (n + 2) / 2;
    U3::LABELS U3Labels; // [0] -> nz ; [1] -> nx ; [2] -> ny

    for (long k = 0; k <= n; k++)
   	{
		U3Labels[U3::NZ] = n - k;
		for (long nx = k; nx >= 0; nx--)
		{
    		U3Labels[U3::NX] = nx;
	    	U3Labels[U3::NY] = (n - U3Labels[U3::NZ] - nx);
		    ShellSPS.push_back(U3Labels);
		}
	}
	GetAllowedSU3xSU2Irreps(n, A, ShellSPS, SpinSU3MULTLabels);
}

void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS, UN::SU3xSU2_VEC& AllowedIrreps)
{
    const unsigned N = (n+1)*(n+2)/2;
    if (2*N < A) {
		std::cout << "Error: only " << 2*N << " states available for " << A << " particles!" << std::endl;
		exit(-1);
    }
    SU2::LABEL S2min = !(A%2) ? 0 : 1; // if A = 2k => S2min = 0; S2min = 1 otherwise
    SU2::LABEL S2max = A;
	unsigned usMult;

    for (SU2::LABEL S2 = S2min; S2 <= S2max; S2+=2)
    {
		UN::LABELS U2Labels(2);
		U2Labels[0] = (A + S2)/2; // Note that this is also equal to the number of fermions with spin up
		U2Labels[1] = (A - S2)/2;

		if (U2Labels[0] > N) {
			continue;
		}

		UN::LABELS UNLabels(N, 0);
		for (size_t i = 0; i < U2Labels[0]; i++)
		{
		    UNLabels[i] = (i < U2Labels[1]) ? 2 : 1;
		}
//	Iterate over U(3) labels and calculate multiplicity Mult (corresponds to
//	alpha in [f] alpha (lm mu)S).  if Mult > 0 => add a given SU(3)xSU(2) irrep into SU3xSU2Labels
    	UN::U3MULT_LIST mU3_mult;
		
		GenerateU3Labels(UNLabels, ShellSPS, mU3_mult);

		UN::U3MULT_LIST::const_iterator citU3 = mU3_mult.begin();
		UN::U3MULT_LIST::const_iterator LAST_U3 = mU3_mult.end();
//	iterate over all [N_{z}, N_{x}, N_{y}] labels spanning U(N) irrep 
       	for(; citU3 != LAST_U3; citU3++)
    	{
    	    U3::LABELS U3Labels(citU3->first);
   	    	if (U3Labels[U3::NZ] >= U3Labels[U3::NX] && U3Labels[U3::NX] >= U3Labels[U3::NY]) // this could be U(3) irrep ... let's check it
    	    {
    			usMult = GetMultiplicity(U3Labels, mU3_mult);
//	Calculate multiplicity alpha in irrep of U(N)
	    		if (usMult) 
				{
				    AllowedIrreps.push_back(UN::SU3xSU2(usMult, SU3::LM(U3Labels), SU3::MU(U3Labels), S2));
    			}
    	    }
    	}
    }
}

void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, UN::SU3xSU2_VEC& AllowedIrreps)
{
	U3::SPS ShellSPS; // set of sps (Nz, Nx, Ny) for harmonic oscillator shell n
	size_t N = (n + 1) * (n + 2) / 2;
    U3::LABELS U3Labels; // [0] -> nz ; [1] -> nx ; [2] -> ny

    for (long k = 0; k <= n; k++)
   	{
		U3Labels[U3::NZ] = n - k;
		for (long nx = k; nx >= 0; nx--)
		{
    		U3Labels[U3::NX] = nx;
	    	U3Labels[U3::NY] = (n - U3Labels[U3::NZ] - nx);
		    ShellSPS.push_back(U3Labels);
		}
	}
	GetAllowedSU3xSU2Irreps(n, A, ShellSPS, AllowedIrreps);
}

void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, const U3::SPS& ShellSPS, std::vector<SU3xSU2::LABELS>& AllowedIrreps)
{
    const unsigned N = (n+1)*(n+2)/2;
    if (2*N < A) {
		std::cout << "Error: only " << 2*N << " states available for " << A << " particles!" << std::endl;
		exit(-1);
    }
    SU2::LABEL S2min = !(A%2) ? 0 : 1; // if A = 2k => S2min = 0; S2min = 1 otherwise
    SU2::LABEL S2max = A;
	unsigned usMult;

    for (SU2::LABEL S2 = S2min; S2 <= S2max; S2+=2)
    {
		UN::LABELS U2Labels(2);
		U2Labels[0] = (A + S2)/2; // Note that this is also equal to the number of fermions with spin up
		U2Labels[1] = (A - S2)/2;

		if (U2Labels[0] > N) {
			continue;
		}

		UN::LABELS UNLabels(N, 0);
		for (size_t i = 0; i < U2Labels[0]; i++)
		{
		    UNLabels[i] = (i < U2Labels[1]) ? 2 : 1;
		}
//	Iterate over U(3) labels and calculate multiplicity Mult (corresponds to
//	alpha in [f] alpha (lm mu)S).  if Mult > 0 => add a given SU(3)xSU(2) irrep into SU3xSU2Labels
    	UN::U3MULT_LIST mU3_mult;
		
		GenerateU3Labels(UNLabels, ShellSPS, mU3_mult);

		UN::U3MULT_LIST::const_iterator citU3 = mU3_mult.begin();
		UN::U3MULT_LIST::const_iterator LAST_U3 = mU3_mult.end();
//	iterate over all [N_{z}, N_{x}, N_{y}] labels spanning U(N) irrep 
       	for(; citU3 != LAST_U3; citU3++)
    	{
    	    U3::LABELS U3Labels(citU3->first);
   	    	if (U3Labels[U3::NZ] >= U3Labels[U3::NX] && U3Labels[U3::NX] >= U3Labels[U3::NY]) // this could be U(3) irrep ... let's check it
    	    {
    			usMult = GetMultiplicity(U3Labels, mU3_mult);
//	Calculate multiplicity alpha in irrep of U(N)
	    		if (usMult) 
				{
				    AllowedIrreps.push_back(SU3xSU2::LABELS(SU3::LM(U3Labels), SU3::MU(U3Labels), S2));
    			}
    	    }
    	}
    }
}


void CUNMaster::GetAllowedSU3xSU2Irreps(const unsigned n, const unsigned A, std::vector<SU3xSU2::LABELS>& AllowedIrreps)
{
	U3::SPS ShellSPS; // set of sps (Nz, Nx, Ny) for harmonic oscillator shell n
	size_t N = (n + 1) * (n + 2) / 2;
    U3::LABELS U3Labels; // [0] -> nz ; [1] -> nx ; [2] -> ny

    for (long k = 0; k <= n; k++)
   	{
		U3Labels[U3::NZ] = n - k;
		for (long nx = k; nx >= 0; nx--)
		{
    		U3Labels[U3::NX] = nx;
	    	U3Labels[U3::NY] = (n - U3Labels[U3::NZ] - nx);
		    ShellSPS.push_back(U3Labels);
		}
	}
	GetAllowedSU3xSU2Irreps(n, A, ShellSPS, AllowedIrreps);
}
