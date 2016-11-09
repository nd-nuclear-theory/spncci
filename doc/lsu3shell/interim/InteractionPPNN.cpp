#include <SU3ME/InteractionPPNN.h>
#include <SU3ME/distr2gamma2omega.h>
#include <LSU3/BroadcastDataContainer.h>

#include <sstream>
#include <stack>

using std::cout;
using cpp0x::array;

struct EQUAL_TO
{
    const SU3xSU2_VEC& to_find;
	const std::vector<char>& m_structure;
	public:
	EQUAL_TO(const SU3xSU2_VEC& labels, const std::vector<char>& structure): to_find(labels), m_structure(structure) {};
    inline bool operator()(const CrmeTable* RMETable) const
    {
        return RMETable->IsTensorEqual(m_structure, to_find);
    }
};

template<typename T>
struct Contains {
	private:
	T labels_;
	public:
	Contains(const T labels): labels_(labels) {};
	bool operator()(const std::pair<T, std::vector<TENSOR_STRENGTH> >& element) const
	{
		return element.first == labels_;
	}
};

//////////////////////////////////////////
/////     CTensorGroups methods    ///////
//////////////////////////////////////////

//	Method Couple returns true if a ket configuration Ket_confs which describes
//	a distribution of fermions over HO shells in ket basis state, is projected
//	by a tensor belonging to a CTensorGroup with vector member variables 
//		CGroupStructure::m_dA
//		CGroupStructure::m_nAnnihilators
//		CGroupStructure::m_ShellsT
//	onto a space of states with distribution of fermions identical to a 
//	bra configuration Bra_confs. 
bool CTensorGroup::Couple(const std::vector<unsigned char>& Bra_confs, const std::vector<unsigned char>& Ket_confs, const std::vector<unsigned char>& hoShells)
{
	assert(Bra_confs.size() == Ket_confs.size());

	std::vector<unsigned char>::const_iterator it;
	size_t index;

	for (size_t i = 0; i < m_ShellsT.size(); ++i)	//	iterate over tensor shells
	{
		it = std::find(hoShells.begin(), hoShells.end(), m_ShellsT[i]);	// does bra/ket configurations have any particle at shell m_ShellsT[i] (i-th single shell tensor HO shell) ?
		if (it == hoShells.end())	//	if not, then the CTensorGroup acts on shell which is not present neither in bra nor in the ket state ==> matrix element vanishes
		{
			return false;
		}
		else
		{
			index = it - hoShells.begin();									//	index, where n = m_ShellsT[i] is stored in array hoShells[index] = n
			if ((Ket_confs[index] + m_dA[i]) != Bra_confs[index]) 			//	Ai + dA != Af
			{
				return false;
			}
			else if (Ket_confs[index] < m_nAnnihilators[i]) // Ai + dA == Af, but (nAnnihilators > Ai) ==> {ta x ta ... ta}^{nAnnihilators} |Ai> == vacuum
			{
				return false;
			}
		}
	}
	return true;
}

bool CTensorGroup::CoupleSingleShellTensor(const std::vector<unsigned char>& Ket_confs, const std::vector<unsigned char>& hoShells)
{
	assert(m_dA.size() == 1 && m_ShellsT.size() == 1 && m_dA[0] == 0);

	std::vector<unsigned char>::const_iterator it = std::find(hoShells.begin(), hoShells.end(), m_ShellsT[0]);	// does bra/ket configurations have any particle at shell m_ShellsT[0] ?
	//	if not, then the CTensorGroup acts on shell which is not present neither in bra nor in the ket state ==> matrix element vanishes
	// if (nAnnihilators > Ai) ==> {ta x ta ... ta}^{nAnnihilators} |Ai> == vacuum
	if (it == hoShells.end() || (Ket_confs[it - hoShells.begin()] < m_nAnnihilators[0]))
	{
		return false;
	}
	return true;
}

CTensorGroup::CTensorGroup(	const std::vector<char>& structure, const std::vector<unsigned char>& Shells, const std::vector<char>& dA, const std::vector<unsigned char> nAnnihilators)
							: m_ShellsT(Shells), m_dA(dA), m_nAnnihilators(nAnnihilators)
{
	n_adta_InSST_.resize(m_ShellsT.size());
	for (size_t i = 0; i < m_ShellsT.size(); ++i)
	{
		int nad = std::count(structure.begin(), structure.end(),  m_ShellsT[i]+1);
		int nta = std::count(structure.begin(), structure.end(), -(m_ShellsT[i]+1));
		n_adta_InSST_[i] = nad + nta;
	}
}

void CTensorGroup::ShowTensorStructure() const
{
	std::cout << "n_adta_InSST_:{";
	for (size_t i = 0; i < n_adta_InSST_.size(); ++i)
	{
		std::cout << (int)n_adta_InSST_[i] << " ";
	}
	std::cout << "}\t\t";
	std::cout << "m_ShellsT:{";
	for (size_t i = 0; i < m_ShellsT.size(); ++i)
	{
		std::cout << (int)m_ShellsT[i] << " ";
	}
	std::cout << "}\t\t m_dA:{";
	for (size_t i = 0; i < m_dA.size(); ++i)
	{
		std::cout << (int)m_dA[i] << " ";
	}
	std::cout << "}\t\t m_nAnnihilators:{";
	for (size_t i = 0; i < m_nAnnihilators.size(); ++i)
	{
		std::cout << (int)m_nAnnihilators[i] << " ";
	}
	cout << "}" << endl;
}

void CTensorGroup::ShowTensors() const
{
	std::cout << "number of tensors in group: " << m_Tensors.size() << std::endl;
	for (size_t i = 0; i < m_Tensors.size(); ++i)
	{
		cout << "\t\t";
		m_Tensors[i].first->ShowLabels();
		cout << "\t";
		size_t nCoeffs = m_Tensors[i].first->GetNCoefficients();
		for (size_t j = 0; j < nCoeffs; ++j)
		{
			std::cout << m_Tensors[i].second[j] << " ";
		}
		std::cout << endl;
	}
}

void CTensorGroup::Show(const CBaseSU3Irreps& BaseSU3Irreps) const
{
	ShowTensorStructure();
	//TODO implement more
}

//	This method is called from CInteractionPPNN::CInteractionPPNN()
//	as pCTensorGroup->AddTensor(SingleShellTensorRmeTables, vOmega, dCoeffs);
void CTensorGroup::AddTensor(const std::vector<CrmeTable*>& SingleShellTensorRmeTables, const std::vector<SU3xSU2::LABELS*> vOmega, WigEckSU3SO3CGTable* pSU3SO3CGTable, const std::vector<COEFF_DOUBLE> dCoeffs)
{
	COEFF_DOUBLE* Coeffs = new COEFF_DOUBLE[dCoeffs.size()];
	std::copy(dCoeffs.begin(), dCoeffs.end(), Coeffs);
	m_Tensors.push_back(std::make_pair(new CTensorStructure(SingleShellTensorRmeTables, vOmega, pSU3SO3CGTable), Coeffs));	//	CTensorGroup allocates CTensorStructure dynamically => it is responsible for cleanup
																											//	see ~CTensorGroup()
}

CTensorGroup::~CTensorGroup()
{
	for (size_t i = 0; i < m_Tensors.size(); ++i)
	{
		delete m_Tensors[i].first;		// 	delete CTensorStructure
		delete []m_Tensors[i].second;	//	delete tensor coefficients
	}
}

size_t CTensorGroup::SelectTensorsByGamma(const	UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket, const std::vector<unsigned char>& Ket_confs,
											const std::vector<unsigned char>& BraKetShells, const int phase, std::vector<std::pair<CRMECalculator*, COEFF_DOUBLE*> >& Tensors)
{
	assert(gamma_bra.size() == gamma_ket.size());
	assert(Ket_confs.size() == BraKetShells.size());

	size_t i;
	std::vector<std::pair<CTensorStructure*, COEFF_DOUBLE*> > Selected;
	Selected.reserve(m_Tensors.size()); // the worst case scenario: all tensors in TensorGroup are to be selected 
										//	==> size of Selected will be at most m_Tensor.size()

//	In this loop we select those elements of m_Tensors vector (i.e.
//	pair<CTensorStructure*, DOUBLE_COEFF*>) that can be coupled by: gamma_ket x
//	Tensor --> gamma_bra
	for (i = 0; i < m_Tensors.size(); ++i)
	{
		if (m_Tensors[i].first->Couple(gamma_bra, gamma_ket, BraKetShells, m_ShellsT))
		{
			Selected.push_back(m_Tensors[i]);
		}
	}

//	For each element of the vector Selected (i.e. pair<CTensorStructure*,
//	DOUBLE_COEFF*>) create CRMECalculator, which is the structure capable of
//	calculating rme for a given vector of intershell coupled SU(3)xSU(2) irreps
//	(omega_bra and omega_ket).
	for(i = 0; i < Selected.size(); ++i)
	{
		CRMECalculator *pRC = Selected[i].first->GetRMECalculator(gamma_bra, gamma_ket, Ket_confs, BraKetShells, m_ShellsT, phase); // GetRMECalculator calls new CRMECalculator ==> I am responsible for clean up!!!
		if (pRC != NULL) // pRC == NULL if all rmes of CTensorStructure are vanishing for a given gamma_bra and gamma_ket
		{
			Tensors.push_back(std::make_pair(pRC, Selected[i].second)); //	push back pair<CRMECalculator*, DOUBLE_COEFF*>
		}
	}
	return Tensors.size();	// user of Tensors is responsible for delete Tensors[i].first!!!!
}

size_t CTensorGroup::SelectTensorsByGamma(const	UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket, const std::vector<unsigned char>& Ket_confs,
											const std::vector<unsigned char>& BraKetShells, const int phase, std::vector<std::pair<CRMECalculator*, const CTensorStructure*> >& Tensors)
{
	assert(gamma_bra.size() == gamma_ket.size());
	assert(Ket_confs.size() == BraKetShells.size());

	size_t i;
	std::vector<CTensorStructure*> Selected;
	Selected.reserve(m_Tensors.size()); // the worst case scenario: all tensors in TensorGroup are to be selected 
										//	==> size of Selected will be at most m_Tensor.size()

//	In this loop we select those elements of m_Tensors vector (i.e.
//	pair<CTensorStructure*, DOUBLE_COEFF*>) that can be coupled by: gamma_ket x
//	Tensor --> gamma_bra
	for (i = 0; i < m_Tensors.size(); ++i)
	{
		if (m_Tensors[i].first->Couple(gamma_bra, gamma_ket, BraKetShells, m_ShellsT))
		{
			Selected.push_back(m_Tensors[i].first);
		}
	}

//	For each element of the vector Selected (i.e. pair<CTensorStructure*,
//	DOUBLE_COEFF*>) create CRMECalculator, which is the structure capable of
//	calculating rme for a given vector of intershell coupled SU(3)xSU(2) irreps
//	(omega_bra and omega_ket).
	for(i = 0; i < Selected.size(); ++i)
	{
		CRMECalculator *pRC = Selected[i]->GetRMECalculator(gamma_bra, gamma_ket, Ket_confs, BraKetShells, m_ShellsT, phase); // GetRMECalculator calls new CRMECalculator ==> I am responsible for clean up!!!
		if (pRC != NULL) // pRC == NULL if all rmes of CTensorStructure are vanishing for a given gamma_bra and gamma_ket
		{
			Tensors.push_back(std::make_pair(pRC, Selected[i])); //	push back pair<CRMECalculator*, DOUBLE_COEFF*>
		}
	}
	return Tensors.size();	// user of Tensors is responsible for delete Tensors[i].first!!!!
}

size_t CTensorGroup::SelectTensorsByGamma(const nucleon::Type type, const	UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket, const std::vector<unsigned char>& Ket_confs, const std::vector<unsigned char>& BraKetShells, const int phase, std::vector<std::pair<CRMECalculator*, COEFF_DOUBLE*> >& Tensors)
{
	assert(gamma_bra.size() == gamma_ket.size());
	assert(Ket_confs.size() == BraKetShells.size());

	size_t i;
	std::vector<std::pair<CTensorStructure*, COEFF_DOUBLE*> > Selected;
	Selected.reserve(m_Tensors.size()); // the worst case scenario: all tensors in TensorGroup are to be selected 
										//	==> size of Selected will be at most m_Tensor.size()

//	In this loop we select those elements of m_Tensors vector (i.e.
//	pair<CTensorStructure*, DOUBLE_COEFF*>) that can be coupled by: gamma_ket x
//	Tensor --> gamma_bra
	for (i = 0; i < m_Tensors.size(); ++i)
	{
		if (m_Tensors[i].first->Couple(gamma_bra, gamma_ket, BraKetShells, m_ShellsT))
		{
			Selected.push_back(m_Tensors[i]);
		}
	}

//	For each element of the vector Selected (i.e. pair<CTensorStructure*,
//	DOUBLE_COEFF*>) create CRMECalculator, which is the structure capable of
//	calculating rme for a given vector of intershell coupled SU(3)xSU(2) irreps
//	(omega_bra and omega_ket).
	for(i = 0; i < Selected.size(); ++i)
	{
		CRMECalculator *pRC = Selected[i].first->GetRMECalculator(gamma_bra, gamma_ket, Ket_confs, BraKetShells, m_ShellsT, phase); // GetRMECalculator calls new CRMECalculator ==> I am responsible for clean up!!!
		if (pRC != NULL) // pRC == NULL if all rmes of CTensorStructure are vanishing for a given gamma_bra and gamma_ket
		{
			COEFF_DOUBLE* pcoeffs = Selected[i].second + type*Selected[i].first->GetNCoefficients()/2;
			Tensors.push_back(std::make_pair(pRC, pcoeffs)); //	push back pair<CRMECalculator*, DOUBLE_COEFF*>
		}
	}
	return Tensors.size();	// user of Tensors is responsible for delete Tensors[i].first!!!!
}



size_t CTensorGroup::SelectTensorsByGamma(const	UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket, const std::vector<unsigned char>& Ket_confs,
											const std::vector<unsigned char>& BraKetShells, std::vector<std::pair<SU3xSU2_VEC, SU3xSU2_VEC> >& Tensors)
{
	assert(gamma_bra.size() == gamma_ket.size());
	assert(Ket_confs.size() == BraKetShells.size());

	size_t i;
	std::vector<std::pair<CTensorStructure*, COEFF_DOUBLE*> > Selected;
	Selected.reserve(m_Tensors.size()); // the worst case scenario: all tensors in TensorGroup are to be selected 
										//	==> size of Selected will be at most m_Tensor.size()

//	In this loop we select those elements of m_Tensors vector (i.e.
//	pair<CTensorStructure*, DOUBLE_COEFF*>) that can be coupled by: gamma_ket x
//	Tensor --> gamma_bra
	for (i = 0; i < m_Tensors.size(); ++i)
	{
		if (m_Tensors[i].first->Couple(gamma_bra, gamma_ket, BraKetShells, m_ShellsT))
		{
			SU3xSU2_VEC Gamma_t, Omega_t;
			m_Tensors[i].first->GetGammaOmega(BraKetShells,  m_ShellsT, Gamma_t, Omega_t);
			Tensors.push_back(std::make_pair(Gamma_t, Omega_t));
		}
	}
	return Tensors.size();	// user of Tensors is responsible for delete Tensors[i].first!!!!
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//    					CInteractionPPNN 
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void CInteractionPPNN::BroadcastInteractionData(std::vector<int8_t>& structures_tot, std::vector<uint16_t>& number_tensors_tot, std::vector<uint8_t>& rholmmu2S_tot, std::vector<TENSOR_STRENGTH>& coeffs_tot, bool bcastRmeTable)
{
	BroadcastDataContainer(structures_tot, "structures_tot");
	BroadcastDataContainer(number_tensors_tot, "number_tensors_tot");
	BroadcastDataContainer(rholmmu2S_tot, "rholmmu2S_tot");
	BroadcastDataContainer(coeffs_tot, "coeffs_tot");

	if (bcastRmeTable)	// ==> broadcast rmeLookUpTableContainer_
	{
		BroadcastDataContainer(rmeLookUpTableContainer_, "rmeLookUpTableContainer_");
	}
}


CInteractionPPNN::~CInteractionPPNN()
{
	for (int i = 0; i < 9; ++i)
	{
		for (int j = 0; j < m_TensorGroups[i].size(); ++j)
		{
			delete m_TensorGroups[i][j];
		}
	}
}

void PrepareTensor(const std::vector<unsigned char>& hoShells, const std::vector<char>& n_adta_InSST, const std::vector<unsigned char>& ShellsT, std::vector<unsigned char>& tensor)
{ 
	std::vector<unsigned char>::const_iterator cit;
	tensor.assign(tensor.size(), 0);	//	tensor is an array whose element [i] hold the number of a+/a operators in HO shell n == hoShells[i]
	for(size_t i = 0; i < hoShells.size(); ++i)	// iterate HO shells with non-zero number of fermions 
	{
		cit = std::find(ShellsT.begin(), ShellsT.end(), hoShells[i]);	// *cit == HO shell n == hoShells[i]
		if (cit != ShellsT.end()) // there is a single-shell tensor acting on HO shell n == hoShells[i]
		{
			tensor[i] = n_adta_InSST[cit - ShellsT.begin()];	// store number of a+/a operators ==> needed for phase calculation
		}
	}
}

// Documentation can be found in devel/docs/CodeSnipets/TensorGroupSelection.tex
size_t  CInteractionPPNN::GetRelevantTensorGroups(	const std::vector<unsigned char>& Bra_confs, const std::vector<unsigned char>& Ket_confs, 
													const std::vector<unsigned char>& hoShells, std::vector<CTensorGroup*>& vTensorGroups, 
													std::vector<int>& vphases) const
{
	int delta = 0;
	for (int i = 0; i < Bra_confs.size(); ++i)
	{
		delta += abs(Bra_confs[i] - Ket_confs[i]);
		if (delta > 4)
		{
			assert(vTensorGroups.empty());
			return vTensorGroups.size(); //	==> 2B interaction could not couple
		}
	}

	std::vector<unsigned char> tensor(hoShells.size());

	if (delta == 4)
	{
		int diff;
		std::vector<unsigned char> nAnnihilators;
		std::vector<char> dA;
		std::vector<unsigned char> shellsT;	//	HO shells with active single-shell tensors.

		dA.reserve(Bra_confs.size());
		shellsT.reserve(Bra_confs.size());
		nAnnihilators.reserve(Bra_confs.size());

		for (int i = 0; i < Bra_confs.size(); ++i)
		{
			diff = Bra_confs[i] - Ket_confs[i];
			if (diff)
			{
				dA.push_back(diff);
				shellsT.push_back(hoShells[i]);
				nAnnihilators.push_back(abs(diff)*(diff < 0));
			}
		}

		int index = delta + dA.size();
		for (int j = 0; j < m_TensorGroups[index].size(); ++j)	// iterate over all TensorGroups in the family ...
		{
			CTensorGroup* pTG = m_TensorGroups[index][j];
			if (pTG->m_dA == dA && pTG->m_ShellsT == shellsT && pTG->m_nAnnihilators == nAnnihilators)
			{
				PrepareTensor(hoShells, pTG->n_adta_InSST(), pTG->ShellsT(), tensor);
				vphases.push_back(rme_uncoupling_phase(Bra_confs, tensor, Ket_confs));	//	calculate and store a phase due to the exchanges of fermion creation/annihilation operators [see relation (5.6)]
				vTensorGroups.push_back(pTG);
				return vTensorGroups.size(); // we are done for delta=4 there exist just a single group of tensors that satisfies m_dA==dA && m_ShellsT == shellsT && m_nAnnihilators == nAnnihilators
			}
		}
	}
	else if (delta == 2)
	{
		for (int n = 2; n <= 3; ++n)
		{
			int index = delta + n;
			for (int j = 0; j < m_TensorGroups[index].size(); ++j)	// iterate over all TensorGroups in the family ...
			{
				CTensorGroup* pTG = m_TensorGroups[index][j];
				if (pTG->Couple(Bra_confs, Ket_confs, hoShells)) // Does TensorGroup connet bra and ket configurations ?
				{
					PrepareTensor(hoShells, pTG->n_adta_InSST(), pTG->ShellsT(), tensor);
					vphases.push_back(rme_uncoupling_phase(Bra_confs, tensor, Ket_confs));	//	calculate and store a phase due to the exchanges of fermion creation/annihilation operators [see relation (5.6)]
					vTensorGroups.push_back(pTG);
				}
			}
		}
	}
	else if (delta == 0)
	{
		int index = delta + 1;
		for (int j = 0; j < m_TensorGroups[index].size(); ++j)	// iterate over all TensorGroups in the family ...
		{
			CTensorGroup* pTG = m_TensorGroups[index][j];
			if (pTG->CoupleSingleShellTensor(Ket_confs, hoShells)) // Does TensorGroup connet bra and ket configurations ?
			{
				PrepareTensor(hoShells, pTG->n_adta_InSST(), pTG->ShellsT(), tensor);
				vphases.push_back(rme_uncoupling_phase(Bra_confs, tensor, Ket_confs));	//	calculate and store a phase due to the exchanges of fermion creation/annihilation operators [see relation (5.6)]
				vTensorGroups.push_back(pTG);
			}
		}

		index = delta + 2;
		for (int j = 0; j < m_TensorGroups[index].size(); ++j)	// iterate over all TensorGroups in the family ...
		{
			CTensorGroup* pTG = m_TensorGroups[index][j];
			if (pTG->Couple(Bra_confs, Ket_confs, hoShells)) // Does TensorGroup connet bra and ket configurations ?
			{
				PrepareTensor(hoShells, pTG->n_adta_InSST(), pTG->ShellsT(), tensor);
				vphases.push_back(rme_uncoupling_phase(Bra_confs, tensor, Ket_confs));	//	calculate and store a phase due to the exchanges of fermion creation/annihilation operators [see relation (5.6)]
				vTensorGroups.push_back(pTG);
			}
		}
	}
	return vTensorGroups.size();
}


bool CInteractionPPNN::adtaBelongToModelSpace(const CTuple<int8_t, 2>& n1n2) const 
{
	int ta_shell;
	int ad_shell;

	if (n1n2[0] > 0)
	{
		ad_shell = n1n2[0] - 1;
		ta_shell = abs(n1n2[1]) - 1;
	}
	else
	{
		ad_shell = n1n2[1] - 1;
		ta_shell = abs(n1n2[0]) - 1;
	}

	if (ad_shell > BaseSU3Irreps_.GetLastShell() || ta_shell > BaseSU3Irreps_.GetLastShell())
	{
		return false;
	}

	if (BaseSU3Irreps_.GetAmax(ta_shell) == 0)
	{
		return false;
	}

	if (BaseSU3Irreps_.GetAmax(ad_shell) == 0)
	{
		return false;
	}
	return true;
}

// Files for one-body operators are small => I read all components and decision whether a given tensor term 
// is not vanishing for a given model space is decided in the AddOneBodyOperators method
void CInteractionPPNN::LoadOneBodyOperatorsIntoDataStructures(const std::vector<std::pair<std::string, TENSOR_STRENGTH> >& fileNameStrength, OneBodyInteractionTensors& unique_Tensors, SU2::LABEL& JJ0, SU2::LABEL& MM0)
{
	int n1, n2, lm0, mu0, S2, ll0, jj0, mm0; // it is assumed jj0 and mm0 are constant for all files in fileNameStrength

	CTuple<int8_t, 2> structure; // two body interaction has structure {na, nb, 0}
	array<uint8_t, 4> labels;
	size_t nTensors, nCoeffs;	// nTensors = # of tensors belonging to a specific tensor group ... obtained from the interaction file

	std::vector<TENSOR_STRENGTH> coeffs;

	coeffs.reserve(256); // buffer with more then plenty of space available

	for (int ifile = 0; ifile < fileNameStrength.size(); ++ifile)
	{
		TENSOR_STRENGTH strength = fileNameStrength[ifile].second;
		
		std::cout << "\t" << fileNameStrength[ifile].first << "\t\tcoeff:" << strength << std::endl; 

		if (Negligible(strength))
		{
			continue;
		}
		std::ifstream interactionFile(fileNameStrength[ifile].first.c_str());
		if (!interactionFile)
		{
			std::cerr << "Could not open file " << fileNameStrength[ifile].first << " with proton-proton & neutron-neutron interaction." << std::endl;
			exit(EXIT_FAILURE);
		}

		interactionFile >> jj0;
		interactionFile >> mm0;

		if (ifile > 0 && (JJ0 != jj0 || MM0 != mm0))
		{
			std::cerr << "Error: trying to load one body operators with different JJ0/MM0 values!" << std::endl;
			std::cerr << "One body operator files '" << fileNameStrength[ifile].first << "' contains tensor with JJ0=" << jj0 << " and MM0=" << mm0;
			std::cerr << ", while the previous one(s) contain(s) JJ0=" << (int)JJ0 << " and MM0:" << (int)MM0 << std::endl;
			exit(EXIT_FAILURE);
		}

		JJ0 = jj0;
		MM0 = mm0;

		while (true)
		{
			interactionFile >> n1 >> n2;
			if (interactionFile.eof()) 
			{
				break;
			}
			structure[0] = n1; structure[1] = n2;
			
			bool bInclude = adtaBelongToModelSpace(structure);

			interactionFile >> nTensors;
			for (int itensor = 0; itensor < nTensors; ++itensor)
			{
				interactionFile >> lm0 >> mu0 >> S2 >> ll0 >> jj0;
				
//				nCoeffs = 2*SU3::kmax(SU3::LABELS(1, lm0, mu0), ll0/2);
				nCoeffs = 2*1;

				coeffs.resize(nCoeffs);
				
				for (int j = 0; j < nCoeffs; ++j)
				{
					interactionFile >> coeffs[j];
				}
				//	multiply each coefficient by coeff

				if (!bInclude)
				{
					continue;
				}

				labels[0] = lm0;
				labels[1] = mu0;
				labels[2] = S2;
				labels[3] = ll0;

				std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(), std::bind2nd(std::multiplies<TENSOR_STRENGTH>(), strength));

				OneBodyInteractionTensors::iterator current_structure = unique_Tensors.find(structure);
				if (current_structure == unique_Tensors.end())
				{
					unique_Tensors.insert(std::make_pair(structure, ONEBODYTENSORS_WITH_STRENGTHS(1, std::make_pair(labels, coeffs))));
				}
				else
				{
					ONEBODYTENSORS_WITH_STRENGTHS::iterator it = std::find_if(current_structure->second.begin(), current_structure->second.end(), Contains<array<uint8_t, 4> >(labels));
					if (it == current_structure->second.end())
					{
						current_structure->second.push_back(std::make_pair(labels, coeffs));
					}
					else
					{
						for (int icoeff = 0; icoeff < it->second.size(); ++icoeff)
						{ 
							it->second[icoeff] += coeffs[icoeff];
						}
					}
				}
			}
		}
	}
}


void CInteractionPPNN::AddOneBodyOperators(const int my_rank, const std::vector<std::pair<std::string, TENSOR_STRENGTH> >& fileNameStrength, SU2::LABEL& jj0, SU2::LABEL& mm0)
{
//	ithe element contains {n1, n2} structure
	std::vector<int8_t> structures_tot;			// .size() == 2*number_records
//	ith element contains number of tensors associated with ith pair {n1 n2}
	std::vector<uint16_t> number_tensors_tot;	// .size() == number_records
	std::vector<uint8_t> rholmmu2S_tot;
	std::vector<TENSOR_STRENGTH> coeffs_tot;

	boost::chrono::system_clock::time_point start;
	boost::chrono::duration<double> duration;

	if (my_rank == 0)
	{
		uint16_t number_non_vanishing_tensors; 
		uint16_t nTensors;
		OneBodyInteractionTensors unique_Tensors;
	

		start = boost::chrono::system_clock::now();

		LoadOneBodyOperatorsIntoDataStructures(fileNameStrength, unique_Tensors, jj0, mm0);

		duration = boost::chrono::system_clock::now() - start;
		cout << "CInteractionPPNN::AddOneBodyOperators\t\t Time for loading ... " << duration << endl;
	

		for (OneBodyInteractionTensors::iterator current_structure = unique_Tensors.begin(); current_structure != unique_Tensors.end(); ++current_structure)
		{
			number_non_vanishing_tensors = 0;
			nTensors = current_structure->second.size();
			for (uint32_t itensor = 0; itensor < nTensors; ++itensor)
			{
				const std::vector<TENSOR_STRENGTH>& coeffs = current_structure->second[itensor].second;
				if (std::count_if(coeffs.begin(), coeffs.end(), Negligible) != coeffs.size())
				{
					const array<uint8_t, 4>& labels = current_structure->second[itensor].first;
					coeffs_tot.insert(coeffs_tot.end(), coeffs.begin(), coeffs.end());

					rholmmu2S_tot.push_back(labels[0]); // lm0
					rholmmu2S_tot.push_back(labels[1]);	// mu0
					rholmmu2S_tot.push_back(labels[2]);	// S2
					rholmmu2S_tot.push_back(labels[3]); // ll0
					number_non_vanishing_tensors++;
				}
			}

			if (number_non_vanishing_tensors > 0)
			{
				structures_tot.insert(structures_tot.end(), current_structure->first.begin(), current_structure->first.end());
				number_tensors_tot.push_back(number_non_vanishing_tensors);
			}
		}
	}

	MPI_Bcast((void*)&jj0, 1, boost::mpi::get_mpi_datatype(jj0), 0, MPI_COMM_WORLD);
	MPI_Bcast((void*)&mm0, 1, boost::mpi::get_mpi_datatype(mm0), 0, MPI_COMM_WORLD);

	bool bcastRmeTensors = false; // ==> rmeLookUpTableContainer_ will not be broadcasted since it is not being prepared in LoadOneBodyOperatorsIntoDataStructures
	BroadcastInteractionData(structures_tot, number_tensors_tot, rholmmu2S_tot, coeffs_tot, bcastRmeTensors);

	uint32_t number_records = number_tensors_tot.size();
	uint16_t nTensors, nCoeffs;
	uint32_t label_index(0), coeffs_index(0);
	std::vector<char> structure(3, 0);
	SU3xSU2_VEC TensorLabels(1);	//	one-body interaction ==> only one SU3xSU2 label is needed 
	int lm0, mu0, S2, ll0;

	for (int irecord = 0; irecord < number_records; ++irecord)
	{
		structure[0] = structures_tot[2*irecord]; 
		structure[1] = structures_tot[2*irecord + 1]; 

		std::vector<std::vector<char> > vSingleShell_structures;	//	array of structures describing construction of the single-shell tensors
		std::vector<unsigned char> ShellsT;							//	array of HO shells of all single-shell tensors present in a general tensor
		char nShellsT;												//	number of shells in tensor => for one-body interaction max(nShellsT) = 2
		std::vector<std::vector<SU3xSU2::LABELS*> > vSingleShellTensorLabels; 	
		std::vector<SU3xSU2::LABELS*> vOmega;

	//	transform structure describing a general one-body tensor and its
	//	SU3xSU2 label into two possible cases:
	//	(A) Input: structure = {1 -3 0}
	//	TensorLabels = {IR0}
	//	Output:
	//	ShellsT = {0, 2}
	//	vSingleShell_structures = [ {1}, {-3}]
	//	vSingleShellTensorLabels = [empty, empty]
	//	vOmega = {&IR0}
	//
	//	(B) Input: structure = {3 -3 0}
	//	TensorLabels = {IR0}
	//	Output:
	//	ShellsT = {2}
	//	vSingleShell_structures = [{3 -3 0}]
	//	vSingleShellTensorLabels = [{&IR0}]
	//	vOmega = empty
		StructureToSingleShellTensors(structure, ShellsT, TensorLabels, vSingleShell_structures, vSingleShellTensorLabels, vOmega);

		nShellsT = ShellsT.size();

		assert(nShellsT == vSingleShell_structures.size());		//	each for each shell in ShellsT[i] there exists one single shell operator with structure stored at vSingleShellTensor[i]
		assert(nShellsT == vSingleShellTensorLabels.size()); 	//	each for each shell in ShellsT[i] there exists one single shell operator with structure stored at vSingleShellTensor[i]

		std::vector<char> dA(ShellsT.size(), 0);				//	case (A): dA = {+1, -1} or {-1, +1}; case (B): dA = {0}
		std::vector<unsigned char> nAnnihilators(nShellsT, 0);	//	case (A): nAnnihilators = {0, 1} or {1, 0}; case (B): nAnnihilators = {1}
		CTensorGroup* pCTensorGroup(NULL);  
		int delta = 0;

		for (int i = 0; i < nShellsT; ++i)
		{
			dA[i] = std::count_if(vSingleShell_structures[i].begin(), vSingleShell_structures[i].end(), std::bind2nd(std::greater<char>(), 0)) // #creations - #annihilations
					- 
					std::count_if(vSingleShell_structures[i].begin(), vSingleShell_structures[i].end(), std::bind2nd(std::less<char>(), 0));

			nAnnihilators[i] = std::count_if(vSingleShell_structures[i].begin(), vSingleShell_structures[i].end(), std::bind2nd(std::less<char>(), 0));
			delta += abs(dA[i]);
		}
		assert(nShellsT == 1 || nShellsT == 2);
		//	STEP #1: use structure, ShellsT, dA, nAnnihilators to create class CTensorGroup
		if (nShellsT == 1)	//	case (B)
		{
			assert(delta == 0); // since we are dealing with a+_i a_i tensor
			int index = delta + nShellsT;
			m_TensorGroups[index].push_back(new CTensorGroup(structure, ShellsT, dA, nAnnihilators));
			pCTensorGroup = m_TensorGroups[index].back();	//	pCTensorGroup is a pointer to just added class CTensorGroup
															//	all tensors that are read from file will be stored
															//	in this tensor group
		}
		else if (nShellsT == 2)	//	case (A) ==> dA = {1, -1} or {-1, 1}
		{
			assert(std::count(dA.begin(), dA.end(), 0) == 0);
			assert(delta == 2); // since we are dealing with a+_i a_j tensors ==> delta = 2;
			int index = delta + nShellsT;
			m_TensorGroups[index].push_back(new CTensorGroup(structure, ShellsT, dA, nAnnihilators));
			pCTensorGroup = m_TensorGroups[index].back();
		}

		nTensors = number_tensors_tot[irecord];

		pCTensorGroup->m_Tensors.reserve(nTensors); // reserve enough of memmory for storage of pointers to rme Tables and coefficients

//		iterate over tensors in the CTensorGroup in the input file 
		for (int itensor = 0; itensor < nTensors; ++itensor)
		{
			lm0 = rholmmu2S_tot[label_index++];
			mu0 = rholmmu2S_tot[label_index++];
			S2  = rholmmu2S_tot[label_index++];
			ll0 = rholmmu2S_tot[label_index++];

			TensorLabels[0] = SU3xSU2::LABELS(1, lm0, mu0, S2);
//			nCoeffs = 2*SU3::kmax(TensorLabels[0], ll0/2);
			nCoeffs = 2*1;

			std::vector<TENSOR_STRENGTH> coeffs(&coeffs_tot[coeffs_index], &coeffs_tot[coeffs_index] + nCoeffs);
			coeffs_index += nCoeffs;

			std::vector<CrmeTable*> SingleShellTensorRmeTables(ShellsT.size(), NULL); // number of single-shell tensors = number of active shells
			//	Generate a vector of pointers to single-shell tensors rme tables
			for (int j = 0; j < ShellsT.size(); ++j)
			{
				//	Find rme table in memory or read it from an input file. If input file does not exist => expection is thrown
				SingleShellTensorRmeTables[j] = rmeLookUpTableContainer_.GetRMETablePointer(ShellsT[j], nAnnihilators[j], dA[j], vSingleShell_structures[j], vSingleShellTensorLabels[j], true, log_is_on(), log_file_);
				assert(SingleShellTensorRmeTables[j]); //	since GetRMETablePointer should never return NULL !
			}
			//	get pointer to table of set of SU3SO3 CG coefficients of
			//	type <* *; IR0 * L0|| * *>_{*}				
			WigEckSU3SO3CGTable* pSU3SO3CGTable = wigEckSU3SO3CGTablesLookUpContainer_.GetSU3SO3CGTablePointer(TensorLabels[0], ll0);
			assert(pSU3SO3CGTable); // since GetSU3SO3CGTablePointer shall never return NULL

			//	stores this tensor in form of CTensorStructure that contains:
			//	(1) pointers to single-shell tensor rme tables
			//	(2) vOmega ... SU(3)xSU2 labels of inter-shell coupling
			pCTensorGroup->AddTensor(SingleShellTensorRmeTables, vOmega, pSU3SO3CGTable, coeffs);
			assert(nCoeffs == (pCTensorGroup->m_Tensors.back()).first->GetNCoefficients());
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//	In some cases we do not need to consider a given group of tensors, since
//	the model space is not big enough. 
//	Couple of annihilation operators ta(n1) x ta(n2) can safely act on state that has
//  at least single nucleon in n1 shell and a single nucleon in n2 shell and the rest of nucleons in core.
//  Such a state spans Nhw space, where if shells n1 or n2 are higher than valence shell then 
//  the lowest allowed state Nhw = (n1 - (Is n1 above core?)*valence shell) + (n2 - (Is n2 above core?)*valence shell)
bool CInteractionPPNN::BelongToModelSpace(const std::vector<char>& structure, const std::vector<unsigned char>& ShellsT, const std::vector<char>& dA, const std::vector<unsigned char>& nAnnihilators) const 
{
	assert(structure[6] == 0);

	int adxad_Nhw(0), taxta_Nhw(0); 
	int nShellsT = dA.size();
	for (int i = 0; i < 6; ++i)
	{
		if (structure[i] > 0) // ==> creation operator
		{
			int shell = (structure[i] - 1);
			if (shell > BaseSU3Irreps_.valence_shell()) // if false => nucleon belongs to the core ==> does not add to Nhw
			{										
				// particle at shell n has n - n_valence HO qanta above the valence configuration
				adxad_Nhw += shell - BaseSU3Irreps_.valence_shell();
			}
		}
		else if (structure[i] < 0) // ==> annihilation operator
		{
			int shell = (abs(structure[i]) - 1);
			if (shell > BaseSU3Irreps_.valence_shell())
			{
				// particle at shell n has n - n_valence HO qanta above the valence configuration
				taxta_Nhw += shell - BaseSU3Irreps_.valence_shell();
			}
		}
	}

	bool bInclude = true;
	if (adxad_Nhw > BaseSU3Irreps_.Nmax() || taxta_Nhw > BaseSU3Irreps_.Nmax())
	{
		bInclude = false;
	}
	else
	{
		for (int i = 0; i < nShellsT; ++i)
		{
			int Amax = BaseSU3Irreps_.GetAmax(ShellsT[i]);	//	shell accomodates at most Amax fermions
			int Ai_min = nAnnihilators[i];

//	First condition:
//	tensor needs at least Ai_min fermions in ket, but shell accomodates at most Amax fermions which is less then Ai_min
//	example:	a-a- ... Ai_min = 2, but for a given shell n there is at most one, i.e. Amax = 1, fermion.
//	second condition:
//	action of a single-shell tensor on Ai_min fermions results in state with Af = Ai_min + dA fermions.
//	If Af > Amax ==> |Af ...> lies beyond the current model space => exclude such a tensor
			if (!(Ai_min <= Amax) || !((Ai_min + dA[i]) <= Amax))
			{
				bInclude = false;
				break;
			}
		}
	}
	return bInclude;
}

void CInteractionPPNN::LoadTwoBodyOperatorsIntoDataStructures(const std::vector<std::pair<std::string, TENSOR_STRENGTH> >& fileNameStrength, PPNNInteractionTensors& unique_Tensors)
{
	std::vector<TENSOR_STRENGTH> coeffs;

	std::vector<char> structure(7);
	SU3xSU2_VEC tensorLabels(3);

	std::vector<unsigned char> active_shells; 
	std::vector<std::vector<char> > single_shell_structures;
	std::vector<std::vector<SU3xSU2::LABELS*> > vSingleShellTensorLabels;
	std::vector<SU3xSU2::LABELS*> vOmega;
	std::vector<char> dA;
	std::vector<unsigned char> nAnnihilators;

//	reserve memory so we do not need to reallocate during run time
	active_shells.reserve(64);
	single_shell_structures.reserve(64);
	vSingleShellTensorLabels.reserve(64);
	vOmega.reserve(64);
	dA.reserve(64);
	nAnnihilators.reserve(64);

	for (int ifile = 0; ifile < fileNameStrength.size(); ++ifile)
	{
		TENSOR_STRENGTH strength = fileNameStrength[ifile].second;
		
		std::cout << "\t" << fileNameStrength[ifile].first << "\t\tcoeff:" << strength << std::endl; 

		if (Negligible(strength))
		{
			continue;
		}
		std::ifstream interactionFile(fileNameStrength[ifile].first.c_str());
		if (!interactionFile)
		{
			std::cerr << "Could not open file " << fileNameStrength[ifile].first << " with proton-proton & neutron-neutron interaction." << std::endl;
			exit(EXIT_FAILURE);
		}

		while (true)
		{
			int n1, n2, n3, n4, n5, n6, nTensors, irho, lm, mu, S2, nCoeffs;

			interactionFile >> n1 >> n2 >> n3 >> n4 >> n5 >> n6;
			if (interactionFile.eof())
			{
				break;
			}
			structure[0] = n1; 
			structure[1] = n2; 
			structure[2] = n3; 
			structure[3] = n4; 
			structure[4] = n5; 
			structure[5] = n6;
			structure[6] = 0;

			active_shells.resize(0);
			single_shell_structures.resize(0);
			vSingleShellTensorLabels.resize(0);
			vOmega.resize(0);

			StructureToSingleShellTensors(structure, active_shells, tensorLabels, single_shell_structures, vSingleShellTensorLabels, vOmega);

			dA.resize(active_shells.size());
			nAnnihilators.resize(active_shells.size());

			for (int ishell = 0; ishell < active_shells.size(); ++ishell)
			{
				int num_creation = std::count_if(single_shell_structures[ishell].begin(), single_shell_structures[ishell].end(), std::bind2nd(std::greater<char>(), 0));
				nAnnihilators[ishell] = std::count_if(single_shell_structures[ishell].begin(), single_shell_structures[ishell].end(), std::bind2nd(std::less<char>(), 0));
				dA[ishell] = num_creation - nAnnihilators[ishell];
			}

			bool bInclude = BelongToModelSpace(structure, active_shells, dA, nAnnihilators);

			interactionFile >> nTensors;
//	create list of all proton {a+a}_p and neutron {a+a}_n labels that are
//	stored in a given interaction file with {n1, n2, 0, n4, n5, 0}
			for (size_t itensor = 0; itensor < nTensors; ++itensor)
			{
				interactionFile >> irho >> lm >> mu >> S2;
				tensorLabels[0] = SU3xSU2::LABELS(irho, lm, mu, S2);

				interactionFile >> irho >> lm >> mu >> S2;
				tensorLabels[1] = SU3xSU2::LABELS(irho, lm, mu, S2);
				
				interactionFile >> irho >> lm >> mu >> S2;
				tensorLabels[2] = SU3xSU2::LABELS(irho, lm, mu, S2);

//	Here we assume that tensors are scalars, i.e. L0 == S0
				int L0 = S2/2;
//				nCoeffs = 2*(irho*SU3::kmax(tensorLabels[2], L0));
				nCoeffs = 2*(irho*1);
				coeffs.resize(nCoeffs);
				for (size_t icoeff = 0; icoeff < nCoeffs; ++icoeff)
				{
					interactionFile >> coeffs[icoeff];
				}

				if (!bInclude)
				{
					continue;
				}

				for (int j = 0; j < active_shells.size(); ++j)
				{
				//	Prepare appropriate rme table, i.e. read it from file and if it does not exist calculate and create such file
					rmeLookUpTableContainer_.PrepareRMETable(active_shells[j], nAnnihilators[j], dA[j], single_shell_structures[j], vSingleShellTensorLabels[j], false, cout);
				}

				std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(), std::bind2nd(std::multiplies<TENSOR_STRENGTH>(), strength));

				PPNNInteractionTensors::iterator current_structure = unique_Tensors.find(structure);
				if (current_structure == unique_Tensors.end())
				{
					std::pair<PPNNInteractionTensors::iterator, bool> result = unique_Tensors.insert(std::make_pair(structure, TENSORS_WITH_STRENGTHS()));
					assert(result.second);
					current_structure = result.first;
				}

				TENSORS_WITH_STRENGTHS::iterator it = std::find_if(current_structure->second.begin(), current_structure->second.end(), Contains<SU3xSU2_VEC>(tensorLabels));
				if (it == current_structure->second.end())
				{
					current_structure->second.push_back(std::make_pair(tensorLabels, coeffs));
				}
				else
				{
					for (int icoeff = 0; icoeff < it->second.size(); ++icoeff)
					{ 
						it->second[icoeff] += coeffs[icoeff];
					}
				}
			}
		}
	}
}
void CInteractionPPNN::LoadTwoBodyOperators(const int my_rank, const std::vector<std::pair<std::string, TENSOR_STRENGTH> >& fileNameStrength)
{
	assert(rmeLookUpTableContainer_.size() == 0);

//	ithe element contains {n1, n2, n3, n4, n5, n6, 0} structure
	std::vector<int8_t> structures_tot;			// .size() == 7*number_records
//	ith element contains number of tensors associated with structures_tot[i]
	std::vector<uint16_t> number_tensors_tot;	// .size() == number_records
	std::vector<uint8_t> rholmmu2S_tot;
	std::vector<TENSOR_STRENGTH> coeffs_tot;

	if (my_rank == 0)
	{
		uint16_t nTensors; 
		uint16_t number_non_vanishing_tensors;
		PPNNInteractionTensors unique_Tensors;

		boost::chrono::system_clock::time_point start = boost::chrono::system_clock::now();

		LoadTwoBodyOperatorsIntoDataStructures(fileNameStrength, unique_Tensors);

		boost::chrono::duration<double> duration = boost::chrono::system_clock::now() - start;
		cout << "Time for loading ... " << duration << endl;
		
		for (PPNNInteractionTensors::iterator current_structure = unique_Tensors.begin(); current_structure != unique_Tensors.end(); ++current_structure)
		{
			number_non_vanishing_tensors = 0;
			nTensors = current_structure->second.size();
			for (uint32_t itensor = 0; itensor < nTensors; ++itensor)
			{
				const std::vector<TENSOR_STRENGTH>& coeffs = current_structure->second[itensor].second;
				if (std::count_if(coeffs.begin(), coeffs.end(), Negligible) != coeffs.size())
				{
					coeffs_tot.insert(coeffs_tot.end(), coeffs.begin(), coeffs.end());
					const SU3xSU2_VEC& tensorLabels = current_structure->second[itensor].first;

					for (uint32_t k = 0; k < 3; ++k)
					{
						rholmmu2S_tot.push_back(tensorLabels[k].rho);
						rholmmu2S_tot.push_back(tensorLabels[k].lm);
						rholmmu2S_tot.push_back(tensorLabels[k].mu);
						rholmmu2S_tot.push_back(tensorLabels[k].S2);
					}
					number_non_vanishing_tensors++;
				}
			}
			
			if (number_non_vanishing_tensors > 0)
			{
				structures_tot.insert(structures_tot.end(), current_structure->first.begin(), current_structure->first.end());
				number_tensors_tot.push_back(number_non_vanishing_tensors);
			}
		}
	}

	bool bcastRmeTensor = true; // ==> broadcast rmeLookUpTableContainer_
	BroadcastInteractionData(structures_tot, number_tensors_tot, rholmmu2S_tot, coeffs_tot, bcastRmeTensor);

	SU3xSU2_VEC tensorLabels(3);
	uint32_t number_records = number_tensors_tot.size();
	uint32_t tensor_index(0), coeffs_index(0);
	std::vector<char> structure(7);
	for (uint32_t irecord = 0; irecord < number_records; ++irecord)
	{
		structure.assign(&structures_tot[irecord*7], &structures_tot[irecord*7] + 7);

		std::vector<std::vector<char> > vSingleShell_structures; 	// 	array of structures describing construction of the single-shell tensors
		std::vector<unsigned char> ShellsT;						//	array of HO shells of all single-shell tensors present in a general tensor
																//	capital T in ShellsT stands for Tensors
		
		// all pointers to SU3xSU2::LABELS in vSingleShellTensorLabels and
		// vOmega vectors point to the elements of vector<SU3xSU2::LABELS>
		// tensorLabels; needed for construction/loading of single-shell tensor
		// rme
		std::vector<std::vector<SU3xSU2::LABELS*> > vSingleShellTensorLabels; 	//	i-th element of vSingleShellTensorLabels contains vector of pointers to 
																				//	SU3xSU2 labels that are needed for 	construction of single-shell tensor			
																				//	It thus has ShellsT.size() elements and it is being filled in StructureToSingleShellTensors 
		std::vector<SU3xSU2::LABELS*> vOmega;	//	contains inter-shell labels of a general tensor 
												//	must have ShellsT.size() - 1 element as {{T^Gamma0 x T^Gamma1}^Omega0 x T^Gamma2}^Omega1 x ... x T^Gamma_{n}}^Omega_{n-1}

	//	transform structure describing a general (i.e. multiple-shell) tensor
	//	and IR0, IR1, and IR2 SU3xSU2 labels into a vector of single-shell
	//	structures and associated inner labels, which are needed for
	//	finding/creating a new table of rmes (class CrmTable), as well as the
	//	labels of inter-shell coupled tensors (vOmega)
		StructureToSingleShellTensors(structure, ShellsT, tensorLabels, vSingleShell_structures, vSingleShellTensorLabels, vOmega);
	//	Example I:
	//	Input:
	//	structure = {1 3 0 -4 -4 0 0}
	//	tensorLabels = {IR0, IR1, IR2}
	//	Output:
	//	ShellsT = {0, 2, 3}
	//	vSingleShell_structures = [ {1}, {3}, {-4, -4, 0}]
	//	vSingleShellTensorLabels = [empty, empty, {&IR1}}
	//	vOmega = {&IR0, &IR2}
	//
	//	Example II:
	//	Input:
	//	structure = {3 3 0 -3 -3 0 0}
	//	TensorLabels = {IR0, IR1, IR2}
	//	Output:
	//	ShellsT = {2}
	//	vSingleShell_structures = [{3 3 0 -3 -3 0 0}]
	//	vSingleShellTensorLabels = [{&IR0, &IR1, &IR2}]
	//	vOmega = empty

// 		number of shells in tensor => for 2B interaction max(nShellsT) = 4	=> It can be represented by char data type								
		int nShellsT = ShellsT.size();

		assert(nShellsT == vSingleShell_structures.size() && nShellsT == vSingleShellTensorLabels.size());

		std::vector<char> dA(ShellsT.size(), 0);
		std::vector<unsigned char> nAnnihilators(nShellsT, 0);
		int delta = 0;
		CTensorGroup* pCTensorGroup(NULL);  

		for (int ishell = 0; ishell < nShellsT; ++ishell)
		{
			int nCreators = std::count_if(vSingleShell_structures[ishell].begin(), vSingleShell_structures[ishell].end(), std::bind2nd(std::greater<char>(), 0));
			nAnnihilators[ishell] = std::count_if(vSingleShell_structures[ishell].begin(), vSingleShell_structures[ishell].end(), std::bind2nd(std::less<char>(), 0));
			dA[ishell] =  nCreators - nAnnihilators[ishell];// #creations - #annihilations
			delta += abs(dA[ishell]);
		}
////////////////////////////////////////////////////////////////////////////		
#ifdef	DEBUG_OUTPUT	// provide more detailed output ... useful for debugging StructureToSingleShellTensors
		ShowVector(structure); std::cout << "\t\t"; ShowVector(ShellsT); 
		cout << endl;
		for (int i = 0; i < nShellsT; ++i)
		{
			std::cout << "n = " << (int)ShellsT[i] << ":\tstructure = {";
			ShowVector(vSingleShell_structures[i]);
			std::cout << "}\t" << "dA = " << (int)dA[i] << "\t" << "nAnnihilators = " << (int)nAnnihilators[i];
			//	print single-shell tensors inner coupling labels
			std::cout << "\t InnerShellCoupling = {";
			if (vSingleShellTensorLabels[i].empty())
			{
				std::cout << "empty";
			}
			else
			{
				for (size_t j = 0; j < vSingleShellTensorLabels[i].size(); ++j)
				{
					if (vSingleShellTensorLabels[i][j] == &tensorLabels[0])
					{
						std::cout << "IR0" << " ";
					}
					else if (vSingleShellTensorLabels[i][j] == &tensorLabels[1])
					{
						std::cout << "IR1" << " ";
					}
					else if (vSingleShellTensorLabels[i][j] == &tensorLabels[2])
					{
						std::cout << "IR2" << " ";
					}
				}
			}
			std::cout << "}" << std::endl;
		}
		std::cout << "vOmega = {";
		if (nShellsT == 1)
		{
			cout << "empty";
		}
		else
		{
			for (i = 0; i < nShellsT - 1; ++i)
			{
				if (vOmega[i]== &tensorLabels[0])
				{
					std::cout << "IR0" << " ";
				}
				else if (vOmega[i] == &tensorLabels[1])
				{
					std::cout << "IR1" << " ";
				}
				else if (vOmega[i] == &tensorLabels[2])
				{
					std::cout << "IR2" << " ";
				}
			}
		}
		std::cout << "}" << std::endl;
#endif		
////////////////////////////////////////////////////////////////////////////		

		//	STEP #1: use structure, ShellsT, dA, nAnnihilators to create class CTensorGroup
		int index = delta + nShellsT;
		m_TensorGroups[index].push_back(new CTensorGroup(structure, ShellsT, dA, nAnnihilators));
		pCTensorGroup = m_TensorGroups[index].back();	//	pCTensorGroup is a pointer to just added class CTensorGroup
														//	all tensors that are read from file will be stored in this tensor group
		int nTensors = number_tensors_tot[irecord];
		pCTensorGroup->m_Tensors.reserve(nTensors); // reserve enough of memmory for storage of pointers to rme Tables and coefficients

		if (my_rank == 0 && log_is_on())
		{
			log_file_ << "IR0\t\tIR1\t\tIR2\n";
		}

//		iterate over tensors in the CTensorGroup in the input file 
		for (int itensor = 0; itensor < nTensors; ++itensor)
		{
			for (int k = 0; k < 3; ++k)
			{
				tensorLabels[k].rho = rholmmu2S_tot[tensor_index++];
				tensorLabels[k].lm  = rholmmu2S_tot[tensor_index++];
				tensorLabels[k].mu  = rholmmu2S_tot[tensor_index++];
				tensorLabels[k].S2  = rholmmu2S_tot[tensor_index++];
			}
			SU3::LABELS IR0(tensorLabels[2]);
			int LL0 = tensorLabels[2].S2;
//			int nCoeffs = 2*tensorLabels[2].rho*SU3::kmax(tensorLabels[2], LL0/2);
			int nCoeffs = 2*tensorLabels[2].rho*1;
			std::vector<TENSOR_STRENGTH> coeffs(&coeffs_tot[coeffs_index], &coeffs_tot[coeffs_index] + nCoeffs);
			coeffs_index += nCoeffs;

			if (std::count_if(coeffs.begin(), coeffs.end(), Negligible) == coeffs.size())
			{
				continue;
			}

//////////////////////////////////////////////////////////
			if (my_rank == 0 && log_is_on())
			{
				log_file_ << tensorLabels[0] << "\t" << tensorLabels[1] << "\t" << tensorLabels[2] << std::endl;
			}
#ifdef DEBUG_OUTPUT
			//	{c_{pp}^{rho0 kappa0}, c_{nn}^{rho0 kappa0}, ... , c_{pp}^{rho0max kappa0max}, c_{nn}^{rho0max kappa0max}}
			//	with index = kappa0*rho0max*2 + rho0*2 + PP/NN
			std::cout << "coeffs = ";
			for (int j = 0; j < nCoeffs; ++j)
			{
				std::cout << coeffs[j] << " ";
			}
			std::cout << std::endl;
#endif
			std::vector<CrmeTable*> SingleShellTensorRmeTables(ShellsT.size(), NULL); // number of single-shell tensors = number of active shells
			//	Generate a vector of pointers to single-shell tensors rme tables
			for (int j = 0; j < ShellsT.size(); ++j)
			{
				//	either find already existing rme table, or generate new
				//	rme table, store it in m_RMETables and return its
				//	address.
				SingleShellTensorRmeTables[j] = rmeLookUpTableContainer_.GetRMETablePointer(ShellsT[j], nAnnihilators[j], dA[j], vSingleShell_structures[j], vSingleShellTensorLabels[j]);
				assert(SingleShellTensorRmeTables[j]); //	since GetRMETablePointer should never return NULL !
			}
			//	get pointer to table of set of SU3SO3 CG coeffiients of
			//	type <* *; IR0 * L0|| * *>_{*}				
			WigEckSU3SO3CGTable* pSU3SO3CGTable = wigEckSU3SO3CGTablesLookUpContainer_.GetSU3SO3CGTablePointer(IR0, LL0);
			assert(pSU3SO3CGTable); // since GetSU3SO3CGTablePointer shall never return NULL

//	This is the best place in code to transform elements stored in
//	dCoeffs.
//	From: (currently)
//	index = k0*rho0max*2 + rho*2 + TYPE /0 == PP 1 == NN/
//	Into:
//	index = iType*k0max*rho0max + k0*rho0max + rho0
//	It may/may not turn to be more fast ...	needs to test it
//	Transform(dCoeffs);

			//	stores this tensor in form of CTensorStructure that contains:
			//	(1) pointers to single-shell tensor rme tables
			//	(2) vOmega ... SU(3)xSU2 labels of inter-shell coupling
			pCTensorGroup->AddTensor(SingleShellTensorRmeTables, vOmega, pSU3SO3CGTable, coeffs);
			assert(nCoeffs == (pCTensorGroup->m_Tensors.back()).first->GetNCoefficients());
		}
	}
}

void CInteractionPPNN::TransformTensorStrengthsIntoPP_NN_structure()
{
	for (size_t i = 0; i < 9; ++i)
	{
		for (size_t j = 0; j < m_TensorGroups[i].size(); ++j)
		{
			std::vector<std::pair<CTensorStructure*, CTensorGroup::COEFF_DOUBLE*> >::iterator tensorStructure_Coeffs = m_TensorGroups[i][j]->m_Tensors.begin();
			std::vector<std::pair<CTensorStructure*, CTensorGroup::COEFF_DOUBLE*> >::iterator last_tensorStructure_Coeffs = m_TensorGroups[i][j]->m_Tensors.end();
			for (; tensorStructure_Coeffs < last_tensorStructure_Coeffs; ++tensorStructure_Coeffs)
			{
				size_t ncoeffsPN = tensorStructure_Coeffs->first->GetNCoefficients();
				size_t ncoeffsP, ncoeffsN;
				ncoeffsP = ncoeffsN = ncoeffsPN/2;	//	 k0max * rho0max
				std::vector<CTensorGroup::COEFF_DOUBLE> oldCoeffs;
				oldCoeffs.assign(tensorStructure_Coeffs->second, tensorStructure_Coeffs->second + ncoeffsPN);

				CTensorGroup::COEFF_DOUBLE* pcoeffsP = tensorStructure_Coeffs->second;
				for (size_t index = 0, k = 0; index < ncoeffsP; k += 2, ++index)
				{
					pcoeffsP[index] = oldCoeffs[k];
				}

				CTensorGroup::COEFF_DOUBLE* pcoeffsN = &tensorStructure_Coeffs->second[ncoeffsP];
				for (size_t index = 0, k = 1; index < ncoeffsN; k += 2, ++index)	//	the first neutron tensor strength has index = 1
				{
					pcoeffsN[index] = oldCoeffs[k];
				}
			}
		}
	}
}

void CInteractionPPNN::GetTensors(std::vector<CTensorStructure*>& tensors) const
{
	for (size_t i = 0; i < 9; ++i)
	{
		for (size_t j = 0; j < m_TensorGroups[i].size(); ++j)
		{
			std::vector<std::pair<CTensorStructure*, CTensorGroup::COEFF_DOUBLE*> >::iterator tensorStructure_Coeffs = m_TensorGroups[i][j]->m_Tensors.begin();
			std::vector<std::pair<CTensorStructure*, CTensorGroup::COEFF_DOUBLE*> >::iterator last_tensorStructure_Coeffs = m_TensorGroups[i][j]->m_Tensors.end();
			for (; tensorStructure_Coeffs < last_tensorStructure_Coeffs; ++tensorStructure_Coeffs)
			{
				tensors.push_back(tensorStructure_Coeffs->first);
			}
		}
	}
}

void CInteractionPPNN::Show() const
{
	for (size_t i = 0; i < 9; ++i)
	{
		for (size_t j = 0; j < m_TensorGroups[i].size(); ++j)
		{
			m_TensorGroups[i][j]->Show(BaseSU3Irreps_);
		}
	}
}
void CInteractionPPNN::ShowTensorGroups() const
{
	std::cout << "structure\t\t\t#nShells\tActiveShells" << "\tdA" << "\t\tnAnnihilators" << endl;
	for (int i = 0; i <= 9; ++i)
	{
		for (int j = 0; j < m_TensorGroups[i].size(); ++j)
		{
			m_TensorGroups[i][j]->ShowTensorStructure();
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
}
void CInteractionPPNN::ShowTensors() const
{
	size_t i, j;
	std::cout << "structure\t\t\t#nShells\tActiveShells" << "\tdA" << "\t\tnAnnihilators" << endl;
	for (i = 0; i < 9; ++i)
	{
		for (j = 0; j < m_TensorGroups[i].size(); ++j)
		{
			m_TensorGroups[i][j]->ShowTensorStructure();
			m_TensorGroups[i][j]->ShowTensors();
		}
		std::cout << std::endl;
	}
}
