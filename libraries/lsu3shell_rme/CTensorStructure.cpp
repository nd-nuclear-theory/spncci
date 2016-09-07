#include <SU3ME/CTensorStructure.h>
#include <SU3ME/distr2gamma2omega.h>
#include <stack>
///////////////////////////////////////////////////////////////////
///						CTensorStructure						///	
///////////////////////////////////////////////////////////////////
CTensorStructure::CTensorStructure(const std::vector<CrmeTable*>& SingleShellTensorRmeTables, const std::vector<SU3xSU2::LABELS*> vOmega, WigEckSU3SO3CGTable* pSU3SO3CGTable)
{
	m_pSU3SO3CGTable = pSU3SO3CGTable;
	m_nShells = SingleShellTensorRmeTables.size();
	m_SingleShellTensors = new CrmeTable*[m_nShells];
	std::copy(SingleShellTensorRmeTables.begin(), SingleShellTensorRmeTables.end(), m_SingleShellTensors);
	if (!vOmega.empty())
	{
		assert(vOmega.size() == (m_nShells - 1));

		m_Omega = new SU3xSU2::LABELS[vOmega.size()];
		for (size_t i = 0; i < vOmega.size(); ++i)
		{
			m_Omega[i] = *vOmega[i];
		}
	}
	else
	{
		m_Omega = NULL; // no inter-shell coupling ==> CTensorStructure decribes a general tensor that is composed by a one single-shell tensor.
	}
}

void CTensorStructure::ShowLabels() const
{
	cout << "structure:[";
//	return;
	for (size_t i = 0; i < m_nShells; ++i)
	{
		const std::vector<char>& structure = m_SingleShellTensors[i]->GetStructure();
		if (structure.size() > 0)
		{
			cout << "{";
			for (size_t j = 0; j < structure.size()-1; ++j)
			{
				std::cout << (int)structure[j] << ", ";
			}
			std::cout << (int)structure.back() << "}";
		}
		else
		{
			std::cout << "{}";
		}
	}
	cout << "]\t\t";


	std::cout << "\tSingleShellTensors = [";
	for (size_t i = 0; i < m_nShells - 1; ++i)
	{
		m_SingleShellTensors[i]->ShowTensorLabels();
		std::cout << ", ";
	}
	m_SingleShellTensors[m_nShells-1]->ShowTensorLabels();
	std::cout << "]\t";

	std::cout << "Omega ={";
	if (m_nShells > 1)
	{
		assert(m_Omega != NULL);
		for (size_t i = 0; i < m_nShells - 2; ++i)
		{
			std::cout << m_Omega[i] << ", ";
		}
		std::cout << m_Omega[m_nShells-2] << "}";
	}
	else
	{
		assert(m_Omega == NULL);
		std::cout << "empty}";
	}
/*
	cout << endl;

	std::vector<char> original_structure;
	SU3xSU2_VEC labels;
	GetStructure(original_structure, labels);
	cout << "original_structure:{";
	for (size_t i = 0; i < original_structure.size(); ++i)
	{
		cout << (int)original_structure[i] << " ";
	}
	cout << "}\t\t\t";
	cout << "original_labels:{";
	for (size_t i = 0; i < labels.size(); ++i)
	{
		cout << labels[i] << " ";
	}
	cout << "}\t\t\t";
*/	
}

size_t CTensorStructure::GetNCoefficients() const
{
	assert(m_nShells > 0);
	if (m_nShells > 1)
	{
		assert(m_Omega != NULL);
	}
	if (m_nShells == 1)
	{
		assert(m_Omega == NULL);
	}

	//	if m_nShells == 1 ===> there is just a single shell tensor and hence Omega_t is empty
	SU3xSU2::LABELS IR0(m_nShells > 1 ? m_Omega[m_nShells - 2] : *m_SingleShellTensors[0]);
//	return 2*(IR0.rho*SU3::kmax(IR0, (m_pSU3SO3CGTable->m_LL0 >> 1))); 
	return 2*(IR0.rho*1);
}

CTensorStructure::~CTensorStructure()
{
	delete []m_SingleShellTensors;
	if (m_Omega)	//	m_Omega could be 
	{
		delete []m_Omega;
	}
}

//	gamma_ket[0]	gamma_ket[1] ...	gamma_ket[max]
//		x				x					x
//	Gamma_t[0]		Gamma_t[1]	...		Gamma_t[max]
//		|				|					|
//		|				|					|
//		|				|					|
//		V				V					V
//	gamma_bra[0]	gamma_bra[1]...		gamma_bra[max]	
//
//	NOTE: if BraKetShell[index] in not in TensorShells ==> Gamma_t[index] = IDENTITY_OPERATOR and gamma_ket[index] == gamma_bra[index]
bool CTensorStructure::Couple(	const UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket,
								const std::vector<unsigned char>& BraKetShells,
								const std::vector<unsigned char>& TensorShells)
{
	size_t i;
	std::vector<unsigned char>::const_iterator it;
	for (i = 0; i < BraKetShells.size(); ++i)
	{
		it = std::find(TensorShells.begin(), TensorShells.end(), BraKetShells[i]);
		if (it == TensorShells.end())	//	==> HO shell BraKetShells[i] is not present in the tensor ==> only the identity operator can yield non-zero 
		{
			if ((gamma_bra[i].mult != gamma_ket[i].mult) || (gamma_bra[i].lm != gamma_ket[i].lm) || (gamma_bra[i].mu != gamma_ket[i].mu) || (gamma_bra[i].S2 != gamma_ket[i].S2))
//			if (!(gamma_bra[i] == gamma_ket[i])) // identity operator implies that gamma_bra and gamma_ket are identical. If not => me's vanish
			{
				return false;
			}
		}
		else	//	There is a single-shell tensor with label IR0 = *m_SingleShellTensors[it - TensorShells.begin()]
				//	which must couple with ket to yield bra
		{
			SU3xSU2::LABELS IR0(*m_SingleShellTensors[it - TensorShells.begin()]);
			if (!SU2::mult(gamma_ket[i].S2, IR0.S2, gamma_bra[i].S2))	//	 gamma_ket.S2 x S0 ---> gamma_bra.S2 ?
			{
				return false;
			}
			if (!SU3::mult(gamma_ket[i], IR0, gamma_bra[i]))			//	 (lm mu)gamma_ket x (lm0 mu0) ---> (lm mu)gamma_bra ?
			{
				return false;
			}
		}
	}
	return true;
}
///////////////////////////////////////////////////////////////////////////////
//	There are two versions of CTensorStructure::GetRMECalculator:
//	(i)		Does not neet to delete RMEs in case pRME == NULL
//			since RMEs (pointers to SU3xSU2::RME) are in the second loop
//			after the first one has been completed.
//	(ii)	Builds RMEs in one loop, but in case pRME == NULL
//			it has to delete all pointers in RMEs vector
///////////////////////////////////////////////////////////////////////////////

//	This method takes as an input:
//	BraKetShells  ... shells occupied either in bra or ket state.
//	TensorsShells ... shells at which a tensor represented by CTensorStructure acts,
//	As the output it returns:
//	Gamma_t: 
//	a set of single-shell tensor labels augmented by
//	the identity operators acting at shells occuring at bra/ket space, but not
//	occuing in CTensorStructure
//	Omega_t:
//	coupling of single-shell operators (see Gamma_t) including coupling with identity operators
void CTensorStructure::GetGammaOmega(const std::vector<unsigned char>& BraKetShells,  const std::vector<unsigned char>& TensorShells, SU3xSU2_VEC& Gamma_t, SU3xSU2_VEC& Omega_t)
{
	assert(Omega_t.empty());
	assert(m_nShells == TensorShells.size());

	std::vector<unsigned char> SSTDistribution(BraKetShells.size(), 0);	// implicitly no SST in shell

/*
 *	looks like the intel compiler does not like SU3xSU2::LABELS::IDENTITY_OPERATOR [error while linking ... missing symbol]
 *	instead of trying to fix it I decided to replace SU3xSU2::LABELS::IDENTITY_OPERATOR with SU3xSU2::LABELS(+1, 0, 0, 0)
	Gamma_t.assign(BraKetShells.size(), SU3xSU2::IDENTITY_OPERATOR);	//	fill the whole vector with identity operator SU(3)xSU(2) labels [i.e. 1(0 0)0]
*/	
	Gamma_t.assign(BraKetShells.size(), SU3xSU2::LABELS(+1, 0, 0, 0));	//	fill the whole vector with identity operator SU(3)xSU(2) labels [i.e. 1(0 0)0]
	Omega_t.resize(m_nShells - 1);	//	if m_nShells == 1 ==> Omega_t = {empty}
	
	//	==> if m_nShells == 1 ==> m_Omega == NULL
	if (!Omega_t.empty())
	{
		Omega_t.assign(m_Omega, m_Omega + m_nShells - 1);	// I need to check that this is correct
	}

	for (size_t i = 0; i < BraKetShells.size(); ++i)
	{
		std::vector<unsigned char>::const_iterator it = std::find(TensorShells.begin(), TensorShells.end(), BraKetShells[i]);
		if (it != TensorShells.end())	// ==> BraKetShells[i] is present in TensorShells ==> find and store pointer to relevant rmes
		{
			size_t index  = it - TensorShells.begin();	//	n = TensorShells[index]
	 		SSTDistribution[i] 	= 1;
			Gamma_t[i] 			= SU3xSU2::LABELS(*m_SingleShellTensors[index]);
		}
	}
	if (Gamma_t.size() > 1) // rule out the case when bra_conf = ket_conf = [0 6 0 0 0 0 ]
	{
		//	This is basically the same transformation as
		AugmentOmegaByVacuumShells(SSTDistribution, Gamma_t, Omega_t);
	}
}

//	This method assumes that gamma_ket[i] x m_SingleShellTensors[i] --> gamma_bra[i]
//	and gamma_ket[j] == gamma_bra[j] if BraKetShells[j] is not in TensorShells
//
//	That is, it is assumed that CTensorStructure::Couple(gamma_bra, gamma_ket ... ) == true
//
//	Construct struct CRMECalculator which performs calculation of rme for
//	multi-shell configurations using SU(3)xSU(2) RME reduction rule.
//
//	to constuct CRMECalculator one needs to
//		(1)	place rme of identity operator on shells which are not in HO shells
//			of CTensorStructure tensor (but are in BraKetShells)
//		(2)	vector of inter-shell coupling which properly takes into account 
//			that CRMECalculator has identity operators
//
//	To get (2) one needs to transform CTensorStructure::m_Omega into
//	CRMECalculator::Omega, which properly takes into account presence of
//	identity operators and their inter-shell coupling
CRMECalculator* CTensorStructure::GetRMECalculator(	const UN::SU3xSU2_VEC& gamma_bra,
													const UN::SU3xSU2_VEC& gamma_ket,
													const std::vector<unsigned char>& Ket_confs,
													const std::vector<unsigned char>& BraKetShells,
													const std::vector<unsigned char>& TensorShells,
													const int phase) const
{
	assert(gamma_bra.size() == gamma_ket.size());
	assert(m_nShells == TensorShells.size());

#ifdef DEBUG_OUTPUT
////////////////////////////////////////////////////////////////	
	std::cout << "CTensorStructure::GetRMECalculator";
	std::cout << std::endl;
////////////////////////////////////////////////////////////////	
#endif	
	std::vector<unsigned char> SSTDistribution(BraKetShells.size(), 0);	// implicitly no SST in shell
	std::vector<SU3xSU2::RME*> RMEs(BraKetShells.size(), NULL);
	std::vector<SU3xSU2::RME::DOUBLE*> RMEsPointers(BraKetShells.size(), NULL);
/*
 *	looks like the intel compiler does not like SU3xSU2::LABELS::IDENTITY_OPERATOR [error while linking ... missing symbol]
 *	instead of trying to fix it I decided to replace SU3xSU2::LABELS::IDENTITY_OPERATOR with SU3xSU2::LABELS(+1, 0, 0, 0)
	SU3xSU2_VEC Gamma_t(BraKetShells.size(), SU3xSU2::LABELS::IDENTITY_OPERATOR);	//	fill the whole vector with identity operator SU(3)xSU(2) labels [i.e. 1(0 0)0]
*/	
	SU3xSU2_VEC Gamma_t(BraKetShells.size(), SU3xSU2::LABELS(+1, 0, 0, 0));	//	fill the whole vector with identity operator SU(3)xSU(2) labels [i.e. 1(0 0)0]
	SU3xSU2_VEC Omega_t(m_nShells - 1);	//	if m_nShells == 1 ==> Omega_t = {empty}
	
	//	==> if m_nShells == 1 ==> m_Omega == NULL
	if (!Omega_t.empty())
	{
		Omega_t.assign(m_Omega, m_Omega + m_nShells - 1);	// I need to check that this is correct
	}

	size_t i, index;
	std::vector<unsigned char>::const_iterator it;

//	This loop finds a pointers to rme of the single-shell operators and stores
//	them in vector RMEsPointers. It also set value 1 in vector SSTDistribution,
//	which describes structure of CTensorStructure in terms of single-shell
//	operators. (1 if there is single-shell operator, 0 if there is identity operator)
//	It also set Gamma_t to be SU3xSU2 label of the single shell operator in case the 
//	latter exist, or it will contain label of identity operator due to initialization:
// 	SU3xSU2_VEC Gamma_t(BraKetShells.size(), SU3xSU2::LABELS::IDENTITY_OPERATOR
// 	If rme of a single-shell operator for a given gamma_ket and gamma_bra does not
// 	exist ==> return NULL.
	for (i = 0; i < BraKetShells.size(); ++i)
	{
		it = std::find(TensorShells.begin(), TensorShells.end(), BraKetShells[i]);
		if (it != TensorShells.end())	// ==> BraKetShells[i] is present in TensorShells ==> find and store pointer to relevant rmes
		{
			index  = it - TensorShells.begin();	//	n = TensorShells[index]
			//	use CrmeTable::GetRME_Ai to find pointer to array or rmes of type:  < gamma_bra[i] ||| T^rho0Gamma[index] ||| gamma_ket[i] >_rhot
			SU3xSU2::RME::DOUBLE* pRME = m_SingleShellTensors[index]->GetRME_Ai(gamma_bra[i], gamma_ket[i], Ket_confs[i]);
			if (pRME == NULL) // ==> <Af = Ai + dA[index]; gamma_bra[i] ||| T[index] ||| Ai = Ket_confs[i]; gamma_ket[i]> is not stored in memory, which means that it is equal to zero
			{
#ifdef DEBUG_OUTPUT
////////////////////////////////////////////////////////////////	
				size_t Ai = Ket_confs[i];
				size_t Af = Ai + m_SingleShellTensors[index]->GetdA();
				std::cout << "rme: " << m_SingleShellTensors[index]->GetRMETableFileName() << " with ";
				cout << "<Af = " << Af << "; " << gamma_bra[i] << " ||| T ||| Ai = " << Ai << "; " << gamma_ket[i] << ">" << endl;
				std::cout << std::endl;
////////////////////////////////////////////////////////////////	
#endif	
				return NULL;	// but in that case rmes of the tensor stored in CTensorStructure vanishes for gamma_bra and gamma_ket ==> return NULL pointer instead of CRMECalculator
			}
			else
			{
	 			SSTDistribution[i] 	= 1;
				Gamma_t[i] 			= SU3xSU2::LABELS(*m_SingleShellTensors[index]);
				RMEsPointers[i] 	= pRME;
			}
		}
	}

	for (i = 0; i < SSTDistribution.size(); ++i)
	{
		if (SSTDistribution[i] == 0)	//	==> HO shell BraKetShells[i] is not present in the tensor ==> only the identity operator can yield non-zero 
		{								//	place rme of identity operator on shells which are not in HO shells of CTensorStructure tensor (but are in BraKetShells)
			assert(gamma_bra[i] == gamma_ket[i]);
			assert(RMEsPointers[i] == NULL);

			RMEs[i] = new SU3xSU2::RME(gamma_bra[i]); 	// this creates rme of identity operator
														// CRMECalculator is responsible for deleting this structure in CRMECalculator::~CRMECalculator()
			//	since Gamma_t[i] contains implicitly SU3xSU2::LABELS::IDENTITY
			//	we do not need to change its value
		}
		else
		{

			int tensor_max = Gamma_t[i].rho;
			//	Notice that RMEs[i] contains information on gamma_bra[i].mult
			//	and gamma_ket[i].mult ==> it knows the true size of array of
			//	rmes at which pRME points
			RMEs[i] = new SU3xSU2::RME(gamma_bra[i], tensor_max, Gamma_t[i], gamma_ket[i], RMEsPointers[i]);	// CRMECalculator is responsible for deleting this structure in CRMECalculator::~CRMECalculator()
		}
	}

#ifdef DEBUG_OUTPUT
////////////////////////////////////////////////////////////////	
	cout << "\t\t\tTransforming: ";
	ShowVector(SSTDistribution);
	cout << "\tGamma_t = ";
	for (i = 0; i < Gamma_t.size(); ++i)
	{
		cout << Gamma_t[i] << "  ";
	}
	cout << "\tOmega_t = ";
	if (Omega_t.empty())
	{
		cout << "empty";
	}
	else
	{
		for (i = 0; i < Omega_t.size(); ++i)
		{
			cout << Omega_t[i] << "  ";
		}
	}
	cout << "\t-->\t";
	cout.flush();
////////////////////////////////////////////////////////////////	
#endif	
	if (Gamma_t.size() > 1) // rule out the case when bra_conf = ket_conf = [0 6 0 0 0 0 ]
	{
		AugmentOmegaByVacuumShells(SSTDistribution, Gamma_t, Omega_t);
	}
#ifdef DEBUG_OUTPUT
////////////////////////////////////////////////////////////////	
	cout << "\tOmega_t = ";
	if (Omega_t.empty())
	{
		cout << "empty";
	}
	else
	{
		for (i = 0; i < Omega_t.size(); ++i)
		{
			cout << Omega_t[i] << "  ";
		}
	}
	cout << endl;
////////////////////////////////////////////////////////////////	
#endif	
	return new CRMECalculator(RMEs, Omega_t, m_pSU3SO3CGTable, phase); // m_pSU3SO3CGTable .... <* *; IR0 * L0 || * *>_*
}

//	This function is used by CInteractionPPNN::LoadTwoBodyOperators in order to
//	be able to get rid of tensors acting on states beyond model space.
void StructureToSingleShellTensors(	
//	input
									const CTuple<char, 7>& structure, 
//	output									
									std::vector<unsigned char>& Shells, //	list of HO shells with a tensor
									std::vector<std::vector<char> >& vSingleShell_structures)

{
	static const char MULTI_SHELL_TENSOR = -128;
	std::stack<char> st;
	std::vector<char> SingleShellTensor;	//structure of single-shell tensor
	char n, left, right;
	
	size_t i = 0;
	size_t indexIR = 0;

	st.push(structure[0]);
	SingleShellTensor.push_back(structure[0]);
	Shells.push_back(abs(structure[0]) - 1);

	for (i = 1; i < 7; ++i)
	{
		n = abs(structure[i])-1;
		if (n != Shells.back() && structure[i] != 0) // second condition make sure we do not put n=-1 into Shells when structure[i] = 0.
		{
			Shells.push_back(n);
		}

		if (structure[i] == 0)
		{
			assert(st.size() > 1);

			left = abs(st.top()); st.pop();
			right = abs(st.top()); st.pop();
			if (left != right)	// this is coupling of two tensors acting on different shells
			{
				vSingleShell_structures.push_back(SingleShellTensor);

				SingleShellTensor.clear();

				st.push(MULTI_SHELL_TENSOR);
			}
			else	// coupling inside of single-shell tensor
			{
				SingleShellTensor.push_back(0);	
				st.push(abs(left));
			}
		}
		else
		{
			if (abs(structure[i]) != abs(st.top()) && st.top() != MULTI_SHELL_TENSOR)
			{
				vSingleShell_structures.push_back(SingleShellTensor);
				SingleShellTensor.clear();
			}
			SingleShellTensor.push_back(structure[i]);
			st.push(structure[i]);
		}
	}
	if (!SingleShellTensor.empty()) // e.g. if input structure = +n +n 0 -n -n 0 0
	{
		vSingleShell_structures.push_back(SingleShellTensor);
	}
}


void StructureToSingleShellTensors(	
//	input
									const std::vector<char> structure, 
//	output									
									std::vector<unsigned char>& Shells, //	list of HO shells with a tensor
									SU3xSU2_VEC& TensorLabels, 
									std::vector<std::vector<char> >& vSingleShell_structures,
									std::vector<std::vector<SU3xSU2::LABELS*> >& vSingleShellTensorLabels,
									std::vector<SU3xSU2::LABELS*>& vOmega)

{
	static const char MULTI_SHELL_TENSOR = -128;
	std::stack<char> st;
	std::vector<char> SingleShellTensor;	//structure of single-shell tensor
	std::vector<SU3xSU2::LABELS*> SingleShellTensorLabel;	// labels of single-shell tensor
	char n, left, right;
	
	size_t i = 0;
	size_t indexIR = 0;

	st.push(structure[0]);
	SingleShellTensor.push_back(structure[0]);
	Shells.push_back(abs(structure[0]) - 1);

	for (i = 1; i < structure.size(); ++i)
	{
		n = abs(structure[i])-1;
		if (n != Shells.back() && structure[i] != 0) // second condition make sure we do not put n=-1 into Shells when structure[i] = 0.
		{
			Shells.push_back(n);
		}

		if (structure[i] == 0)
		{
			assert(st.size() > 1);

			left = abs(st.top()); st.pop();
			right = abs(st.top()); st.pop();
			if (left != right)	// this is coupling of two tensors acting on different shells
			{
				vOmega.push_back(&TensorLabels[indexIR++]);
				vSingleShell_structures.push_back(SingleShellTensor);
				vSingleShellTensorLabels.push_back(SingleShellTensorLabel);

				SingleShellTensor.clear();
				SingleShellTensorLabel.clear();

				st.push(MULTI_SHELL_TENSOR);
			}
			else	// coupling inside of single-shell tensor
			{
				SingleShellTensorLabel.push_back(&TensorLabels[indexIR++]);
				SingleShellTensor.push_back(0);	
				st.push(abs(left));
			}
		}
		else
		{
			if (abs(structure[i]) != abs(st.top()) && st.top() != MULTI_SHELL_TENSOR)
			{
				vSingleShellTensorLabels.push_back(SingleShellTensorLabel);
				vSingleShell_structures.push_back(SingleShellTensor);
				SingleShellTensor.clear();
				SingleShellTensorLabel.clear();
			}
			SingleShellTensor.push_back(structure[i]);
			st.push(structure[i]);
		}
	}
	if (!SingleShellTensor.empty()) // e.g. if input structure = +n +n 0 -n -n 0 0
	{
		vSingleShell_structures.push_back(SingleShellTensor);
		vSingleShellTensorLabels.push_back(SingleShellTensorLabel);
	}
	assert(vOmega.size() == Shells.size()-1);
	assert(vSingleShell_structures.size() == vSingleShellTensorLabels.size());
}




