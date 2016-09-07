#include <SU3ME/CInteractionPN.h>
#include <LSU3/BroadcastDataContainer.h>

#include <numeric>

using cpp0x::tuple;
using cpp0x::make_tuple;
using cpp0x::get;

CTensorGroup_ada::~CTensorGroup_ada()
{
	for (size_t i = 0; i < m_Tensors.size(); ++i)
	{
		delete m_Tensors[i];		// 	delete CTensorStructure
	}
}

bool CTensorGroup_ada::Couple(const std::vector<unsigned char>& Bra_confs, const std::vector<unsigned char>& Ket_confs, const std::vector<unsigned char>& hoShells)
{
	assert(Bra_confs.size() == Ket_confs.size());

	size_t i, index;
	std::vector<unsigned char>::const_iterator it;
	std::vector<unsigned char> Ket_confs_transformed(Ket_confs);	//	configuration resulting from the action of tensor with structure described
																	//	by member vector vaciables of CTensorGroup
	for (i = 0; i < m_ShellsT.size(); ++i)	//	iterate over tensor shells
	{
		it = std::find(hoShells.begin(), hoShells.end(), m_ShellsT[i]);	// does bra/ket configurations have any particle at shell m_ShellsT[i] (i-th single shell tensor HO shell) ?
		if (it == hoShells.end())	//	if not, then the CTensorGroup acts on shell which is not present neither in bra nor in the ket state ==> matrix element vanishes
		{
			return false;
		}
		else
		{
			index = it - hoShells.begin();									//	index, where n = m_ShellsT[i] is stored in array hoShells[index] = n
			Ket_confs_transformed[index] = (m_dA[i] + Ket_confs[index]); 	//	Action of the single shell tensor on Ai fermions in the shell m_Shell[i]
			if (Ket_confs_transformed[index] != Bra_confs[index]) 			//	Ai + dA != Af
			{
				return false;
			}
			else if (((char)Ket_confs[index] - (char)m_nAnnihilators[i]) < 0) // Ai + dA == Af, but (nAnnihilators > Ai) ==> {ta x ta ... ta}^{nAnnihilators} |Ai> == vacuum
			{
				return false;
			}
		}
	}
	return (Ket_confs_transformed == Bra_confs); // T x |ket_confs> = |ket_confs>
}

CTensorGroup_ada::CTensorGroup_ada(	const std::vector<char>& structure, const std::vector<unsigned char>& Shells, const std::vector<char>& dA, const std::vector<unsigned char> nAnnihilators, const unsigned int firstLabelId)
							:m_nAnnihilators(nAnnihilators), m_dA(dA), m_ShellsT(Shells), firstLabelId_(firstLabelId)
{
	assert(structure.size() == 3 && structure[2] == 0);

	structure_[0] = structure[0];
	structure_[1] = structure[1];

	n_adta_InSST_.resize(m_ShellsT.size());
	for (size_t i = 0; i < m_ShellsT.size(); ++i)
	{
		int nad = std::count(structure.begin(), structure.end(),  m_ShellsT[i]+1);
		int nta = std::count(structure.begin(), structure.end(), -(m_ShellsT[i]+1));
		n_adta_InSST_[i] = nad + nta;
	}
}

void CTensorGroup_ada::ShowTensorStructure()
{
	size_t i;

	std::cout << "CTensorGroup_ada" << std::endl; 
	std::cout << "structure: " << (int)structure_[0] << " " << (int)structure_[1] << std::endl;

	std::cout << "#active shells: " << (int)m_ShellsT.size() << "\t\t";
	std::cout << "HO shells: ";
	for (i = 0; i < m_ShellsT.size(); ++i)
	{
		std::cout << (int)m_ShellsT[i] << " ";
	}
	
	std::cout << "\t\t";
	std::cout << "#number of a+ta tensors in each HO shell: ";
	for (i = 0; i < n_adta_InSST_.size(); ++i)
	{
		std::cout << (int)n_adta_InSST_[i] << " ";
	}
	std::cout << "\t\t";
	std::cout << "dA: ";
	for (i = 0; i < m_dA.size(); ++i)
	{
		std::cout << (int)m_dA[i] << " ";
	}
	std::cout << "\t\t";
	std::cout << "#annihilators: ";
	for (i = 0; i < m_nAnnihilators.size(); ++i)
	{
		std::cout << (int)m_nAnnihilators[i] << " ";
	}
	std::cout << std::endl;
}

void CTensorGroup_ada::ShowTensors()
{
	std::cout << m_Tensors.size() << std::endl;
	for (size_t i = 0; i < m_Tensors.size(); ++i)
	{
		cout << "\t\t";
		m_Tensors[i]->ShowLabels();
		std::cout << "\tID: " << firstLabelId_ + i << endl;
	}
}

void CTensorGroup_ada::Show(const CBaseSU3Irreps& BaseSU3Irreps)
{
	ShowTensorStructure();
	//TODO implement more
}

//	This method is called from CInteractionPPNN::CInteractionPPNN()
//	as pCTensorGroup->AddTensor(SingleShellTensorRmeTables, vOmega, dCoeffs);
void CTensorGroup_ada::AddTensor(const std::vector<CrmeTable*>& SingleShellTensorRmeTables, const std::vector<SU3xSU2::LABELS*> vOmega)
{	
//	CTensorGroup_ada allocates CTensorStructure dynamically => it is responsible for cleanup
	m_Tensors.push_back(new CTensorStructure(SingleShellTensorRmeTables, vOmega, NULL));
}

void CTensorGroup_ada::SelectTensorsByGamma(	const UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket, const std::vector<unsigned char>& Ket_confs,
												const std::vector<unsigned char>& hoShells, 
												const int phase, std::vector<std::pair<CRMECalculator*, unsigned int> >& selected_tensors_pn)
{
	assert(gamma_bra.size() == gamma_ket.size());
	assert(Ket_confs.size() == hoShells.size());

	size_t i;
	std::vector<std::pair<CTensorStructure*, unsigned int> > Selected;
	Selected.reserve(m_Tensors.size()); // the worst case scenario: all tensors in TensorGroup are to be selected 
										//	==> size of Selected will be at most m_Tensor.size()

//	In this loop we select those elements of m_Tensors vector (i.e.
//	pair<CTensorStructure*, DOUBLE_COEFF*>) that can be coupled by: gamma_ket x
//	Tensor --> gamma_bra
	for (int itensor = 0; itensor < m_Tensors.size(); ++itensor)
	{
		if (m_Tensors[itensor]->Couple(gamma_bra, gamma_ket, hoShells, m_ShellsT))
		{
			Selected.push_back(std::make_pair(m_Tensors[itensor], firstLabelId_ + itensor));
		}
	}

//	For each element of the vector Selected (i.e. pair<CTensorStructure*,
//	DOUBLE_COEFF*>) create CRMECalculator, which is the structure capable of
//	calculating rme for a given vector of intershell coupled SU(3)xSU(2) irreps
//	(omega_bra and omega_ket).
	for(i = 0; i < Selected.size(); ++i)
	{
		CRMECalculator *pRC = Selected[i].first->GetRMECalculator(gamma_bra, gamma_ket, Ket_confs, hoShells, m_ShellsT, phase); // GetRMECalculator calls new CRMECalculator ==> I am responsible for clean up!!!
		if (pRC != NULL) // pRC == NULL if all rmes of CTensorStructure are vanishing for a given gamma_bra and gamma_ket
		{
			selected_tensors_pn.push_back(std::make_pair(pRC, Selected[i].second)); //	push back pair<CRMECalculator*, DOUBLE_COEFF*>
		}
	}
}

void CTensorGroup_ada::SelectTensorsByGamma(	const UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket, const std::vector<unsigned char>& Ket_confs,
												const std::vector<unsigned char>& hoShells, 
												const int phase, std::vector<std::pair<CRMECalculator*, const CTensorStructure*>>& selected_tensors_pn)
{
	assert(gamma_bra.size() == gamma_ket.size());
	assert(Ket_confs.size() == hoShells.size());

	size_t i;
	std::vector<const CTensorStructure*> Selected;
	Selected.reserve(m_Tensors.size()); // the worst case scenario: all tensors in TensorGroup are to be selected 
										//	==> size of Selected will be at most m_Tensor.size()

//	In this loop we select those elements of m_Tensors vector (i.e.
//	pair<CTensorStructure*, DOUBLE_COEFF*>) that can be coupled by: gamma_ket x
//	Tensor --> gamma_bra
	for (size_t itensor = 0; itensor < m_Tensors.size(); ++itensor)
	{
		if (m_Tensors[itensor]->Couple(gamma_bra, gamma_ket, hoShells, m_ShellsT))
		{
			Selected.push_back(m_Tensors[itensor]);
		}
	}

//	For each element of the vector Selected (i.e. pair<CTensorStructure*,
//	DOUBLE_COEFF*>) create CRMECalculator, which is the structure capable of
//	calculating rme for a given vector of intershell coupled SU(3)xSU(2) irreps
//	(omega_bra and omega_ket).
	for(i = 0; i < Selected.size(); ++i)
	{
		CRMECalculator *pRC = Selected[i]->GetRMECalculator(gamma_bra, gamma_ket, Ket_confs, hoShells, m_ShellsT, phase); // GetRMECalculator calls new CRMECalculator ==> I am responsible for clean up!!!
		if (pRC != NULL) // pRC == NULL if all rmes of CTensorStructure are vanishing for a given gamma_bra and gamma_ket
		{
			selected_tensors_pn.push_back(std::make_pair(pRC, Selected[i])); //	push back pair<CRMECalculator*, DOUBLE_COEFF*>
		}
	}
}

size_t CTensorGroup_ada::SelectTensorsByGamma(const	UN::SU3xSU2_VEC& gamma_bra, const UN::SU3xSU2_VEC& gamma_ket, const std::vector<unsigned char>& Ket_confs,
											const std::vector<unsigned char>& BraKetShells, const int phase, std::vector<CRMECalculator*>& Tensors)
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
		if (m_Tensors[i]->Couple(gamma_bra, gamma_ket, BraKetShells, m_ShellsT))
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
		CRMECalculator *pRC = Selected[i]->GetRMECalculator(gamma_bra, gamma_ket, Ket_confs, BraKetShells, m_ShellsT, phase); // GetRMECalculator calls new CRMECalculator ==> I am responsible for clean up!!!
		if (pRC != NULL) // pRC == NULL if all rmes of CTensorStructure are vanishing for a given gamma_bra and gamma_ket
		{
			Tensors.push_back(pRC); //	push back pair<CRMECalculator*, DOUBLE_COEFF*>
		}
	}
	return Tensors.size();	// user of Tensors is responsible for delete Tensors[i].first!!!!
}


struct Contains {
	private:
	tuple<SU3xSU2::LABELS, SU3xSU2::LABELS, SU3xSU2::LABELS> tensorLabels_;
	public:
	Contains(const tuple<SU3xSU2::LABELS, SU3xSU2::LABELS, SU3xSU2::LABELS>& tensorLabels): tensorLabels_(tensorLabels) {};
	bool operator()(const std::pair<tuple<SU3xSU2::LABELS, SU3xSU2::LABELS, SU3xSU2::LABELS>, std::vector<TENSOR_STRENGTH> >& element) const
	{
		return element.first == tensorLabels_;
	}
};

std::pair<uint32_t, uint32_t> CInteractionPN::LoadPNOperatorsIntoDataStructures(const std::vector<std::pair<std::string, TENSOR_STRENGTH> >& fileNameStrength, PNInteractionTensors& unique_Tensors)
{
	enum SizeOfArrays {kCoeffs = 0, kLabels = 1};

	std::pair<uint32_t, uint32_t> arraySizes(0, 0);

	CTuple<int8_t, 4> structure;	// 4 indices: n1 n2 n4 n5, since structure of PN interaction is n1 n2 0 n3 n4 0
	tuple<SU3xSU2::LABELS, SU3xSU2::LABELS, SU3xSU2::LABELS> tensorLabels;
	std::vector<TENSOR_STRENGTH> coeffs;
	bool bInclude;
	unsigned char max_hoshell = BaseSU3Irreps_.GetLastShell();

	coeffs.reserve(256); // allocate a buffer large enough for any realistic case

	for (int ifile = 0; ifile < fileNameStrength.size(); ++ifile)
	{
		TENSOR_STRENGTH strength = fileNameStrength[ifile].second;

		if (Negligible(strength)) {
			continue;
		}
		
		std::cout << "\t" << fileNameStrength[ifile].first << "\t\tcoeff:" << strength << std::endl; 

		std::ifstream interactionFile(fileNameStrength[ifile].first.c_str());
		if (!interactionFile)
		{
			std::cerr << "Could not open file " << fileNameStrength[ifile].first << " with proton-proton & neutron-neutron interaction." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 0);
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
			structure[2] = n4;
			structure[3] = n5;

			interactionFile >> nTensors;
//	create list of all proton {a+a}_p and neutron {a+a}_n labels that are
//	stored in a given interaction file with {n1, n2, 0, n4, n5, 0}
			for (size_t itensor = 0; itensor < nTensors; ++itensor)
			{
				interactionFile >> irho >> lm >> mu >> S2;
				get<0>(tensorLabels) = SU3xSU2::LABELS(irho, lm, mu, S2);

				interactionFile >> irho >> lm >> mu >> S2;
				get<1>(tensorLabels) = SU3xSU2::LABELS(irho, lm, mu, S2);
				
				interactionFile >> irho >> lm >> mu >> S2;
				get<2>(tensorLabels) = SU3xSU2::LABELS(irho, lm, mu, S2);

//	Here we assume that tensors are scalars, i.e. L0 == S0
//				nCoeffs = irho*SU3::kmax(get<2>(tensorLabels), S2/2);
				nCoeffs = irho*1;
				coeffs.resize(nCoeffs);
				for (size_t icoeff = 0; icoeff < nCoeffs; ++icoeff)
				{
					interactionFile >> coeffs[icoeff];
				}

				bInclude = !((abs(n1) - 1) > max_hoshell || (abs(n2) - 1) > max_hoshell || (abs(n4) - 1) > max_hoshell || (abs(n5) - 1) > max_hoshell);

				if (!bInclude)
				{
					continue;
				}

				std::transform(coeffs.begin(), coeffs.end(), coeffs.begin(), std::bind2nd(std::multiplies<TENSOR_STRENGTH>(), strength));

				PNInteractionTensors::iterator current_structure = unique_Tensors.find(structure);
				if (current_structure == unique_Tensors.end())
				{
					std::pair<PNInteractionTensors::iterator, bool> result = unique_Tensors.insert(std::make_pair(structure, TENSORS_WITH_STRENGTHS()));
					assert(result.second);
					current_structure = result.first;
				}

				TENSORS_WITH_STRENGTHS::iterator it = std::find_if(current_structure->second.begin(), current_structure->second.end(), Contains(tensorLabels));
				if (it == current_structure->second.end())
				{
					current_structure->second.push_back(std::make_pair(tensorLabels, coeffs));

					std::get<kCoeffs>(arraySizes) += coeffs.size();
					std::get<kLabels>(arraySizes) += 10; // lmp mup S2p lmn mun S2n rho0 lm0 mu0 S20: 10 elements
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
	return arraySizes;
}


void CInteractionPN::AddOperators(int my_rank, const std::vector<std::pair<std::string, TENSOR_STRENGTH> >& fileNameStrength)
{
	enum SizeOfArrays {kCoeffs = 0, kLabels = 1};
//	array of {n1, n2, n4, n5} structures
	std::vector<int8_t> structures_tot;			// .size() == 4 * #{n1 n2 n4 n5}
//	ith element contains number of tensors associated with {n1 n2 n4 n4} structure
	std::vector<uint16_t> number_tensors_tot;	// .size() == number_records
	std::vector<uint8_t> quantum_labels_tot;
	std::vector<TENSOR_STRENGTH> coeffs_tot;

	boost::chrono::system_clock::time_point start;
	boost::chrono::duration<double> duration;

	if (my_rank == 0)
	{ 
		uint16_t nTensors, number_non_vanishing_tensors;
		PNInteractionTensors unique_Tensors;

		start = boost::chrono::system_clock::now();
		std::pair<uint32_t, uint32_t> sizeOfCoeffsAndLabels = LoadPNOperatorsIntoDataStructures(fileNameStrength, unique_Tensors);
		duration = boost::chrono::system_clock::now() - start;
		cout << "CInteractionPN::AddOperators\t\tTime for loading PN interaction ... " << duration << endl;
		
		coeffs_tot.reserve(get<kCoeffs>(sizeOfCoeffsAndLabels));
		quantum_labels_tot.reserve(get<kLabels>(sizeOfCoeffsAndLabels));
		structures_tot.reserve(4*unique_Tensors.size()); // unique_Tensors.size(): # of different {n1 n2 n4 n5} structures
		number_tensors_tot.reserve(unique_Tensors.size());

		for (PNInteractionTensors::iterator current_structure = unique_Tensors.begin(); current_structure != unique_Tensors.end(); ++current_structure)
		{
			number_non_vanishing_tensors = 0;
			nTensors = current_structure->second.size();
			for (uint32_t itensor = 0; itensor < nTensors; ++itensor)
			{
				const std::vector<TENSOR_STRENGTH>& coeffs = current_structure->second[itensor].second;
				if (std::count_if(coeffs.begin(), coeffs.end(), Negligible) != coeffs.size())
				{
					coeffs_tot.insert(coeffs_tot.end(), coeffs.begin(), coeffs.end());
					const tuple<SU3xSU2::LABELS, SU3xSU2::LABELS, SU3xSU2::LABELS>& tensorLabels = current_structure->second[itensor].first;

					assert(						 get<0>(tensorLabels).rho == 1); // (n 0) x (0 n') --> rhomax=1 
					quantum_labels_tot.push_back(get<0>(tensorLabels).lm);
					quantum_labels_tot.push_back(get<0>(tensorLabels).mu);
					quantum_labels_tot.push_back(get<0>(tensorLabels).S2);

					assert(						 get<1>(tensorLabels).rho == 1); // (n 0) x (0 n') --> rhomax=1 
					quantum_labels_tot.push_back(get<1>(tensorLabels).lm);
					quantum_labels_tot.push_back(get<1>(tensorLabels).mu);
					quantum_labels_tot.push_back(get<1>(tensorLabels).S2);

					quantum_labels_tot.push_back(get<2>(tensorLabels).rho);
					quantum_labels_tot.push_back(get<2>(tensorLabels).lm);
					quantum_labels_tot.push_back(get<2>(tensorLabels).mu);
					quantum_labels_tot.push_back(get<2>(tensorLabels).S2);

					number_non_vanishing_tensors++;
				}
			}

			if (number_non_vanishing_tensors > 0)
			{
				structures_tot.insert(structures_tot.end(), current_structure->first.begin(), current_structure->first.end());
				number_tensors_tot.push_back(nTensors);
			}
		}
	}

	BroadcastDataContainer(structures_tot, "structures_tot");
	BroadcastDataContainer(number_tensors_tot, "number_tensors_tot");
	BroadcastDataContainer(quantum_labels_tot, "quantum_labels_tot");
	BroadcastDataContainer(coeffs_tot, "coeffs_tot");

	start = boost::chrono::system_clock::now();
//	
//	Step 1: Read input files and store all unique combinations {a+a} and their
//	labels rho(lm mu)S0 in map "structure_with_irreps" structure_with_irreps[na nb] 
//	= {rho (lm mu)S, .... , rho'(lm' mu')S'}, where each element is result
//	of coupling creation and annihilation operators in shells na and nb.
//
	uint32_t number_records = number_tensors_tot.size();
	uint16_t nTensors;
	CTuple<char, 2> structure_proton, structure_neutron;
	SU3xSU2::LABELS adta_p, adta_n;
	SU3xSU2_VEC tensors_proton, tensors_neutron;
	std::map<CTuple<char, 2>, SU3xSU2_VEC> structure_with_irreps;

	tensors_proton.reserve(1024);
	tensors_neutron.reserve(1024);

	for (uint32_t irecord = 0, quantum_index = 0; irecord < number_records; ++irecord)
	{
		int8_t* structure = &structures_tot[4*irecord];

		structure_proton[0]  = structure[0];
		structure_proton[1]  = structure[1];
		structure_neutron[0] = structure[2]; 
		structure_neutron[1] = structure[3];

		tensors_proton.resize(0); 
		tensors_neutron.resize(0);
		nTensors = number_tensors_tot[irecord];
//	create list of all proton {a+a}_p and neutron {a+a}_n labels that are
//	stored in a given interaction file with {n1, n2, 0, n4, n5, 0}
		for (size_t itensor = 0; itensor < nTensors; ++itensor)
		{
			adta_p.rho = 1;
			adta_p.lm  = quantum_labels_tot[quantum_index++];
			adta_p.mu  = quantum_labels_tot[quantum_index++];
			adta_p.S2  = quantum_labels_tot[quantum_index++];

			adta_n.rho = 1;
			adta_n.lm  = quantum_labels_tot[quantum_index++];
			adta_n.mu  = quantum_labels_tot[quantum_index++];
			adta_n.S2  = quantum_labels_tot[quantum_index++];

			quantum_index += 4;	//Important: we need to skip quantum labels of the proton-neutron coupling!

			tensors_proton.push_back(adta_p);
			tensors_neutron.push_back(adta_n);
		}
//	not all the tensors labels are unique ==> std::sort and use resize & unique
//	to squeeze non unique labels out of tensor_proton and tensor_neutron
		std::sort(tensors_proton.begin(), tensors_proton.end());
		std::sort(tensors_neutron.begin(), tensors_neutron.end());
		tensors_proton.resize(std::unique(tensors_proton.begin(), tensors_proton.end()) - tensors_proton.begin());
		tensors_neutron.resize(std::unique(tensors_neutron.begin(), tensors_neutron.end()) - tensors_neutron.begin());

		std::map<CTuple<char, 2>, SU3xSU2_VEC>::iterator it = structure_with_irreps.find(structure_proton);
//	is {a+a} [as described by structure_proton] already in the map ?
		if (it == structure_with_irreps.end()) // ==> add it into map together with all the proton_irreps
		{
			structure_with_irreps.insert(std::make_pair(structure_proton, tensors_proton));
		}
		else	//	a+a is already stored ==> for each tensor {a+a}_p label in tensor_proton test whether it is already stored in structure_with_irreps[na nb] 
		{
			SU3xSU2_VEC& stored_array = it->second;	//	stored_array == vector of SU3xSU2 irreps stored in map
			for (size_t itensor = 0; itensor < tensors_proton.size(); ++itensor)
			{	
				//	store irreps in tensor_proton only if they are not stored in memory in structure stored_array
				if (std::find(stored_array.begin(), stored_array.end(), tensors_proton[itensor]) == stored_array.end())
				{
					stored_array.push_back(tensors_proton[itensor]);
				}
			}
		}
//	is {a+a} [as described by structure_neutron] already in the map ?
		it = structure_with_irreps.find(structure_neutron);
		if (it == structure_with_irreps.end())
		{
			structure_with_irreps.insert(std::make_pair(structure_neutron, tensors_neutron));
		}
		else	//	a+a is already stored ==> for each tensor in tensor_neutron test whether it is also stored
		{
			SU3xSU2_VEC& stored_array = it->second;	//	stored_array == vector of SU3xSU2 irreps stored in map
			for (size_t itensor = 0; itensor < tensors_neutron.size(); ++itensor)
			{	
				if (std::find(stored_array.begin(), stored_array.end(), tensors_neutron[itensor]) == stored_array.end())
				{
					stored_array.push_back(tensors_neutron[itensor]);
				}
			}
		}
	}

	if (my_rank == 0)
	{
		duration = boost::chrono::system_clock::now() - start;
		cout << "CInteractionPN::AddOperators\t\tLoop #1 ... " << duration << endl;
	}
	start = boost::chrono::system_clock::now();


//	
//	Step 2:	use data stored in structure_with_irreps to create a vector of CTensorGroup_ada* whose elements are then filled with CTensorStructure* 	
//
	unsigned int firstLabelId(0);

//	{na nb} ---> {ir1, ir2, ir3, ...}
	std::map<CTuple<char, 2>, SU3xSU2_VEC>::iterator it = structure_with_irreps.begin();
	for (; it != structure_with_irreps.end(); ++it)
	{
		SU3xSU2_VEC tensorLabels(1);	//	a+a is a one-body interaction ==> only one SU3xSU2 label is needed 
		std::vector<char> structure(3);
		std::vector<std::vector<char> > vSingleShell_structures;	//	array of structures describing construction of the single-shell tensors
		std::vector<unsigned char> ShellsT;							//	array of HO shells of all single-shell tensors present in a general tensor
		char nShellsT;												//	number of shells in tensor => for one-body interaction max(nShellsT) = 2
		std::vector<std::vector<SU3xSU2::LABELS*> > vSingleShellTensorLabels; 	
		std::vector<SU3xSU2::LABELS*> vOmega;

		structure[0] = it->first[0];
		structure[1] = it->first[1];
		structure[2] = 0;
	
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
		StructureToSingleShellTensors(structure, ShellsT, tensorLabels, vSingleShell_structures, vSingleShellTensorLabels, vOmega);
		nShellsT = ShellsT.size();
		assert(nShellsT == 1 || nShellsT == 2);

		std::vector<char> dA(ShellsT.size(), 0);				//	case (A): dA = {+1, -1} or {-1, +1}; case (B): dA = {0}
		std::vector<unsigned char> nAnnihilators(nShellsT, 0);	//	case (A): nAnnihilators = {0, 1} or {1, 0}; case (B): nAnnihilators = {1}
		CTensorGroup_ada* pCTensorGroup(NULL);  

		for (int i = 0; i < nShellsT; ++i)
		{
			nAnnihilators[i] = std::count_if(vSingleShell_structures[i].begin(), vSingleShell_structures[i].end(), std::bind2nd(std::less<char>(), 0));
			dA[i] = std::count_if(vSingleShell_structures[i].begin(), vSingleShell_structures[i].end(), std::bind2nd(std::greater<char>(), 0)) // #creations - #annihilations
					- nAnnihilators[i];
		}
		assert(std::accumulate(dA.begin(), dA.end(), 0) == 0);

		//	use structure, ShellsT, dA, nAnnihilators to create class CTensorGroup_ada
		tensorGroups_.push_back(new CTensorGroup_ada(structure, ShellsT, dA, nAnnihilators, firstLabelId));
		pCTensorGroup = tensorGroups_.back();

		nTensors = it->second.size();
//	At this moment, vectors of labels from structure_with_irreps are not
//	sorted!!!!  It is absolutly imperative to sort them, otherwise we could not
//	do a binary search which is needed for fast step #3 implementation
		std::sort(it->second.begin(), it->second.end());
		pCTensorGroup->m_Tensors.reserve(nTensors);
		for (size_t itensor = 0; itensor < nTensors; ++itensor)
		{
			//	either vOmega has one element ---> &tensorLabels[0]	[if na != nb] && vSingleShell_structures = {na irrep, nb irrep}
			//	or vSingleShell_structures has one element ----> &tensorLabels[0]	[if na == nb] && vOmega={empty}
			tensorLabels[0] = it->second[itensor];
			std::vector<CrmeTable*> SingleShellTensorRmeTables(ShellsT.size(), NULL); // number of single-shell tensors = number of active shells
			//	Generate a vector of pointers to single-shell tensors rme tables
			for (int j = 0; j < ShellsT.size(); ++j)
			{
//	use global class rmeLookUpTableContainer_ to pointer to rme for a given single shell tensor
				SingleShellTensorRmeTables[j] = rmeLookUpTableContainer_.GetRMETablePointer(ShellsT[j], nAnnihilators[j], dA[j], vSingleShell_structures[j], vSingleShellTensorLabels[j], generate_missing_rme_files(), my_rank == 0, std::cout);
				assert(SingleShellTensorRmeTables[j]); //	since GetRMETablePointer should never return NULL !
			}
			pCTensorGroup->AddTensor(SingleShellTensorRmeTables, vOmega);
		}
		firstLabelId += nTensors;
	}
//	sort CTensorGroup_ada according to structure_ arrays ==> I am able to use binary search to find a particular structure	
	std::sort(tensorGroups_.begin(), tensorGroups_.end(), CTensorGroup_adaLess());

	if (my_rank == 0)
	{
		duration = boost::chrono::system_clock::now() - start;
		cout << "CInteractionPN::AddOperators\t\tLoop #2 ... " << duration << endl;
	}
	start = boost::chrono::system_clock::now();



//	Step 3:	
//	read input files and create a map of <index_proton, index_neutron>
//	associated to tensor strenght {a1, .... } and tensor labels {IR0, IR0',	...} 
//	at the very end copy the content of this map into the vector whose
//	elements are <index_proton, index_neutron>, {IR0, IR0', ... }, and {a1, ... } 
//	sorted according to <index_proton, index_neutron> pairs
//
	uint8_t irho, lm, mu, S2;
	uint16_t nCoeffs;

	std::map<std::pair<unsigned int, unsigned int>, std::pair<std::vector<TENSOR_STRENGTH>, SU3xSU2_VEC> > unique_pnInteractionData;
	std::map<std::pair<unsigned int, unsigned int>, std::pair<std::vector<TENSOR_STRENGTH>, SU3xSU2_VEC> >::iterator pnTensor;
	
	unsigned int index_proton, index_neutron;
	
	CTensorGroup_ada* proton_tensor_group;
	CTensorGroup_ada* neutron_tensor_group;

	SU3xSU2_VEC* proton_labels;
	SU3xSU2_VEC* neutron_labels;
	std::vector<TENSOR_STRENGTH> coeffs;
	coeffs.reserve(256); // allocate buffer large enough for any (l0 m0)S0 irrep

	for (uint32_t irecord = 0, quantum_index = 0, coeffs_index = 0; irecord < number_records; ++irecord)
	{
		int8_t* structure = &structures_tot[4*irecord];

		structure_proton[0]  = structure[0];
		structure_proton[1]  = structure[1];

		structure_neutron[0] = structure[2]; 
		structure_neutron[1] = structure[3];

//	find CTensorGroup_ada that has structure_ == structure_proton				
		std::vector<CTensorGroup_ada*>::iterator it = std::lower_bound(tensorGroups_.begin(), tensorGroups_.end(), structure_proton, CTensorGroup_adaLess());
		assert(it != tensorGroups_.end() && (*it)->structure_ == structure_proton);
//	proton_tensor_group_pointer = pointer to CTensorGroup_ada that has structure_ == structure_proton				
		proton_tensor_group = *it;
//	proton_labels == pointer to SU3xSU2_VEC associated with structure_proton
		proton_labels = &((structure_with_irreps.find(structure_proton))->second);

//	the same steps for neutrons
		it = std::lower_bound(tensorGroups_.begin(), tensorGroups_.end(), structure_neutron, CTensorGroup_adaLess());
		assert(it != tensorGroups_.end() && (*it)->structure_ == structure_neutron);
		neutron_tensor_group = *it;
		neutron_labels = &((structure_with_irreps.find(structure_neutron))->second);

		nTensors = number_tensors_tot[irecord];
//	initialize indices to the highest possible value
		for (size_t itensor = 0; itensor < nTensors; ++itensor)
		{
//	read proton tensor SU3xSU2 labels				
			lm   = quantum_labels_tot[quantum_index++];
			mu   = quantum_labels_tot[quantum_index++];
			S2   = quantum_labels_tot[quantum_index++];
//	if we are about to store a given set of proton-neutron tensors ==>
//	calculate index for a given {a+a}_{proton} irho(lm mu)S tensor
//	index_proton <-- index, where index is such that (*proton_labels)[index] == (lm mu)S2
			index_proton = std::lower_bound(proton_labels->begin(), proton_labels->end(), SU3xSU2::LABELS(1, lm, mu, S2)) - proton_labels->begin();
//	to get the total index_proton one needs to the first Id of the first tensor label					
			index_proton += proton_tensor_group->firstLabelId_;
//	same steps for neutron ...				
			lm   = quantum_labels_tot[quantum_index++];
			mu   = quantum_labels_tot[quantum_index++];
			S2   = quantum_labels_tot[quantum_index++];
		
			index_neutron  = std::lower_bound(neutron_labels->begin(), neutron_labels->end(), SU3xSU2::LABELS(1, lm, mu, S2)) - neutron_labels->begin();
			index_neutron += neutron_tensor_group->firstLabelId_;

			irho = quantum_labels_tot[quantum_index++];
			lm   = quantum_labels_tot[quantum_index++];
			mu   = quantum_labels_tot[quantum_index++];
			S2   = quantum_labels_tot[quantum_index++];
			SU3xSU2::LABELS ir0(irho, lm, mu, S2);

//			nCoeffs = irho*SU3::kmax(ir0, S2/2);
			nCoeffs = irho*1;

			coeffs.resize(0);
			coeffs.insert(coeffs.end(), &coeffs_tot[coeffs_index], &coeffs_tot[coeffs_index] + nCoeffs);
			coeffs_index += nCoeffs;

			pnTensor = unique_pnInteractionData.find(std::make_pair(index_proton, index_neutron));
			if (pnTensor == unique_pnInteractionData.end())
			{
				unique_pnInteractionData.insert(std::make_pair(std::make_pair(index_proton, index_neutron), std::make_pair(coeffs, SU3xSU2_VEC(1, ir0))));
			}
			else
			{
				pnTensor->second.first.insert(pnTensor->second.first.end(), coeffs.begin(), coeffs.end());
				pnTensor->second.second.push_back(ir0);
			}
		}
	}

	if (my_rank == 0)
	{
		duration = boost::chrono::system_clock::now() - start;
		cout << "CInteractionPN::AddOperators\t\tLoop #3 ... " << duration << endl;
	}

	pnInteraction_.reserve(unique_pnInteractionData.size());
	pnTensor = unique_pnInteractionData.begin();
	for (; pnTensor != unique_pnInteractionData.end(); ++pnTensor)
	{
		pnInteraction_.push_back(PNInteractionData(pnTensor->first, pnTensor->second.first, pnTensor->second.second));
	}
}

void CInteractionPN::GetTensors(std::vector<CTensorStructure*>& tensors) const
{
	for (size_t itensor_group = 0; itensor_group < tensorGroups_.size(); ++itensor_group)
	{
		tensorGroups_[itensor_group]->GetTensors(tensors);
	}
}

CInteractionPN::~CInteractionPN()
{
	for (size_t j = 0; j < tensorGroups_.size(); ++j)
	{
		delete tensorGroups_[j];
	}
}


void CInteractionPN::ShowTensorGroupsAndTheirTensorStructures() const
{
	for (size_t itensor_group = 0; itensor_group < tensorGroups_.size(); ++itensor_group)
	{
		tensorGroups_[itensor_group]->ShowTensorStructure();
		tensorGroups_[itensor_group]->ShowTensors();
	}
}

void CInteractionPN::ShowInteractionStructures() const
{
	for (std::vector<PNInteractionData>::const_iterator cit = pnInteraction_.begin(); cit != pnInteraction_.end(); ++cit)
	{
		std::cout << "(" << cit->id_.first << " " << cit->id_.second << ")\t{";
		
		for (size_t i = 0; i < cit->ncoeffs_ - 1; ++i)
		{
			std::cout << cit->coeffs_[i] << ", ";
		}
		std::cout << cit->coeffs_[cit->ncoeffs_ - 1] << "}\t{";

		for (size_t i = 0; i < cit->nTensors_ - 1; ++i)
		{
			std::cout << cit->tensor_labels_[i] << ", ";
		}
		std::cout << cit->tensor_labels_[cit->nTensors_ - 1] << "}\t{";

		for (size_t i = 0; i < cit->nTensors_ - 1; ++i)
		{
			std::cout << cit->su3so3cgTables_[i] << ", ";
		}
		std::cout << cit->su3so3cgTables_[cit->nTensors_ - 1] << "}" << std::endl;
	}
}


size_t CInteractionPN::GetTensorsCoefficientsLabels(const unsigned int index_proton, const unsigned int index_neutron, 
													SU3xSU2::LABELS*& tensor_labels, WigEckSU3SO3CGTable**& su3so3cgTables, TENSOR_STRENGTH*& coeffs) const
{
	std::vector<PNInteractionData>::const_iterator it = std::lower_bound(pnInteraction_.begin(), pnInteraction_.end(), std::make_pair(index_proton, index_neutron));
	if (it == pnInteraction_.end() || (it->id_.first != index_proton || it->id_.second != index_neutron))
	{
		return 0;
	}
	else
	{
		tensor_labels = it->tensor_labels_;
		su3so3cgTables = it->su3so3cgTables_;
		coeffs = it->coeffs_;
		return it->nTensors_;
	}
}

size_t CInteractionPN::GetRelevantTensorGroups(	const std::vector<unsigned char>& Bra_confs, 
												const std::vector<unsigned char>& Ket_confs, 
												const std::vector<unsigned char>& hoShells, 
												std::vector<CTensorGroup_ada*>& tensorGroups, 
												std::vector<int>& phase_tensor_group) const
{
	std::vector<unsigned char> tensor(hoShells.size());
	size_t i, j;
	int dA = 0;
//	two-body interation => sum |bra_confs - ket_confs| <= 2	
	for (i = 0; i < Bra_confs.size(); ++i)
	{
		dA += abs(Bra_confs[i] - Ket_confs[i]);
		if (dA > 2)
		{
			assert(tensorGroups.empty());
			return tensorGroups.size(); //	==> 2B interaction could not couple
		}
	}
//	==> bra and ket configurations differ at most on 2 shells.  Right now I
//	trivially iterate over all families of TensorGroups, which is not very
//	economical. 
//
//	Iterate over all families of TensorGroups (.i.e. SingleShell, TwoShellsdAZero, TwoShells, ThreeShells, FourShells)
	for (j = 0; j < tensorGroups_.size(); ++j)	// iterate over all TensorGroups in the family ...
	{
		CTensorGroup_ada* pTG = tensorGroups_[j];
		if (pTG->Couple(Bra_confs, Ket_confs, hoShells)) // Does TensorGroup connet bra and ket configurations ?
		{
			tensor.assign(tensor.size(), 0);	//	tensor is an array whose element [i] hold the number of a+/a operators in HO shell n == hoShells[i]
			const std::vector<char>& n_adta_InSST = pTG->n_adta_InSST(); 
			const std::vector<unsigned char>& ShellsT = pTG->ShellsT(); 
			for(size_t i = 0; i < hoShells.size(); ++i)	// iterate HO shells with non-zero number of fermions 
			{
				std::vector<unsigned char>::const_iterator cit = std::find(ShellsT.begin(), ShellsT.end(), hoShells[i]);	// *cit == HO shell n == hoShells[i]
				if (cit != ShellsT.end()) // there is a single-shell tensor acting on HO shell n == hoShells[i]
				{
					size_t index = cit - ShellsT.begin();
					tensor[i] = n_adta_InSST[index];	// store number of a+/a operators ==> needed for phase calculation
				}
			}
			phase_tensor_group.push_back(rme_uncoupling_phase(Bra_confs, tensor, Ket_confs));	//	calculate and store a phase due to the exchanges of fermion creation/annihilation operators [see relation (5.6)]
			tensorGroups.push_back(pTG);
		}
	}
	return tensorGroups.size();
}


//////////////////////////////////////////////////////////////////////////////////////	
/////////////////     P N I n t e r a c t i o n D a t a    ///////////////////////////	
//////////////////////////////////////////////////////////////////////////////////////	
PNInteractionData::PNInteractionData(const std::pair<unsigned int, unsigned int>& id, const std::vector<TENSOR_STRENGTH>& coeffs, const SU3xSU2_VEC& tensor_labels): id_(id)
{
	CWigEckSU3SO3CGTablesLookUpContainer wigEckSU3SO3CGTablesLookUpContainer;
	nTensors_ = tensor_labels.size();
	ncoeffs_ = coeffs.size();

	tensor_labels_ = new SU3xSU2::LABELS[nTensors_];
	su3so3cgTables_ = new WigEckSU3SO3CGTable*[nTensors_];
	coeffs_ = new TENSOR_STRENGTH[coeffs.size()];

	std::copy(tensor_labels.begin(), tensor_labels.end(), tensor_labels_);
	std::copy(coeffs.begin(), coeffs.end(), coeffs_);

	for (size_t i = 0; i < nTensors_; ++i)
	{
		su3so3cgTables_[i] = wigEckSU3SO3CGTablesLookUpContainer.GetSU3SO3CGTablePointer(tensor_labels_[i], tensor_labels_[i].S2);
	}
}

void PNInteractionData::operator=(const PNInteractionData& pndata)
{
	id_ = pndata.id_; 
	nTensors_ = pndata.nTensors_; 
	ncoeffs_ = pndata.ncoeffs_;

	tensor_labels_ = new SU3xSU2::LABELS[nTensors_];
	coeffs_ = new TENSOR_STRENGTH[ncoeffs_];
	su3so3cgTables_ = new WigEckSU3SO3CGTable*[nTensors_];

	std::copy(pndata.tensor_labels_, pndata.tensor_labels_ + nTensors_, tensor_labels_);
	std::copy(pndata.coeffs_, pndata.coeffs_ + ncoeffs_, coeffs_);
	std::copy(pndata.su3so3cgTables_, pndata.su3so3cgTables_ + nTensors_, su3so3cgTables_);
}

PNInteractionData::PNInteractionData(const PNInteractionData& pndata):id_(pndata.id_), nTensors_(pndata.nTensors_), ncoeffs_(pndata.ncoeffs_)
{
	tensor_labels_ = new SU3xSU2::LABELS[nTensors_];
	su3so3cgTables_ = new WigEckSU3SO3CGTable*[nTensors_];
	coeffs_ = new TENSOR_STRENGTH[ncoeffs_];

	std::copy(pndata.tensor_labels_, pndata.tensor_labels_ + nTensors_, tensor_labels_);
	std::copy(pndata.coeffs_, pndata.coeffs_ + ncoeffs_, coeffs_);
	std::copy(pndata.su3so3cgTables_, pndata.su3so3cgTables_ + nTensors_, su3so3cgTables_);
}


