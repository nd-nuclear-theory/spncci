#include <SU3ME/SU3InteractionRecoupler.h>
//TESTED
namespace recoupler_nonscalar

{
struct CNegligible
{
	const float eps_;
	CNegligible(const float eps): eps_(eps) {assert(eps_ > 0.0);}
	inline bool operator() (const double d) const
	{
		return fabs(d) < eps_; // true ==> d is negligible
	}
};

void SU3InteractionRecoupler::RemoveTensorsWithAmplitudesLess(const float eps)
{
	CNegligible IsNegligible(eps);
	for (int Type = TENSOR::PPNN; Type <= TENSOR::PN; ++Type)
	{
		for (size_t n = 0; n <= 3; ++n)
		{
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>& nShellsTensors = (Type == TENSOR::PN) ? m_PN[n] : m_PPNN[n];
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::iterator it = nShellsTensors.begin();
	
			if (it == nShellsTensors.end()) {
				continue;
			}
			for (; it != nShellsTensors.end(); ++it)
			{
				TENSOR_COMPONENTS& TensorComponents = it->second;
				//	TensorsComponents is a map of TENSOR_LABELS and associated
				//	coefficients. Since I am erasing elements of this map,
				//	position of the end of map could change. Therefore, i need
				//	to call tensorsComponents.end() to make sure that it has
				//	the correct value 
				for (TENSOR_COMPONENTS::iterator itTensor = TensorComponents.begin(); itTensor != TensorComponents.end();)
				{
					const std::vector<double>& Coeffs = itTensor->second;
					size_t nCoeffs = itTensor->first.CoeffsSize((Type == TENSOR::PN) ? TENSOR::PN : TENSOR::PPNN);
					
					assert(nCoeffs == Coeffs.size());

					if (std::count_if(Coeffs.begin(), Coeffs.end(), IsNegligible) == nCoeffs) // if |coeff| < 0.00001 ==> erase
					{
						it->second.erase(itTensor++);
					}
					else
					{
						++itTensor;
					}
				}
			}
		}
	}

	for (int Type = TENSOR::PPNN; Type <= TENSOR::PN; ++Type)
	{
		for (size_t n = 0; n <= 3; ++n)
		{
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>& nShellsTensors = (Type == TENSOR::PN) ? m_PN[n] : m_PPNN[n];
			for (std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::iterator it = nShellsTensors.begin(); it != nShellsTensors.end();)
			{
				if (it->second.empty())
				{
					nShellsTensors.erase(it++);
				}
				else
				{
					++it;
				}	
			}
		}
	}
}

void SU3InteractionRecoupler::PrintInfo()
{
	for (int Type = TENSOR::PPNN; Type <= TENSOR::PN; ++Type)
	{
		if (Type == TENSOR::PN) 
		{
			std::cout << "PN interaction:" << std::endl;
		} 
		else 
		{
			std::cout << "PPNN interaction:" << std::endl;
		}
		size_t nCTuples[4] = {0, 0, 0, 0};
		size_t nTensors[4] = {0, 0, 0, 0};

		for (size_t n = 0; n <= 3; ++n)
		{
			nCTuples[n] = (Type == TENSOR::PN) ? m_PN[n].size() : m_PPNN[n].size();

			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator cit = (Type == TENSOR::PN) ? m_PN[n].begin() : m_PPNN[n].begin();
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator LAST =(Type == TENSOR::PN) ? m_PN[n].end() : m_PPNN[n].end();

			for (; cit != LAST; ++cit)
			{
				nTensors[n] += cit->second.size();
			}
		}
		for (size_t n = 0; n <= 3; ++n)
		{
			std::cout << n+1 << " shells:" << std::endl;
			std::cout << "\tnumber of n1 n2 n3 n4: \t" << nCTuples[n] << std::endl;
			std::cout << "\tnumber of all tensors: \t" << nTensors[n] << std::endl;
			std::cout << std::endl;
		}
	}
}



void SU3InteractionRecoupler::Save(const std::string& sOutputFileName, const TENSOR::Type Type)
{
	std::ofstream Outputfile(sOutputFileName.c_str());

	if (!Outputfile)
	{
		std::cerr << "Can not open file " << sOutputFileName << " for the output!" << std::endl;
		exit(EXIT_FAILURE);
	}	

	Outputfile.precision(10);
	
	for (size_t n = 0; n <= 3; ++n)
	{
		std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator cit = (Type == TENSOR::PN) ? m_PN[n].begin() : m_PPNN[n].begin();
		std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator LAST =(Type == TENSOR::PN) ? m_PN[n].end() : m_PPNN[n].end();
		for (; cit != LAST; ++cit)
		{
			CTuple<char, 6> structure(cit->first);

			for (size_t i = 0; i < 6; ++i)
			{
				Outputfile << (int)structure[i] << " ";
			}
			Outputfile << std::endl;
			Outputfile << cit->second.size() << std::endl;

			TENSOR_COMPONENTS::const_iterator citTensor = cit->second.begin();
			TENSOR_COMPONENTS::const_iterator LAST_TENSOR = cit->second.end();
			for (; citTensor != LAST_TENSOR; ++citTensor)
			{
				const TENSOR_LABELS& Labels = citTensor->first;
				Outputfile << Labels.IR1 << "\t" << Labels.IR2 << "\t" << Labels.IR0 << std::endl;
				const std::vector<double>& Coeffs = citTensor->second;
				int rho0max = Labels.IR0.rho;
		//		int k0max   = SU3::kmax(Labels.IR0, Labels.IR0.S2/2);
				int k0max   = 1;
				size_t nCoeffs = Labels.CoeffsSize(Type);
				if (Type == TENSOR::PN)
				{
					for (size_t k0 = 0, i = 0; k0 < k0max; ++k0)
					{
						for(size_t rho0 = 0; rho0 < rho0max; ++rho0)
						{
							if (Negligible6(Coeffs[i]))
							{
								Outputfile << 0 << " " << std::endl;
							}
							else
							{
								Outputfile << Coeffs[i] << " " << std::endl;
							}
							i++;
						}
					}
				}
				else
				{
					for (size_t k0 = 0, i = 0; k0 < k0max; ++k0)
					{
						for(size_t rho0 = 0; rho0 < rho0max; ++rho0, i += 2)
						{
							if (Negligible6(Coeffs[i + PP]))
							{
								Outputfile << 0 << " ";
							}
							else
							{
								Outputfile << Coeffs[i + PP] << " "; 
							}

							if (Negligible6(Coeffs[i + NN]))
							{
								Outputfile << 0 << " " << std::endl;
							}
							else
							{
								Outputfile << Coeffs[i + NN] << " " << std::endl;
							}
						}
					}
				}
			}
		}
	}
}

void SU3InteractionRecoupler::Show2()
{
	for (int Type = TENSOR::PPNN; Type <= TENSOR::PN; ++Type)
	{
		if (Type == TENSOR::PN) 
		{
			std::cout << "PN interaction:" << std::endl;
		} 
		else 
		{
			std::cout << "PPNN interaction:" << std::endl;
		}

		for (size_t n = 0; n <= 3; ++n)
		{
			std::cout << n+1 << " shells:" << std::endl;
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator cit = (Type == TENSOR::PN) ? m_PN[n].begin() : m_PPNN[n].begin();
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator LAST =(Type == TENSOR::PN) ? m_PN[n].end() : m_PPNN[n].end();
	
			if (cit == LAST) 
			{
				std::cout << "\t\tEmpty" << std::endl;
				continue;
			}
//	std::cout << m_Data[TensorType].size() << std::endl;
			for (; cit != LAST; ++cit)
			{
				CTuple<char, 6> structure(cit->first);
				std::cout << "\t";
				for (size_t i = 0; i < 6; ++i)
				{
					std::cout << (int)structure[i] << " ";
				}
				std::cout << std::endl;

				TENSOR_COMPONENTS::const_iterator citTensor = cit->second.begin();
				TENSOR_COMPONENTS::const_iterator LAST_TENSOR = cit->second.end();
				for (; citTensor != LAST_TENSOR; ++citTensor)
				{
					const TENSOR_LABELS& Labels = citTensor->first;
					std::cout << "\t\t";
					std::cout << Labels.IR1 << " " << Labels.IR2 << " " << Labels.IR0 << "\t";

					const std::vector<double>& Coeffs = citTensor->second;
					int rho0max = Labels.IR0.rho;
//					int k0max   = SU3::kmax(Labels.IR0, Labels.IR0.S2/2);
					int k0max   = 1;
					size_t nCoeffs = Labels.CoeffsSize((Type == TENSOR::PN) ? TENSOR::PN : TENSOR::PPNN);
					if (Type == TENSOR::PN)
					{
						for (size_t k0 = 0, i = 0; k0 < k0max; ++k0)
						{
							for(size_t rho0 = 0; rho0 < rho0max; ++rho0)
							{
								std::cout << Coeffs[i++] << " ";
							}
						}
						std::cout << std::endl;
					}
					else
					{
						for (size_t k0 = 0, i = 0; k0 < k0max; ++k0)
						{
							for(size_t rho0 = 0; rho0 < rho0max; ++rho0, i += 2)
							{
								std::cout << Coeffs[i+PP] << " ";
							}
						}
						std::cout << std::endl;
					}
				}
			}
		}
	}
}


void SU3InteractionRecoupler::Show()
{
	for (int Type = TENSOR::PPNN; Type <= TENSOR::PN; ++Type)
	{
		if (Type == TENSOR::PN) 
		{
			std::cout << "PN interaction:" << std::endl;
		} 
		else 
		{
			std::cout << "PPNN interaction:" << std::endl;
		}

		for (size_t n = 0; n <= 3; ++n)
		{
			std::cout << n+1 << " shells:" << std::endl;
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator cit = (Type == TENSOR::PN) ? m_PN[n].begin() : m_PPNN[n].begin();
			std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::const_iterator LAST =(Type == TENSOR::PN) ? m_PN[n].end() : m_PPNN[n].end();
	
			if (cit == LAST) 
			{
				std::cout << "\t\tEmpty" << std::endl;
				continue;
			}
//	std::cout << m_Data[TensorType].size() << std::endl;
			for (; cit != LAST; ++cit)
			{
				CTuple<char, 6> structure(cit->first);
				TENSOR_COMPONENTS::const_iterator citTensor = cit->second.begin();
				TENSOR_COMPONENTS::const_iterator LAST_TENSOR = cit->second.end();
				for (; citTensor != LAST_TENSOR; ++citTensor)
				{
					const TENSOR_LABELS& Labels = citTensor->first;
					for (size_t i = 0; i < 6; ++i)
					{
						std::cout << (int)structure[i] << " ";
					}
					std::cout << "\t";
					std::cout << Labels.IR1 << " " << Labels.IR2 << " " << Labels.IR0 << std::endl;

					const std::vector<double>& Coeffs = citTensor->second;
					int rho0max = Labels.IR0.rho;
//					int k0max   = SU3::kmax(Labels.IR0, Labels.IR0.S2/2);
					int k0max   = 1;
					size_t nCoeffs = Labels.CoeffsSize((Type == TENSOR::PN) ? TENSOR::PN : TENSOR::PPNN);
					if (Type == TENSOR::PN)
					{
						for (size_t k0 = 0, i = 0; k0 < k0max; ++k0)
						{
							for(size_t rho0 = 0; rho0 < rho0max; ++rho0)
							{
								std::cout << "k0 = " << k0 << " rho0 = " << rho0 << "\t\t" << Coeffs[i++] << std::endl;
							}
						}
					}
					else
					{
						for (size_t k0 = 0, i = 0; k0 < k0max; ++k0)
						{
							for(size_t rho0 = 0; rho0 < rho0max; ++rho0, i += 2)
							{
								std::cout << "k0 = " << k0 << " rho0 = " << rho0 << "\t\t" << Coeffs[i+PP] << " " << Coeffs[i + NN] << std::endl;
							}
						}
					}
				}
			}
		}
	}
}

//TESTED
void SU3InteractionRecoupler::Add(const TENSOR::Type Type, const char n, const char* structure, const TENSOR_LABELS& TensorLabels, const std::vector<double>& TensorCoeffs)
{
//	if (std::count_if(TensorCoeffs.begin(), TensorCoeffs.end(), Negligible) == TensorCoeffs.size())
//	{
//		return;
//	}

	std::map<CTuple<char, 6>, TENSOR_COMPONENTS>& nShellsTensors = (Type == TENSOR::PN) ? m_PN[n-1] : m_PPNN[n-1];
	std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::iterator it 	=	nShellsTensors.find(structure);
	std::map<CTuple<char, 6>, TENSOR_COMPONENTS>::iterator LAST =	nShellsTensors.end();

	if (it == LAST)
	{
		TENSOR_COMPONENTS m;
		m.insert(make_pair(TensorLabels, TensorCoeffs));
		nShellsTensors.insert(make_pair(CTuple<char, 6>(structure), m));
	}
	else
	{
		TENSOR_COMPONENTS::iterator itTensor = (it->second).find(TensorLabels);
		if (itTensor == it->second.end())
		{
			(it->second).insert(make_pair(TensorLabels, TensorCoeffs));
		} 
		else 
		{
			int nCoeffs = TensorLabels.CoeffsSize(Type);
			std::vector<double>& NewCoeffs = itTensor->second;
			
			assert(nCoeffs);	// if nCoeffs == 0 ==> L0 does not exist in the basis of TensorLabels.IR0 
			assert(nCoeffs == NewCoeffs.size()); 

			for (size_t i = 0; i < nCoeffs; ++i)
			{
				NewCoeffs[i] += TensorCoeffs[i];
			}
		}
	}
}

//	Implement recoupling formula [{U x V} x {X x Z}] ---> [{U x Z} x {V x X}]  ... see (A.57)
void SU3InteractionRecoupler::UV_XZtoUZ_VX(
    const TENSOR::Type Type, const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2,
    const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4, const TENSOR_LABELS& TensorLabels,
    const std::vector<double>& TensorCoeffs,
    std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorComponents) {
/*
//	This sequence of commands performs A.57 in two steps: 
//	#1)A.52 #2) A.56

	std::vector<aPP_aNN_aPN> Coeffs(TensorCoeffs);
	UV_XZtoUV_ZX(ir3, ir4, TensorLabel, n1n2n3n4New, Coeffs);
	UV_XZtoUX_VZ(ir1, ir2, ir3, ir4, TensorLabel, Coeffs, n1n2n3n4New, TTY1PNTensorComponents);
*/	
	static double dExchangePhase = +1.0;


	const SU3xSU2::LABELS& ir12 = TensorLabels.IR1;
	const SU3xSU2::LABELS& ir34 = TensorLabels.IR2;
	const SU3xSU2::LABELS& ir0 	= TensorLabels.IR0;
	const SU2::LABEL& S12 = ir12.S2;
	const SU2::LABEL& S34 = ir34.S2;
	const SU2::LABEL& S0  = ir0.S2;

	SU2::LABEL S14, S23;

	double dPhase = dExchangePhase*MINUSto((ir3.lm + ir3.mu) + (ir4.lm + ir4.mu) - (ir34.lm + ir34.mu) + 1 - S34/2)*sqrt((double)((S12 + 1)*(S34 + 1)));

	TENSOR_LABELS NewLabel;
	NewLabel.IR0 = TensorLabels.IR0;

	int rho1234max = abs(ir0.rho); // ir0.rho could be equal to -1 in case it was constructed with SU3xSU2::LABELS(lm, mu, S2)
	int rho1423max;

	double dsu2;
	size_t i, irho0p, irho0;
//	size_t k0, k0max = SU3::kmax(TensorLabels.IR0, TensorLabels.IR0.S2/2);
	size_t k0, k0max = 1;
	size_t index, index_k0_rho0p, index_k0_rho0;

	std::vector<SU3::LABELS>::const_iterator ir14, LAST_14_IRREP;
	std::vector<SU3::LABELS>::const_iterator ir23, LAST_23_IRREP;

	std::vector<SU3::LABELS> Irreps14;	
	std::vector<SU3::LABELS> Irreps23;	

	SU3::Couple(ir1, ir4, Irreps14);	// ir1 x ir4 -> ir14
	SU3::Couple(ir2, ir3, Irreps23); 	// ir2 x ir3 -> ir23

	LAST_14_IRREP = Irreps14.end();
	LAST_23_IRREP = Irreps23.end();
	for (ir14 = Irreps14.begin(); ir14 != LAST_14_IRREP; ++ir14)
	{
		NewLabel.IR1.rho = ir14->rho;
		NewLabel.IR1.lm  = ir14->lm;
		NewLabel.IR1.mu  = ir14->mu;
		for (ir23 = Irreps23.begin(); ir23 != LAST_23_IRREP; ++ir23)
		{
			rho1423max = SU3::mult(*ir14, *ir23, ir0); 	// 	Does ir14 x ir23 yields ir0 ?
			if (rho1423max == 0) 						// 	if NOT => continue;
			{
				continue;								
			}
// All tensors must have the same SU(3) label (lm0 mu0) = ir0. However,
// multiplicity is given by SU3::multu(ir14, ir23, ir0) = rho1423max ==
// rho0'_max
			NewLabel.IR0.rho = rho1423max;	

			NewLabel.IR2.rho = ir23->rho;
			NewLabel.IR2.lm  = ir23->lm;
			NewLabel.IR2.mu  = ir23->mu;
	
			size_t nCoeffs = NewLabel.CoeffsSize(Type);
			std::vector<double> NewCoeffs_SU3Part(nCoeffs, 0.0); 
			std::vector<double> NewCoeffs(nCoeffs); 
			std::vector<double> su39lm(rho1234max*rho1423max, 0.0);
			u9lm_look_up_table_.Get9lm(ir1, ir2, ir12, ir4, ir3, ir34, *ir14, *ir23, ir0, &su39lm[0]);

			NewCoeffs_SU3Part.assign(nCoeffs, 0.0); // set all coefficients to 0.0

			if (Type == TENSOR::PN)
			{
//				the two very first loops over k0 and rho0p enable us to
// 				calculate (lm14 mu14)S14, (lm23 mu23)S23 coefficients
// 				for all k0 and rho0p
				for (k0 = 0; k0 < k0max; ++k0)
				{
					for (irho0p = 0, index_k0_rho0p = k0*rho1423max; irho0p < rho1423max; ++irho0p, ++index_k0_rho0p)
					{ 
//						index_k0_rho0p = k0*rho1423max + irho0p; 
						for (irho0 = 0, index_k0_rho0  = k0*rho1234max, index = irho0p*rho1234max; irho0 < rho1234max; ++irho0, ++index_k0_rho0, ++index)
						{
//							index_k0_rho0  = k0*rho1234max + irho0;	
// 							index = irho0p*rho1234max + irho0;
//	index = rho12 + rho34*n1 + rho1234*n2 + rho14*n3 + rho23*n4 + rho1423*n5;
//	Since rho12 = rho34 = rho14 = rho23 = 0 and rho1234=rho0, rho1423=rho0',n2=1, n5=rho1234max, we have index = irho0 + irho0p*rho1234max
							NewCoeffs_SU3Part[index_k0_rho0p] += TensorCoeffs[index_k0_rho0] * su39lm[index];
						}
					}
				}
			}
			else
			{
				for (k0 = 0; k0 < k0max; ++k0)
				{
					for (irho0p = 0, index_k0_rho0p = k0*rho1423max*2; irho0p < rho1423max; ++irho0p, index_k0_rho0p += 2)
//					for (irho0p = 0; irho0p < rho1423max; ++irho0p)
					{ 
//						index_k0_rho0p = k0*rho1423max*2 + irho0p*2; 	
						for (irho0 = 0, index_k0_rho0 = k0*rho1234max*2, index = irho0p*rho1234max; irho0 < rho1234max; ++irho0, index_k0_rho0 += 2, ++index)
//						for (irho0 = 0; irho0 < rho1234max; ++irho0)
						{
//							index_k0_rho0  = k0*rho1234max*2 + irho0*2;	
//							index = irho0p*rho1234max + irho0;			
//	index = rho12 + rho34*n1 + rho1234*n2 + rho14*n3 + rho23*n4 + rho1423*n5;
//	Since rho12 = rho34 = rho14 = rho23 = 0 and rho1234=rho0, rho1423=rho0',n2=1, n5=rho1234max, we have index = irho0 + irho0p*rho1234max
							NewCoeffs_SU3Part[index_k0_rho0p + PP] += TensorCoeffs[index_k0_rho0 + PP] * su39lm[index];
							NewCoeffs_SU3Part[index_k0_rho0p + NN] += TensorCoeffs[index_k0_rho0 + NN] * su39lm[index];
						}
					}
				}
			}
			for (S14 = 0; S14 <= 2; S14 += 2)
			{
				NewLabel.IR1.S2 = S14;
				for (S23 = 0; S23 <= 2; S23 += 2)
				{
					if (SU2::mult(S14, S23, S0) == 0)
					{
						continue;
					}
					NewLabel.IR2.S2 = S23;
					dsu2 = dPhase * sqrt((double)((S14 + 1)*(S23 + 1)))*wigner9j(1, 1, S12, 1, 1, S34, S14, S23, S0);
//std::cout << "wigner9j(1, 1," <<  (int)S12 << ", 1, 1, " << (int)S34 <<", " << (int)S14 << ", " << (int)S23 << ", " << (int)S0 << ") = " << wigner9j(1, 1, S12, 1, 1, S34, S14, S23, S0) << std::endl;
					if (Negligible8(dsu2)) {
						continue;
					}
					for (i = 0; i < nCoeffs; ++i) 
					{
						NewCoeffs[i] = dsu2*NewCoeffs_SU3Part[i]; // multiply all the resulting coefficients by dsu2
					}
					TensorComponents.push_back(make_pair(NewLabel, NewCoeffs));
				}
			}
		}
	}
}

// this facilitates the two transformations: UV_XZ -> UV_ZX and UV_XZ -> VU_XZ, where IT IS ASSUMED THAT U, V, X, and Z ARE ALL COMPOSED BY THE SINGLE CREATION/ANNIHILATION OPERATOR
void SU3InteractionRecoupler::FlipTensors(const SU3::LABELS& irA, const SU3::LABELS& irB,  const SU3xSU2::LABELS& irAB, std::vector<double>& Coeffs)
{
	static double dExchangePhase = -1.0;
	double dPhase = dExchangePhase*MINUSto(irA.lm + irA.mu + irB.lm + irB.mu - irAB.lm - irAB.mu + 1 - irAB.S2/2);
	for (size_t i = 0; i < Coeffs.size(); ++i)
	{
		Coeffs[i] *= dPhase;
	}
}

void SU3InteractionRecoupler::UV_XZtoVU_XZ(char* structure, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ)
{
	SU3xSU2::LABELS ir1(structure[U]);
	SU3xSU2::LABELS ir2(structure[V]);
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::iterator it = TensorsUV_XZ.begin();
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::iterator LAST = TensorsUV_XZ.end();
	for (; it != LAST; ++it)
	{
		FlipTensors(ir1, ir2, it->first.IR1, it->second);
	}
	std::swap(structure[U], structure[V]);
}
void SU3InteractionRecoupler::UV_XZtoVU_XZ(char* structure, const SU3xSU2::LABELS& IR12, std::vector<double>& Coeffs)
{
	FlipTensors(SU3xSU2::LABELS(structure[U]), SU3xSU2::LABELS(structure[V]), IR12, Coeffs);
	std::swap(structure[U], structure[V]);
}

void SU3InteractionRecoupler::UV_XZtoUV_ZX(char* structure, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ)
{
	SU3xSU2::LABELS ir3(structure[X]);
	SU3xSU2::LABELS ir4(structure[Z]);
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::iterator it = TensorsUV_XZ.begin();
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::iterator LAST = TensorsUV_XZ.end();
	for (; it != LAST; ++it)
	{
		FlipTensors(ir3, ir4, it->first.IR2, it->second);
	}
	std::swap(structure[X], structure[Z]);
}
void SU3InteractionRecoupler::UV_XZtoUV_ZX(char* structure, const SU3xSU2::LABELS& IR34, std::vector<double>& Coeffs)
{
	FlipTensors(SU3xSU2::LABELS(structure[X]), SU3xSU2::LABELS(structure[Z]), IR34, Coeffs);
	std::swap(structure[X], structure[Z]);
}


//	(A.53)
void SU3InteractionRecoupler::UV_XZtoU__V_XZ(	const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4, 
												const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
												std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU__V_XZ)
{
	static double dExchangePhase = +1.0;

	CSU3CGMaster SU3Lib;

	int rho0max = TensorUV_XZ.IR0.rho;
	SU2::LABEL S12 = TensorUV_XZ.IR1.S2;
	SU2::LABEL S34 = TensorUV_XZ.IR2.S2;
	SU2::LABEL  S0 = TensorUV_XZ.IR0.S2;
	SU2::LABEL S234;

	const SU3xSU2::LABELS& ir12 = TensorUV_XZ.IR1;
	const SU3xSU2::LABELS& ir34 = TensorUV_XZ.IR2;
	const SU3xSU2::LABELS& ir0 	= TensorUV_XZ.IR0;

	TENSOR_LABELS NewLabel;
	NewLabel.IR1 = ir34;
	NewLabel.IR2.rho = 1; 	//	rho2_34max ... is equal to 1 due to the one-by-two coupling of creation/annihilation operators [Vx{XZ}], 
							//	since there are no multiplicities rho > 1 when coupling (lm mu) x (n 0) and (lm mu) x (0 n) irreps
	NewLabel.IR0.rho = 1; 	// 	rho0pmax ...   is equal to 1 due to the one-by-one-by-two coupling of creation/annihilation operators [Ux{Vx{XZ}}]
							//	since there are no multiplicities rho > 1 when coupling (lm mu) x (n 0) and (lm mu) x (0 n) irreps
	NewLabel.IR0.lm = ir0.lm;
	NewLabel.IR0.mu = ir0.mu;
	NewLabel.IR0.S2 = ir0.S2;


//	size_t k0, k0max = SU3::kmax(TensorUV_XZ.IR0, S0/2);
	size_t k0, k0max = 1;
	size_t nU6lm;

	double dPhase = dExchangePhase*MINUSto(1 + (S34 + S0)/2)*sqrt((double)(S12 + 1));
	double dSU2_S234;
	size_t irho0; 
	size_t i, index_k0_rho0, index_k0_rho0p;

	std::vector<SU3::LABELS> Irreps2_34;

	SU3::Couple(ir2, ir34, Irreps2_34);

	std::vector<SU3::LABELS>::const_iterator ir2_34 = Irreps2_34.begin();
	std::vector<SU3::LABELS>::const_iterator LAST_2_34_IRREP = Irreps2_34.end();
	for (; ir2_34 != LAST_2_34_IRREP; ++ir2_34)
	{
//	U[ ir1 ir2  ir0 ir3;  ir12 rho12   rho12_3   ir23  rho23  rho1_23 ]
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |     |     |       |         |  rho2_34 rho0p                     
//      |	|	 |	 |     |     |       |         |     |      | 
//      |	|	 |	 |     |     |       |         |     |      | 
//      V   V    V   V     V     V       V         V     V      V 
//	U[ ir1 ir2  ir0 ir34; ir12   1     rho0     ir2_34   1      1  ]

//		rho1_23max 	= SU3::mult(ir1,  ir23,   ir1_23);
//			|					 |		|  	    |			
//			|					 |		|	    |
//			V					 V		V	    V
		int rho0pmax = SU3::mult(ir1, *ir2_34,  ir0);
		if (rho0pmax == 0) {
			continue;
		}
		assert(rho0pmax == 1); 		//	if rho0pmax > 1 ==> SU3::mult does not work properly!!!
		assert(ir2_34->rho == 1);	//	and similarly, rho2_34max == 1

		NewLabel.IR2.lm = ir2_34->lm;
		NewLabel.IR2.mu = ir2_34->mu;

//		nU6lm = rho12max * rho12_3max * rho23max   * rho1_23max;
//				   |             |          |          |
//				   |             |          |          |
//				   |             |       rho2_34max rho0pmax 
//				   |             |          |          |
//				   |             |          |          |
//				   V             V          V          V
		nU6lm =    1     *  rho0max   *     1     *    1;
		std::vector<double> U6lm(nU6lm, 0.0);
		SU3Lib.Get6lm(CSU3CGMaster::U6LM, ir1, ir2, ir0, ir34, ir12, *ir2_34, U6lm);
////////////////////////////////////////////
//		int na = rho12max;
//					|
//					V
//		int na 	  = 1;
//
//		int nb = rho12_3max*na;
//                  |        |
//                  V        V
//		int nb = rho0max   * 1;
//
//		int nc = rho23max*nb;
//                  |
//                  V
//		int nc =    1*nb;
////////////////////////////////////////////
		size_t nCoeffs = NewLabel.CoeffsSize(TENSOR::PPNN);
		assert(nCoeffs == 2*k0max); // since rho0pmax = 1 ==> number of tensor components = 2*k0max*rho0pmax = 2*k0max
		std::vector<double> dSU3_234(nCoeffs, 0.0); 

//	iterate over all k0 components of the new tensor [Ux{Vx{XZ}IR1}IR2}IR0
		for (k0 = 0; k0 < k0max; ++k0)
		{
//			index_k0_rho0p = k0*rho0pmax*2 + rho0p*2 + Type
//			                     |
//                               |
//                               V
//                               1
//
			index_k0_rho0p = k0*2;
			for (irho0 = 0, index_k0_rho0 = k0*rho0max*2; irho0 < rho0max; ++irho0, index_k0_rho0 += 2)
			{
//				index =   rho12 + rho12_3 * na + rho23    * nb + rho1_23 * nc; 
//                 	     	|       |       |      |         |     |		|	
//                 	     	|       |       |      |         |     |		|	
//                 	     	|       |       V  irho2_34 = 0  V  irho0p = 0  V
//                 	     	|       |       1      |         1             rho0max
//                 	     	|       |        irho2_34*nb == 0	irho0p*nc == 0	
//			        	    |       |              |               | 
//			            	V       V              V               V
//				index =     0   + irho0   * 1 +    0*nb +          0*nc  ==> index = irho0
				dSU3_234[index_k0_rho0p + PP] += Coeffs[index_k0_rho0 + PP]*U6lm[irho0];
				dSU3_234[index_k0_rho0p + NN] += Coeffs[index_k0_rho0 + NN]*U6lm[irho0];
			}
		}
//	Bugfix: 10/16/2010: 
//	Array of NewCoeffs must be initialized 
//	inside of loop for S234=1/2 as well as S234=3/2
//		std::vector<double> NewCoeffs(dSU3_234);
		for (S234 = 1; S234 <= 3; S234 += 2) // 1/2 x {1/2 x 1/2} ---> S234 = {1/2, 3/2}
		{
			std::vector<double> NewCoeffs(dSU3_234);
			dSU2_S234 = wigner6j(1, 1, S12, S34, S0, S234);
			if (Negligible8(dSU2_S234)) {
				continue;
			}
			dSU2_S234 *= dPhase*sqrt((double)(S234 + 1));
			for (i = 0; i < nCoeffs; ++i)
			{
				NewCoeffs[i] *= dSU2_S234;
			}
			// Store [U x {V x {X x Z}IR1}IR2]IR0 tensor labels and strengths
			// of its components
			NewLabel.IR2.S2 = S234;
			TensorsU__V_XZ.push_back(make_pair(NewLabel, NewCoeffs));
		}
	}
}

//	(A.53)
void SU3InteractionRecoupler::UV_XZtoU__V_XZ(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
												std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU__V_XZ)
{
	UV_XZtoU__V_XZ(	SU3xSU2::LABELS(structure[U]), SU3xSU2::LABELS(structure[V]), SU3xSU2::LABELS(structure[X]), SU3xSU2::LABELS(structure[Z]), 
						TensorUV_XZ, Coeffs, TensorsU__V_XZ);

	// UV0XZ0 ----> UVXZ00	
	structure[V+1] = structure[X];
	structure[X] = structure[Z];
	structure[Z] = 0;
}

//	(A.54)
void SU3InteractionRecoupler::UV_XZtoUV_X__Z(	const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4, 
												const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
												std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z)
{
	static double dExchangePhase = +1.0;

	CSU3CGMaster SU3Lib;

	int rho0max = TensorUV_XZ.IR0.rho;
	SU2::LABEL S12 = TensorUV_XZ.IR1.S2;
	SU2::LABEL S34 = TensorUV_XZ.IR2.S2;
	SU2::LABEL  S0 = TensorUV_XZ.IR0.S2;
	SU2::LABEL S123;

	const SU3xSU2::LABELS& ir12 = TensorUV_XZ.IR1;
	const SU3xSU2::LABELS& ir34 = TensorUV_XZ.IR2;
	const SU3xSU2::LABELS& ir0 	= TensorUV_XZ.IR0;

	TENSOR_LABELS NewLabel;
	NewLabel.IR1 = ir12;
	NewLabel.IR2.rho = 1; 	//	rho12_3max ... is equal to 1 due to the one-by-two coupling of creation/annihilation operators {UV}xX, 
							//	since there are no multiplicities rho > 1 when coupling (lm mu) x (n 0) and (lm mu) x (0 n) irreps
	NewLabel.IR0.rho = 1; 	// 	rho0pmax ...   is equal to 1 due to the one-by-one-by-two coupling of creation/annihilation operators [Ux{Vx{XZ}}]
							//	since there are no multiplicities rho > 1 when coupling (lm mu) x (n 0) and (lm mu) x (0 n) irreps
	NewLabel.IR0.lm = ir0.lm;
	NewLabel.IR0.mu = ir0.mu;
	NewLabel.IR0.S2 = ir0.S2;

//	size_t k0, k0max = SU3::kmax(TensorUV_XZ.IR0, S0/2);
	size_t k0, k0max = 1;
	size_t nU6lm;

	double dPhase = dExchangePhase*MINUSto(1 + (S12 + S0)/2)*sqrt((double)(S34 + 1));
	double dSU2_S123;
	size_t irho0; 
	size_t i, index_k0_rho0, index_k0_rho0p;

	std::vector<SU3::LABELS> Irreps12_3;

	SU3::Couple(ir12, ir3, Irreps12_3);

	std::vector<SU3::LABELS>::const_iterator ir12_3 = Irreps12_3.begin();
	std::vector<SU3::LABELS>::const_iterator LAST_12_3_IRREP = Irreps12_3.end();
	for (; ir12_3 != LAST_12_3_IRREP; ++ir12_3)
	{
//	U[ ir1 ir2  ir0 ir3;  ir12 rho12   rho12_3   ir23  rho23  rho1_23 ]
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |     |  rho12_3  rho0p       |   rho34    |                   
//      |	|	 |	 |     |     |       |         |     |      | 
//      |	|	 |	 |     |     |       |         |     |      | 
//      V   V    V   V     V     V       V         V     V      V 
//	U[ ir12 ir3  ir0 ir4; ir12_3 1       1       ir34    1     rho0  ]

//	 	 rho12_3max  = SU3::mult(ir12,    ir3, ir12_3);
//			   |				  |		   |     |			
//			   |				  |		   |     |
//			   V				  V		   V	 V
		int rho0pmax = SU3::mult(*ir12_3, ir4,  ir0);
		if (rho0pmax == 0) {
			continue;
		}
		assert(rho0pmax == 1); 		//	if rho0pmax > 1 ==> SU3::mult does not work properly!!!
		assert(ir12_3->rho == 1);	//	and similarly, rho2_34max == 1

		NewLabel.IR2.lm = ir12_3->lm;
		NewLabel.IR2.mu = ir12_3->mu;

//		nU6lm = rho12max * rho12_3max * rho23max   * rho1_23max;
//				   |             |          |          |
//				   |             |          |          |
//			    rho12_3      rho0pmax     rho34max   rho0max 
//				   |             |          |          |
//				   |             |          |          |
//				   V             V          V          V
//		nU6lm =    1     *       1   *      1     *  rho0max;
		nU6lm =    rho0max;
		std::vector<double> U6lm(nU6lm, 0.0);
		SU3Lib.Get6lm(CSU3CGMaster::U6LM, ir12, ir3, ir0, ir4, *ir12_3, ir34, U6lm);
////////////////////////////////////////////
//		int na = rho12max;
//					|
//					V
//		int na 	  = 1;
//
//		int nb = rho12_3max*na;
//                  |        |
//                  V        V
//		int nb =    1    *   1;
//
//		int nc = rho23max*nb;
//                  |
//                  V
//		int nc =    1*1;
////////////////////////////////////////////
		size_t nCoeffs = NewLabel.CoeffsSize(TENSOR::PPNN);
		assert(nCoeffs == 2*k0max); // since rho0pmax = 1 ==> number of tensor components = 2*k0max*rho0pmax = 2*k0max
		std::vector<double> dSU3_234(nCoeffs, 0.0); 

//	iterate over all k0 components of the new tensor [Ux{Vx{XZ}IR1}IR2}IR0
		for (k0 = 0; k0 < k0max; ++k0)
		{
//			index_k0_rho0p = k0*rho0pmax*2 + rho0p*2 + Type
//			                     |
//                               |
//                               V
//                               1
//
			index_k0_rho0p = k0*2;
			for (irho0 = 0, index_k0_rho0 = k0*rho0max*2; irho0 < rho0max; ++irho0, index_k0_rho0 += 2)
			{
//				index =   rho12 + rho12_3 * na + rho23    * nb + rho1_23 * nc; 
//                 	     	|       |       |      |         |     |		|	
//                 	     	|       |       |      |         |     |		|	
//                 	    rho12_3  rho0p      V    rho34       V    rho0      V 
//                 	     	|       |       1      |         1     |        1
//                 	     	|       |              |               |
//			        	    |       |              |               | 
//			            	V       V              V               V
//				index =     0   +   0   * na +     0*nb +        irho0*1; 
				dSU3_234[index_k0_rho0p + PP] += Coeffs[index_k0_rho0 + PP]*U6lm[irho0];
				dSU3_234[index_k0_rho0p + NN] += Coeffs[index_k0_rho0 + NN]*U6lm[irho0];
			}
		}

//	Bugfix: 10/16/2010: 
//	Array of NewCoeffs must be initialized 
//	inside of loop for S234=1/2 as well as S234=3/2
//		std::vector<double> NewCoeffs(dSU3_234);
		for (S123 = 1; S123 <= 3; S123 += 2) // {1/2 x 1/2} x 1/2 ---> S123 = {1/2, 3/2}
		{
			std::vector<double> NewCoeffs(dSU3_234);
			dSU2_S123 = wigner6j(S12, 1, S123, 1, S0, S34);
			if (Negligible8(dSU2_S123)) {
				continue;
			}
			dSU2_S123 *= dPhase*sqrt((double)(S123 + 1));
			for (i = 0; i < nCoeffs; ++i)
			{
				NewCoeffs[i] *= dSU2_S123;
			}
			NewLabel.IR2.S2 = S123;
			// Store [{{U x V}IR1 x X}IR2 x Z]IR0 tensor labels and strengths of its components
			TensorsUV_X__Z.push_back(make_pair(NewLabel, NewCoeffs));
		}
	}
}

void SU3InteractionRecoupler::UV_XZtoUV_X__Z(char* structure, const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z)
{
	SU3xSU2::LABELS ir1(structure[U]), ir2(structure[V]), ir3(structure[X]), ir4(structure[Z]);

	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator cit  = TensorsUV_XZ.begin();
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator LAST = TensorsUV_XZ.end();
	for (; cit != LAST; ++cit)
	{
		UV_XZtoUV_X__Z(ir1, ir2, ir3, ir4,  cit->first, cit->second, TensorsUV_X__Z);
	}

	// UV0XZ0 ----> UV0X0Z	
	std::swap(structure[Z], structure[Z+1]);
}

//	(A.54)
void SU3InteractionRecoupler::UV_XZtoUV_X__Z(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
												std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z)
{
	UV_XZtoUV_X__Z(	SU3xSU2::LABELS(structure[U]), SU3xSU2::LABELS(structure[V]), SU3xSU2::LABELS(structure[X]), SU3xSU2::LABELS(structure[Z]), 
					TensorUV_XZ, Coeffs, TensorsUV_X__Z);
	// UV0XZ0 ----> UV0X0Z	
	std::swap(structure[Z], structure[Z+1]);
}


//	(A.55)
void SU3InteractionRecoupler::UV_XZtoU_XZ__V(	const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4, 
												const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
												std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU_XZ__V)
{
	static double dExchangePhase = +1.0;

	CSU3CGMaster SU3Lib;

	int rho0max = TensorUV_XZ.IR0.rho;
	SU2::LABEL S12 = TensorUV_XZ.IR1.S2;
	SU2::LABEL S34 = TensorUV_XZ.IR2.S2;
	SU2::LABEL  S0 = TensorUV_XZ.IR0.S2;
	SU2::LABEL S134;

	const SU3xSU2::LABELS& ir12 = TensorUV_XZ.IR1;
	const SU3xSU2::LABELS& ir34 = TensorUV_XZ.IR2;
	const SU3xSU2::LABELS& ir0 	= TensorUV_XZ.IR0;

	TENSOR_LABELS NewLabel;
	NewLabel.IR1 = ir34;
	NewLabel.IR2.rho = 1; 	//	rho1_34max ... is equal to 1 due to the one-by-two coupling of creation/annihilation operators {UV}xX, 
							//	since there are no multiplicities rho > 1 when coupling (lm mu) x (n 0) and (lm mu) x (0 n) irreps
	NewLabel.IR0.rho = 1; 	// 	rho0pmax ...   is equal to 1 due to the one-by-one-by-two coupling of creation/annihilation operators [Ux{Vx{XZ}}]
							//	since there are no multiplicities rho > 1 when coupling (lm mu) x (n 0) and (lm mu) x (0 n) irreps
	NewLabel.IR0.lm = ir0.lm;
	NewLabel.IR0.mu = ir0.mu;
	NewLabel.IR0.S2 = ir0.S2;


//	size_t k0, k0max = SU3::kmax(TensorUV_XZ.IR0, S0/2);
	size_t k0, k0max = 1;
	size_t nZ6lm;

	double dPhase = dExchangePhase*MINUSto((S34 + S12)/2)*sqrt((double)(S12 + 1));
	double dSU2_S134;
	size_t irho0; 
	size_t i, index_k0_rho0, index_k0_rho0p;

	size_t nCoeffs = NewLabel.CoeffsSize(TENSOR::PPNN);
	assert(nCoeffs == 2*k0max); // since rho0pmax = 1 ==> number of tensor components = 2*k0max*rho0pmax = 2*k0max
	std::vector<double> dSU3_134(nCoeffs);

	std::vector<SU3::LABELS> Irreps1_34;

	SU3::Couple(ir1, ir34, Irreps1_34);

	std::vector<SU3::LABELS>::const_iterator ir1_34 = Irreps1_34.begin();
	std::vector<SU3::LABELS>::const_iterator LAST_1_34_IRREP = Irreps1_34.end();
	for (; ir1_34 != LAST_1_34_IRREP; ++ir1_34)
	{
//	Z[ ir2 ir1  ir0 ir3;  ir12 rho12   rho12_3   ir13  rho13  rho13_2 ] // see Jutta's dissertation page 127, equation C.33
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |     |  rho12    rho0        |  rho1_34 rho0p                      
//      |	|	 |	 |     |     |       |         |     |      | 
//      |	|	 |	 |     |     |       |         |     |      | 
//      V   V    V   V     V     V       V         V     V      V 
//	Z[ ir2 ir1  ir0 ir34; ir12   1     rho0      ir1_34  1      1 ]

//	 	 rho13_2max  = SU3::mult(ir13,    ir2, ir0);
//			   |				  |		   |     |			
//			   |				  |		   |     |
//			   V				  V		   V	 V
		int rho0pmax = SU3::mult(*ir1_34, ir2,  ir0);
		if (rho0pmax == 0) {
			continue;
		}
		assert(rho0pmax == 1); 		//	if rho0pmax > 1 ==> SU3::mult does not work properly!!!
		assert(ir1_34->rho == 1);	//	and similarly, rho2_34max == 1

		NewLabel.IR2.lm = ir1_34->lm;
		NewLabel.IR2.mu = ir1_34->mu;

//		nZ6lm = rho12max * rho12_3max * rho13max   * rho13_2max;
//				   |             |          |          |
//				   |             |          |          |
//			    rho12          rho0     rho1_34max   rho0pmax 
//				   |             |          |          |
//				   |             |          |          |
//				   V             V          V          V
//		nU6lm =    1     *    rho0max   *   1     *    1;
		nZ6lm =    rho0max;
		std::vector<double> Z6lm(nZ6lm, 0);
		SU3Lib.Get6lm(CSU3CGMaster::Z6LM, ir2, ir1, ir0, ir34, ir12, *ir1_34, Z6lm);
////////////////////////////////////////////
//		int na = rho12max;
//					|
//			   	  rho12	
//					|
//					V
//		int na 	  = 1;
//
//		int nb = rho12_3max*na;
//                  |        |
//                rho0max    |
//                  |        |
//                  V        V
//		int nb =  rho0max *  1;
//
//		int nc = rho13max*nb;
//                  |
//              rho1_34max  
//                  |
//                  V
//		int nc =    1*rho0max;
////////////////////////////////////////////

		dSU3_134.assign(nCoeffs, 0.0);
//	iterate over all k0 components of the new tensor
		for (k0 = 0; k0 < k0max; ++k0)
		{
//			index_k0_rho0p = k0*rho0pmax*2 + rho0p*2 + Type
//			                     |
//                               |
//                               V
//                               1
//
			index_k0_rho0p = k0*2;
			for (irho0 = 0, index_k0_rho0 = k0*rho0max*2; irho0 < rho0max; ++irho0, index_k0_rho0 += 2)
			{
//				index =   rho12 + rho12_3 * na + rho13  * nb + rho13_2 * nc; 
//                 	     	|       |       |      |      |     | 		 |	
//                 	     	|       |       |      |      |     |		 |	
//                 	      rho12   rho0      V    rho1_34  V    rho0p     V 
//                 	     	|       |       1      |      1     |        1
//                 	     	|       |              |            |
//			        	    |       |              |            | 
//			            	V       V              V            V
//				index =     0   +  irho0 * na +    0   *  nb +  0    *   1; 
				dSU3_134[index_k0_rho0p + PP] += Coeffs[index_k0_rho0 + PP]*Z6lm[irho0];
				dSU3_134[index_k0_rho0p + NN] += Coeffs[index_k0_rho0 + NN]*Z6lm[irho0];
			}
		}
//	Bugfix: 10/16/2010: 
//	Array of NewCoeffs must be initialized 
//	inside of loop for S234=1/2 as well as S234=3/2
//		std::vector<double> NewCoeffs(dSU3_134);
		for (S134 = 1; S134 <= 3; S134 += 2) //   1/2 x {1/2 x 1/2} ---> S134 = {1/2, 3/2}
		{
			std::vector<double> NewCoeffs(dSU3_134);
			dSU2_S134 = wigner6j(1, 1, S12, S34, S0, S134);
			if (Negligible8(dSU2_S134)) {
				continue;
			}
			dSU2_S134 *= dPhase*MINUSto((1 + S134)/2)*sqrt((double)(S134 + 1));
			for (i = 0; i < nCoeffs; ++i)
			{
				NewCoeffs[i] *= dSU2_S134;
			}
			NewLabel.IR2.S2 = S134;
			// Store [{{U x V}IR1 x X}IR2 x Z]IR0 tensor labels and strengths of its components
			TensorsU_XZ__V.push_back(make_pair(NewLabel, NewCoeffs));
		}
	}
}
//	(A.55)
void SU3InteractionRecoupler::UV_XZtoU_XZ__V(char* structure, const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU_XZ__V)
{
	SU3xSU2::LABELS ir1(structure[U]);
	SU3xSU2::LABELS ir2(structure[V]);
	SU3xSU2::LABELS ir3(structure[X]);
	SU3xSU2::LABELS ir4(structure[Z]);

	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator cit  = TensorsUV_XZ.begin();
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator LAST = TensorsUV_XZ.end();
	for (; cit != LAST; ++cit)
	{
		UV_XZtoU_XZ__V(ir1, ir2, ir3, ir4, cit->first, cit->second, TensorsU_XZ__V);// (A.54)   {b c} x {a a} --> {b {a a}}c
	}

	// UV0XZ0 ----> UXZ00V
	structure[Z+1] = structure[V];
	structure[V] = structure[X];
	structure[V+1] = structure[Z];
	structure[X] = 0;
	structure[Z] = 0;
}
//	(A.55)
void SU3InteractionRecoupler::UV_XZtoU_XZ__V(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
												std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU_XZ__V)
{
	UV_XZtoU_XZ__V(	SU3xSU2::LABELS(structure[U]), SU3xSU2::LABELS(structure[V]), SU3xSU2::LABELS(structure[X]), SU3xSU2::LABELS(structure[Z]), 
					TensorUV_XZ, Coeffs, TensorsU_XZ__V);

	// UV0XZ0 ----> UXZ00V
	structure[Z+1] = structure[V];
	structure[V] = structure[X];
	structure[V+1] = structure[Z];
	structure[X] = 0;
	structure[Z] = 0;
}

//	(A.56)
void SU3InteractionRecoupler::UV_XZtoUX_VZ(const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4,
											const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& TensorCoeffs,
											std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUX_VZ)
{
	static double dExchangePhase = -1.0;

	const SU3xSU2::LABELS& ir12 = TensorUV_XZ.IR1;
	const SU3xSU2::LABELS& ir34 = TensorUV_XZ.IR2;
	const SU3xSU2::LABELS& ir0 	= TensorUV_XZ.IR0;
	const SU2::LABEL& S12 = ir12.S2;
	const SU2::LABEL& S34 = ir34.S2;
	const SU2::LABEL& S0  = ir0.S2;
	SU2::LABEL S13, S24;

	double dPhase = dExchangePhase*sqrt((double)((S12 + 1)*(S34 + 1)));

	TENSOR_LABELS NewLabel;
	NewLabel.IR0 = TensorUV_XZ.IR0;

	int rho1234max = abs(ir0.rho); // rho0max .... ir0.rho could be equal to -1 in case it was constructed with SU3xSU2::LABELS(lm, mu, S2)
	int rho1324max; // rho0pmax

	double dsu2;
	size_t i, irho0p, irho0;
//	size_t k0, k0max = SU3::kmax(TensorUV_XZ.IR0, TensorUV_XZ.IR0.S2/2);
	size_t k0, k0max = 1;
	size_t index, index_k0_rho0p, index_k0_rho0;

	std::vector<SU3::LABELS>::const_iterator ir13, LAST_13_IRREP;
	std::vector<SU3::LABELS>::const_iterator ir24, LAST_24_IRREP;

	std::vector<SU3::LABELS> Irreps13;	
	std::vector<SU3::LABELS> Irreps24;	

	SU3::Couple(ir1, ir3, Irreps13);	// ir1 x ir3 -> ir13
	SU3::Couple(ir2, ir4, Irreps24); 	// ir2 x ir4 -> ir24

	LAST_13_IRREP = Irreps13.end();
	LAST_24_IRREP = Irreps24.end();
	for (ir13 = Irreps13.begin(); ir13 != LAST_13_IRREP; ++ir13)
	{
		NewLabel.IR1.rho = ir13->rho;
		NewLabel.IR1.lm  = ir13->lm;
		NewLabel.IR1.mu  = ir13->mu;
		for (ir24 = Irreps24.begin(); ir24 != LAST_24_IRREP; ++ir24)
		{
			rho1324max = SU3::mult(*ir13, *ir24, ir0); 	// 	Does ir13 x ir24 yields ir0 ?
			if (rho1324max == 0) 						// 	if NOT => continue;
			{
				continue;								
			}
// All tensors must have the same SU(3) label (lm0 mu0) = ir0. However,
// multiplicity is given by SU3::multu(ir13, ir24, ir0) = rho1324max ==
// rho0'_max
			NewLabel.IR0.rho = rho1324max;	

			NewLabel.IR2.rho = ir24->rho;
			NewLabel.IR2.lm  = ir24->lm;
			NewLabel.IR2.mu  = ir24->mu;
	
			size_t nCoeffs = NewLabel.CoeffsSize(TENSOR::PPNN); // nCoefss depends on rho0pmax and hence needs to be evaluated for each rho0pmax ir0
			std::vector<double> NewCoeffs_SU3Part(nCoeffs, 0.0); 
			std::vector<double> NewCoeffs(nCoeffs); 
			std::vector<double> su39lm(rho1234max*rho1324max, 0.0);
			u9lm_look_up_table_.Get9lm(ir1, ir2, ir12, ir3, ir4, ir34, *ir13, *ir24, ir0, &su39lm[0]);

			NewCoeffs_SU3Part.assign(nCoeffs, 0.0); // set all coefficients to 0.0

			for (k0 = 0; k0 < k0max; ++k0)
			{
				for (irho0p = 0, index_k0_rho0p = k0*rho1324max*2; irho0p < rho1324max; ++irho0p, index_k0_rho0p += 2)
//				for (irho0p = 0; irho0p < rho1423max; ++irho0p)
				{ 
//					index_k0_rho0p = k0*rho1423max*2 + irho0p*2; 	
					for (irho0 = 0, index_k0_rho0 = k0*rho1234max*2, index = irho0p*rho1234max; irho0 < rho1234max; ++irho0, index_k0_rho0 += 2, ++index)
//					for (irho0 = 0; irho0 < rho1234max; ++irho0)
					{
//						index_k0_rho0  = k0*rho1234max*2 + irho0*2;	
//						index = irho0p*rho1234max + irho0;			
//	index = rho12 + rho34*n1 + rho1234*n2 + rho14*n3 + rho23*n4 + rho1324*n5;
//	Since rho12 = rho34 = rho13 = rho24 = 0 and rho1234=rho0, rho1324=rho0',n2=1, n5=rho1234max, we have index = irho0 + irho0p*rho1234max
						NewCoeffs_SU3Part[index_k0_rho0p + PP] += TensorCoeffs[index_k0_rho0 + PP] * su39lm[index];
						NewCoeffs_SU3Part[index_k0_rho0p + NN] += TensorCoeffs[index_k0_rho0 + NN] * su39lm[index];
					}
				}
			}
			for (S13 = 0; S13 <= 2; S13 += 2)
			{
				NewLabel.IR1.S2 = S13;
				for (S24 = 0; S24 <= 2; S24 += 2)
				{
					if (SU2::mult(S13, S24, S0) == 0)
					{
						continue;
					}
					NewLabel.IR2.S2 = S24;
					dsu2 = dPhase * sqrt((double)((S13 + 1)*(S24 + 1)))*wigner9j(1, 1, S12, 1, 1, S34, S13, S24, S0);
//std::cout << "wigner9j(1, 1," <<  (int)S12 << ", 1, 1, " << (int)S34 <<", " << (int)S13 << ", " << (int)S24 << ", " << (int)S0 << ") = " << wigner9j(1, 1, S12, 1, 1, S34, S13, S24, S0) << std::endl;
					if (Negligible8(dsu2)) {
						continue;
					}
					for (i = 0; i < nCoeffs; ++i) 
					{
						NewCoeffs[i] = dsu2*NewCoeffs_SU3Part[i]; // multiply all the resulting coefficients by dsu2
					}
					TensorsUX_VZ.push_back(make_pair(NewLabel, NewCoeffs));
				}
			}
		}
	}
}

//	(A.56)
void SU3InteractionRecoupler::UV_XZtoUX_VZ(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
											std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUX_VZ)
{
	UV_XZtoUX_VZ(SU3xSU2::LABELS(structure[U]), SU3xSU2::LABELS(structure[V]), SU3xSU2::LABELS(structure[X]), SU3xSU2::LABELS(structure[Z]), TensorUV_XZ, Coeffs, TensorsUX_VZ);
	// UV0XZ0 ----> UX0VZ0
	std::swap(structure[V], structure[X]);
}

//	(A.50)
void SU3InteractionRecoupler::UV_XZtoXZ_UV(const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsXZ_UV)
{
	static double dExchangePhase = +1;
	CSU3CGMaster SU3Lib;

	SU3xSU2::LABELS ir1 = TensorUV_XZ.IR1;
	SU3xSU2::LABELS ir2 = TensorUV_XZ.IR2;
	SU3xSU2::LABELS ir0 = TensorUV_XZ.IR0;

	double dPhase = dExchangePhase*MINUSto((ir1.S2 + ir2.S2 - ir0.S2)/2); // (-)^{S_1+S_2-S_0}
	size_t irho0p, irho0, index_k0_rho0p, index_k0_rho0, index;
//	size_t k0, k0max = SU3::kmax(ir0, ir0.S2/2);
	size_t k0, k0max = 1;

	int rho0pmax, rho0max;
	rho0pmax = rho0max = ir0.rho;

	std::vector<double> NewCoeffs(Coeffs.size(), 0.0);
//	Z[ ir2 ir1  ir0 ir3;  ir12 rho12   rho12_3   ir13  rho13  rho13_2 ] // see Jutta's dissertation page 127, equation C.33
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |	   |     |       |         |     | 	    |
//      |   |    |   |     |   rho1    rho0        |  rho2    rho0p                      
//      |	|	 |	 |     |     |       |         |     |      | 
//      |	|	 |	 |     |     |       |         |     |      | 
//      V   V    V   V     V     V       V         V     V      V 
//	Z[ ir1 (0 0) ir0 ir2; ir1    1     rho0      ir2     1     rho0p ] == Phi_{rho rho'}

//	nZ6lm = rho12max * rho12_3max * rho13max   * rho13_2max;
//			   |             |          |          |
//			   |             |          |          |
//		     rho1          rho0       rho2    rho0pmax 
//			   |             |          |          |
//			   |             |          |          |
//			   V             V          V          V
//	nZ6lm =    1     *    rho0max   *   1   *  rho0pmax;
	std::vector<double> Phi(rho0max*rho0pmax);
	SU3Lib.Get6lm(CSU3CGMaster::Z6LM, ir1, SU3::LABELS(0, 0), ir0, ir2, ir1, ir2, Phi);

////////////////////////////////////////////
//		int na = rho12max;
//					|
//			   	  rho1max	
//					|
//					V
//		int na 	  = 1;
//
//		int nb = rho12_3max*na;
//                  |        |
//                rho0max    |
//                  |        |
//                  V        V
//		int nb =  rho0max *  1;
//
//		int nc = rho13max*nb;
//                  |
//               rho2max  
//                  |
//                  V
//		int nc =    1*rho0max;
////////////////////////////////////////////

	for (k0 = 0; k0 < k0max; ++k0)
	{
	    index_k0_rho0p = k0*rho0pmax*2;
		for (irho0p = 0; irho0p < rho0pmax; ++irho0p, index_k0_rho0p += 2)
		{	
//			index_k0_rho0p = k0 * rho0pmax * 2 + irho0p * 2 + PP/NN;
			index = irho0p*rho0max;			
			index_k0_rho0 = k0 * rho0max * 2;
			for (irho0 = 0; irho0 < rho0max; ++irho0, index_k0_rho0 += 2, ++index)
			{
//				index_k0_rho0 = k0*rho0max*2 + rho0*2 + PP/NN;				
//				index =   rho12 + rho12_3 * na + rho13  * nb + rho13_2 * nc; 
//                 	     	|       |       |      |      |     | 		 |	
//                 	     	|       |       |      |      |     |		 |	
//                 	      rho1     rho0     V    rho2     V    rho0p     V 
//                 	     	|       |       1      |  rho0max   |      rho0max  
//                 	     	|       |              |            |
//			        	    |       |              |            | 
//			            	V       V              V            V
//				index =     0   +  irho0 * na +    0   *  nb +  rho0p*rho0max; // index = irho0 + irho0p*rho0max;
				NewCoeffs[index_k0_rho0p + PP] += Coeffs[index_k0_rho0 + PP]*Phi[index];
				NewCoeffs[index_k0_rho0p + NN] += Coeffs[index_k0_rho0 + NN]*Phi[index];
			}
			NewCoeffs[index_k0_rho0p + PP] *= dPhase;
			NewCoeffs[index_k0_rho0p + NN] *= dPhase;
		}
	}
	TensorsXZ_UV.push_back(std::make_pair(TENSOR_LABELS(ir2, ir1, ir0), NewCoeffs));
}

//	[{{aa}a}b] --> [b{{aa}a}]
//	structure = {+a +a  0 -a  0 +b}
//			[]=   0  1  2  3  4  5
//			goes to
//	structure = {+b +a +a  0 -a  0}
//			[]=   0  1  2  3  4  5
void SU3InteractionRecoupler::UV_X__ZtoZ__UV_X(char* structure, const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsZ__UV_X)
{
	double dExchangePhase = -1;
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator cit  = TensorsUV_X__Z.begin();
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator LAST = TensorsUV_X__Z.end();

	for (; cit != LAST; ++cit)
	{
		const TENSOR_LABELS& TensorUV_X__Z = cit->first;
		std::vector<double> NewCoeffs(cit->second);

		SU3xSU2::LABELS ir1 = TensorUV_X__Z.IR2;			// [{{aa}^IR1 x a}^IR2 x b]^IR0
		SU3xSU2::LABELS ir2 = SU3xSU2::LABELS(structure[5]); // b
		SU3xSU2::LABELS ir0 = TensorUV_X__Z.IR0;

		double dPhase = dExchangePhase*MINUSto(ir1.lm + ir1.mu + ir2.lm + ir2.mu - ir0.lm - ir0.mu + (ir1.S2 + ir2.S2 - ir0.S2)/2);
		for (size_t i = 0; i < NewCoeffs.size(); ++i)
		{
			NewCoeffs[i] *= dPhase;
		}
		TensorsZ__UV_X.push_back(std::make_pair(TensorUV_X__Z, NewCoeffs));
	}

	char b = structure[5];
	for (size_t i = 5; i > 0; --i)
	{
		structure[i] = structure[i-1];
	}
	structure[0] = b;
}



//	{UV}{XZ} -> {XZ}{UV} [using Phi coefficients]
void SU3InteractionRecoupler::UV_XZtoXZ_UV(char* structure, const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ,  std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsXZ_UV)
{
	//	test whether one can perform this operation
	//	for example {a+ b+}{c a} IS NOT EQUAL TO {c a}{a+ b+} since {a+_{i},a_{j}}=\delta_{ij} - a_{j}a+_{i} !!!
	assert(abs(structure[U]) != abs(structure[X]) && abs(structure[U]) != abs(structure[Z]));
	assert(abs(structure[V]) != abs(structure[X]) && abs(structure[V]) != abs(structure[Z]));

	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator cit  = TensorsUV_XZ.begin();
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator LAST = TensorsUV_XZ.end();
	for (; cit != LAST; ++cit)
	{
		UV_XZtoXZ_UV(cit->first, cit->second, TensorsXZ_UV);
	}

	std::swap(structure[U], structure[X]);
	std::swap(structure[V], structure[Z]);
}

//	(A.50)
void SU3InteractionRecoupler::UV_XZtoXZ_UV(char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsXZ_UV)
{
	//	test whether one can perform this operation
	//	for example {a+ b+}{c a} IS NOT EQUAL TO {c a}{a+ b+} since {a+_{i},a_{j}}=\delta_{ij} - a_{j}a+_{i} !!!
	assert(abs(structure[U]) != abs(structure[X]) && abs(structure[U]) != abs(structure[Z]));
	assert(abs(structure[V]) != abs(structure[X]) && abs(structure[V]) != abs(structure[Z]));

	UV_XZtoXZ_UV(TensorUV_XZ, Coeffs, TensorsXZ_UV);

	std::swap(structure[U], structure[X]);
	std::swap(structure[V], structure[Z]);
}



void SU3InteractionRecoupler::RecoupleTwoShellTensor(char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TransformedTensors)
{
	if (abs(structure[U]) == abs(structure[V]))  // {a a} x {? ?}
	{
		if (abs(structure[X]) == abs(structure[Z]))// {a a} x {b b} has two cases (1) a < b (2) a > b
		{								
			int a = abs(structure[U]);
			int b = abs(structure[X]);
			if (a < b)																				//	{a a}{b b} with a < b ==> no recoupling needed
			{
				TransformedTensors.push_back(std::make_pair(TensorUV_XZ, Coeffs));
				return; 														// no recoupling needed
			}
			else if (a > b)																			//	{a a} x {b b} with a > b ==> not implemented
			{
				UV_XZtoXZ_UV(structure, TensorUV_XZ, Coeffs, TransformedTensors);//	(A.50)
				return; 														
			}
		}
		else if (abs(structure[V]) == abs(structure[X])) 											//	{a a} x {a b}
		{																							//	 U V     X Z
			UV_XZtoUV_X__Z(structure, TensorUV_XZ, Coeffs, TransformedTensors);	// (A.54)
			return;
		}
		else if (abs(structure[V]) == abs(structure[Z]))											//	{a a} x {b a}
		{																							//	 U V     X Z
			std::vector<double> CoeffsUV_ZX(Coeffs);
			UV_XZtoUV_ZX(structure, TensorUV_XZ.IR2, CoeffsUV_ZX);

			std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUV_X__Z;
			UV_XZtoUV_X__Z(structure, TensorUV_XZ, CoeffsUV_ZX, TensorsUV_X__Z);
			UV_X__ZtoZ__UV_X(structure, TensorsUV_X__Z, TransformedTensors);
			return; 									
		}
	}
	else // {ab}{??} or {ba}{??}
	{
		if (abs(structure[U]) == abs(structure[X]) && abs(structure[V]) == abs(structure[Z]))			//	{a b} x {a b}
		{ 																								//	 U V     X Z
			UV_XZtoUX_VZ(structure, TensorUV_XZ, Coeffs, TransformedTensors); // (A.56)
			return;
		}
		else if (abs(structure[U]) != abs(structure[X]) && abs(structure[V]) == abs(structure[Z]))		//	{b a} x {a a}
		{																								//	 U V     X Z
			UV_XZtoU__V_XZ(structure, TensorUV_XZ, Coeffs, TransformedTensors); //	(A.53)
			return;
		}
		else if (abs(structure[U]) == abs(structure[X]) && abs(structure[V]) != abs(structure[Z]))		//	{a b} x {a a}
		{																								//	 U V     X Z
			UV_XZtoU_XZ__V(structure, TensorUV_XZ, Coeffs, TransformedTensors); //	(A.55)
			return;
		}
	}	
}

void SU3InteractionRecoupler::RecoupleThreeShellTensor(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
														std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TransformedTensors)
{
	if (abs(structure[U]) == abs(structure[V])) // {a a} x {b c}
	{
		int a = abs(structure[U]);
		int b = abs(structure[X]);
		int c = abs(structure[Z]);

		if (a < b && b < c)
		{
			UV_XZtoUV_X__Z(structure, TensorUV_XZ, Coeffs, TransformedTensors);	// (A.54)
			return;
		}
		else if (b < a && a < c)
		{
			std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsXZ_UV;
			UV_XZtoXZ_UV(structure, TensorUV_XZ, Coeffs, TensorsXZ_UV); // (A.50)
			UV_XZtoU_XZ__V(structure, TensorsXZ_UV, TransformedTensors);// (A.55)
			return;	
		}
		else if (b < c && c < a)
		{
			UV_XZtoXZ_UV(structure, TensorUV_XZ, Coeffs, TransformedTensors);//	(A.50)
			return;
		}
	}
	else if (abs(structure[X]) == abs(structure[Z])) // {b c} x {a a}
	{
		int a = abs(structure[X]);
		int b = abs(structure[U]);
		int c = abs(structure[V]);

		if (b < c && c < a)
		{
			TransformedTensors.push_back(std::make_pair(TensorUV_XZ, Coeffs)); // no recoupling needed
			return;
		}
		else if (b < a && a < c)
		{
			UV_XZtoU_XZ__V(structure, TensorUV_XZ, Coeffs, TransformedTensors); // (A.55)
			return;
		}
		else if (a < b && b < c)
		{
			std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsXZ_UV;
			UV_XZtoXZ_UV(structure, TensorUV_XZ, Coeffs, TensorsXZ_UV);		//	(A.50)
			UV_XZtoUV_X__Z(structure, TensorsXZ_UV, TransformedTensors);	//	(A.54)
			return;
		}
	}
	else if (abs(structure[U]) == abs(structure[X])) // {a b} x {a c} 
	{
		int b = abs(structure[V]);
		int c = abs(structure[Z]);

		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUX_VZ;
		UV_XZtoUX_VZ(structure, TensorUV_XZ, Coeffs, TensorsUX_VZ); // (A.56) ===> {a a} {b c}
		if (b > c)
		{
			UV_XZtoUV_ZX(structure, TensorsUX_VZ); // (A.52)
		}
		UV_XZtoUV_X__Z(structure, TensorsUX_VZ, TransformedTensors);// (A.54)
		return;
	}
	else if (abs(structure[V]) == abs(structure[Z]))	//	{b a} x {c a}
	{
		int b = abs(structure[U]);
		int c = abs(structure[X]);

		UV_XZtoUX_VZ(structure, TensorUV_XZ, Coeffs, TransformedTensors);	//	(A.56) ==> {b c} {a a}
		if (b > c)
		{
			UV_XZtoVU_XZ(structure, TransformedTensors);	//	(A.51)
		}
		return;
	}
	else if (abs(structure[V]) == abs(structure[X])) 	//	{b a} x {a c}
	{													//	 U V     X Z
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUZ_VX;
		UV_XZtoUZ_VX(TENSOR::PPNN, structure, TensorUV_XZ, Coeffs, TensorsUZ_VX); 	//	(A.57) ---> {b c} x {a a}
		UV_XZtoU_XZ__V(structure, TensorsUZ_VX, TransformedTensors);				//	(A.55) ---> {b {a a}}c
		return;
	}
	else if (abs(structure[U]) == abs(structure[Z])) //	{a b} x {c a}
	{												 //	 U V     X Z
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUZ_VX;
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsVX_UZ;
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsXV_UZ = TensorsVX_UZ;

		UV_XZtoUZ_VX(TENSOR::PPNN, structure, TensorUV_XZ, Coeffs, TensorsUZ_VX); 	//	(A.57) ---> {a a} x {b c}
		UV_XZtoXZ_UV(structure, TensorsUZ_VX, TensorsVX_UZ); 						//	(A.50) ---> {b c} x {a a}
		UV_XZtoVU_XZ(structure, TensorsXV_UZ);										//	(A.51) ---> {c b} x {a a}
		UV_XZtoU_XZ__V(structure, TensorsXV_UZ, TransformedTensors);				//	(A.55) ---> {c {a a}} b
		return;
	}
}

void SU3InteractionRecoupler::RecoupleFourShellTensor(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
														std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TransformedTensors)
{
	int a = abs(structure[U]);
	int b = abs(structure[V]);
	int c = abs(structure[X]);
	int d = abs(structure[Z]);

//	cout << "Condition: " << (c < a && a < b && b < d) << endl;

	//	{a b} x {c d}
	if (a < b && b < c && c < d)
	{
		UV_XZtoUV_X__Z(structure, TensorUV_XZ, Coeffs, TransformedTensors); // (A.54) ---> {{a b} x c} x d
		return;
	}
	else if (a < c && c < b && b < d)
	{
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUX_VZ;
		UV_XZtoUX_VZ(structure, TensorUV_XZ, Coeffs, TensorsUX_VZ); // (A.56) ---> {a c} x {b d}
		UV_XZtoUV_X__Z(structure, TensorsUX_VZ, TransformedTensors);// (A.54) ---> {{a c} x b} x d
		return;
	}
	else if (a < c && c < d && d < b)
	{
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUX_VZ;
		UV_XZtoUX_VZ(structure, TensorUV_XZ, Coeffs, TensorsUX_VZ); // (A.56) ---> {a c} x {b d}
		UV_XZtoUV_ZX(structure, TensorsUX_VZ); 						// (A.52) ---> {a c} x {d b}
		UV_XZtoUV_X__Z(structure, TensorsUX_VZ, TransformedTensors);// (A.54) ---> {{a c} x d} x b
		return;
	}
	else if (c < a && a < b && b < d)
	{
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUX_VZ;
		UV_XZtoUX_VZ(structure, TensorUV_XZ, Coeffs, TensorsUX_VZ); // (A.56) ---> {a c} {b d}
		UV_XZtoVU_XZ(structure, TensorsUX_VZ); 						// (A.51) ---> {c a} {b d}
		UV_XZtoUV_X__Z(structure, TensorsUX_VZ, TransformedTensors);// (A.54) ---> {{c a} x b} x d
		return;
	}
	else if (c < a && a < d && d < b)
	{
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsUX_VZ;
		UV_XZtoUX_VZ(structure, TensorUV_XZ, Coeffs, TensorsUX_VZ); // (A.56) ---> {a c} {b d} 
		UV_XZtoVU_XZ(structure, TensorsUX_VZ); 						// (A.51) ---> {c a} {b d}
		UV_XZtoUV_ZX(structure, TensorsUX_VZ); 						// (A.52) ---> {c a} {d b}
		UV_XZtoUV_X__Z(structure, TensorsUX_VZ, TransformedTensors);// (A.54) ---> {{c a} x d} x b
		return;
	}
	else if (c < d && d < a && a < b)
	{
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TensorsXZ_UV;
		UV_XZtoXZ_UV(structure, TensorUV_XZ, Coeffs, TensorsXZ_UV);		//	(A.50)
		UV_XZtoUV_X__Z(structure, TensorsXZ_UV, TransformedTensors);	//	(A.54)
		return;
	}
}

bool SU3InteractionRecoupler::Insert_adad_aa_Tensor(const char* n1n2n3n4, const SU3xSU2::LABELS& IR1, const SU3xSU2::LABELS& IR2, const SU3xSU2::LABELS& IR0, const std::vector<double>& TensorCoeffs) 
{
	const int CREATION = +1; 
	const int ANNIHILATION = -1;
//	assert(SU3::kmax(IR0, IR0.S2/2)); // check whether L0 == S0 is allowed is allowed

	std::set<char>	shells;
	shells.insert(abs(n1n2n3n4[0])); shells.insert(abs(n1n2n3n4[1])); shells.insert(abs(n1n2n3n4[2])); shells.insert(abs(n1n2n3n4[3]));
	size_t nActiveShells = shells.size(); // number of different shells in the operator 

	size_t nCoeffsPPNN = OPERATOR::CoeffsSizePPNN(IR0);
	size_t nCoeffsPN   = OPERATOR::CoeffsSizePN(IR0);
	size_t nCoeffsTotal = nCoeffsPN + nCoeffsPPNN;

	assert(nCoeffsTotal == TensorCoeffs.size());
	
	std::vector<std::pair<TENSOR_LABELS, std::vector<double> > > TransformedTensorsPN, TransformedTensorsPPNN;

	size_t nZerosPN(0), nZerosPPNN(0);

	char structure_initial[6]; // +(n1+1) +(n2+1) 0 -(n3+1) -(n4+1) 0
	structure_initial[U] = CREATION * (n1n2n3n4[0] + 1);
	structure_initial[V] = CREATION * (n1n2n3n4[1] + 1);
	structure_initial[2] = 0;
	structure_initial[X] = ANNIHILATION * (n1n2n3n4[2] + 1);
	structure_initial[Z] = ANNIHILATION * (n1n2n3n4[3] + 1);
	structure_initial[5] = 0;

	assert(SU3::mult(SU3xSU2::LABELS(structure_initial[U]), SU3xSU2::LABELS(structure_initial[V]), IR1));
	assert(SU3::mult(SU3xSU2::LABELS(structure_initial[X]), SU3xSU2::LABELS(structure_initial[Z]), IR2));
	assert(SU3::mult(IR1, IR2, IR0));

	for (size_t index = 0; index < nCoeffsTotal; index += 3)
	{
		nZerosPPNN 	+= (Negligible(TensorCoeffs[index + PP]) == true) + (Negligible(TensorCoeffs[index + NN]) == true);
		nZerosPN 	+= (Negligible(TensorCoeffs[index + PN]) == true);
	}

	if (nZerosPN < nCoeffsPN)
	{
		char structure[6];
		memcpy(structure, structure_initial, sizeof(structure_initial));

		std::vector<double> Coeffs(nCoeffsPN);
		for (size_t i = 0, index = 0; i < nCoeffsPN; ++i, index += 3)
		{
			Coeffs[i] = TensorCoeffs[index + PN];
		}

		//	Type == TENSOR::PN needs to be specified in order to correctly transform from Coeffs[a_{pp}, a_{nn}, a_{pn}, a'_{pp}, a'_{nn}, a'_{pn}, ... ]
		//	into TransformedTensorsPN which stores {a_{pn}, a'_{pn}, a''_{pn} ...} with the
		UV_XZtoUZ_VX(TENSOR::PN, structure, TENSOR_LABELS(IR1, IR2, IR0), Coeffs, TransformedTensorsPN); // {pn}{np} -> {pp}{nn}

		if (abs(structure[U]) > abs(structure[V]))
		{
			// This function just multiplies all coefficients in TransformedTensorsPN by a phase ===> Type (TENSOR::PPNN or TENSOR::PN doesn't need be specified)
			UV_XZtoVU_XZ(structure, TransformedTensorsPN); // [{p_{1} > p_{2}} x {n_{1} < n_{2}}] ----> [{p_{2} < p_{1}} x {n_{1} n_{2}}]
		}
		if (abs(structure[X]) > abs(structure[Z]))
		{
			UV_XZtoUV_ZX(structure, TransformedTensorsPN); // [{p_{1} < p_{2}}{n_{1} > n_{2}] ----> [{p_{1} < p_{2}} x {n_{1} < n_{2}}]
		}
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator cit  = TransformedTensorsPN.begin();
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator LAST = TransformedTensorsPN.end();
		for (; cit != LAST; ++cit)
		{
			Add(TENSOR::PN, nActiveShells, structure, cit->first, cit->second);
		}
	}

	if (nZerosPPNN < nCoeffsPPNN)
	{
		char structure[6];
		memcpy(structure, structure_initial, sizeof(structure_initial));

		std::vector<double> Coeffs(nCoeffsPPNN);
		for (size_t i = 0, index = 0; i < nCoeffsPPNN; index += 3)
		{
			Coeffs[i++] = TensorCoeffs[index + PP];
			Coeffs[i++] = TensorCoeffs[index + NN];
		}

		if (abs(structure[U]) > abs(structure[V]))
		{
			UV_XZtoVU_XZ(structure, IR1, Coeffs);
		}
		if (abs(structure[X]) > abs(structure[Z]))
		{
			UV_XZtoUV_ZX(structure, IR2, Coeffs);
		}
		switch (nActiveShells)
		{
			case 1: TransformedTensorsPPNN.push_back(make_pair(TENSOR_LABELS(IR1, IR2, IR0), Coeffs)); break;
			case 2: RecoupleTwoShellTensor(structure, TENSOR_LABELS(IR1, IR2, IR0), Coeffs, TransformedTensorsPPNN);break;
			case 3: RecoupleThreeShellTensor(structure, TENSOR_LABELS(IR1, IR2, IR0), Coeffs, TransformedTensorsPPNN);break;
			case 4: RecoupleFourShellTensor(structure, TENSOR_LABELS(IR1, IR2, IR0), Coeffs, TransformedTensorsPPNN);break;
		}
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator cit  = TransformedTensorsPPNN.begin();
		std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >::const_iterator LAST = TransformedTensorsPPNN.end();
		for (; cit != LAST; ++cit)
		{
			Add(TENSOR::PPNN, nActiveShells, structure, cit->first, cit->second);
		}
	}
	return ((TransformedTensorsPPNN.size() > 0) || (TransformedTensorsPN.size() > 0));
}

} //namespace recoupler_nonscalar
