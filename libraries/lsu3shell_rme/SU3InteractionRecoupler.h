#ifndef  	SU3INTERACTIONRECOUPLER_H
#define 	SU3INTERACTIONRECOUPLER_H
#include <UNU3SU3/UNU3SU3Basics.h>
#include <SU3NCSMUtils/clebschGordan.h>
#include <SU3ME/global_definitions.h>
#include <UNU3SU3/CSU3Master.h>
#include <LookUpContainers/CSU39lm.h>
#include <iostream>
#include <set>
#include <cassert>

namespace TENSOR
{
	enum Type {PPNN = 0, PN = 1};
}

namespace OPERATOR
{
/*	
	inline size_t CoeffsSizePN(SU3xSU2::LABELS IR0)   	{return IR0.rho * SU3::kmax(IR0, IR0.S2/2);}
	inline size_t CoeffsSizePPNN(SU3xSU2::LABELS IR0) 	{return 2*IR0.rho * SU3::kmax(IR0, IR0.S2/2);}
	inline size_t CoeffsSize(SU3xSU2::LABELS IR0) 		{return 3*IR0.rho * SU3::kmax(IR0, IR0.S2/2);}
*/
	inline size_t CoeffsSizePN(SU3xSU2::LABELS IR0)   	{return IR0.rho * 1;}
	inline size_t CoeffsSizePPNN(SU3xSU2::LABELS IR0) 	{return 2*IR0.rho * 1;}
	inline size_t CoeffsSize(SU3xSU2::LABELS IR0) 		{return 3*IR0.rho * 1;}
}

//	index_k0_rho0 = k0*rho0max*2 + rho0*2 + PP/NN 	if TENSOR::TYPE == PPNN
//	index_k0_rho0 = k0*rho0max + rho0 				if TENSOR::TYPE == PN
struct TENSOR_LABELS
{
	public:
	SU3xSU2::LABELS IR1, IR2, IR0;
	public:
	TENSOR_LABELS(const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir0):IR1(ir1), IR2(ir2), IR0(ir0) {assert(ir0.rho >= 0);}
	TENSOR_LABELS():IR1(1, 0, 0, 0), IR2(1, 0, 0, 0), IR0(1, 0, 0, 0) {};

	inline size_t CoeffsSize(const TENSOR::Type Type) const { return (Type == TENSOR::PN) ? OPERATOR::CoeffsSizePN(IR0): OPERATOR::CoeffsSizePPNN(IR0); }

	inline bool operator<(const TENSOR_LABELS& R) const
	{
		return 	(IR1 < R.IR1)
				|| (
				(IR1 == R.IR1 && IR2 < R.IR2)
				|| (
				(IR1 == R.IR1 && IR2 == R.IR2 && IR0 < R.IR0)));
	}
};

class SU3InteractionRecoupler
{
	private:
	enum {U = 0, V = 1, X = 3, Z = 4};
	public:
	CSU39lm<double> u9lm_look_up_table_;
//	TENSOR_COMPONENT:
//	TENSOR_LABELS      ------------------------->  {app, ann, app, ann, app, ann} or {apn, apn, apn, apn}
//	index = k0*rho0max*2 + rho0*2 + TYPE, where {TYPE PP == 0, and TYPE NN == 1};
//	index = k0*rho0max + rho0 for PPNN
//
	typedef std::map<TENSOR_LABELS, std::vector<double> >  TENSOR_COMPONENTS;
	std::vector<std::map<CTuple<char, 6>, TENSOR_COMPONENTS> > m_PPNN; 
	std::vector<std::map<CTuple<char, 6>, TENSOR_COMPONENTS> > m_PN; 

	public:
	// Please keep in mind that we can accomodate 30 x 10^6 U6/Z6 symbols. 
	// In case we need more, code will crash with segmentation fault.
	// Possible fix: 
	// (a) increase size
	// (b) compute U9 on the fly
	// (c) replace HashFixed in CSU39lm with LRUhash
	SU3InteractionRecoupler(): u9lm_look_up_table_(30000000,30000000), m_PPNN(4, std::map<CTuple<char, 6>, TENSOR_COMPONENTS>()), m_PN(4, std::map<CTuple<char, 6>, TENSOR_COMPONENTS>()) {};

	bool Insert_adad_aa_Tensor(const char* n1n2n3n4, const SU3xSU2::LABELS& IR1, const SU3xSU2::LABELS& IR2, const SU3xSU2::LABELS& IR0, const std::vector<double>& TensorCoeffs);

/////////////////////////////////////////////////////////////////////////////////////////////
////////               Methods that handle the internal datastructures               ////////               
/////////////////////////////////////////////////////////////////////////////////////////////
	void Add(const TENSOR::Type Type, const char nShells, const char* n1n2n3n4, const TENSOR_LABELS& TensorLabels, const std::vector<double>& TensorCoeffs);
	void RemoveTensorsWithAmplitudesLess(const float eps);
	inline void clear() { m_PN[0].clear();m_PPNN[0].clear(); m_PN[1].clear();m_PPNN[1].clear(); m_PN[2].clear();m_PPNN[2].clear(); m_PN[3].clear();m_PPNN[3].clear();}
	void Show();
	void Show2();
	void Save(const std::string& sOutputFileName, const TENSOR::Type Type);
	void PrintInfo();
	void PrintSameShellTensors() {};
/////////////////////////////////////////////////////////////////////////////////////////////////
////////               R E C O U P L I N G    T R A N S F O R M A T I O N S              ////////               
/////////////////////////////////////////////////////////////////////////////////////////////////

//	(A.57)
	void UV_XZtoUZ_VX(	const TENSOR::Type Type, const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4,
						const TENSOR_LABELS& TensorLabel, const std::vector<double>& TensorCoeffs,
						std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& PNTensorComponents);
	inline void UV_XZtoUZ_VX(	const TENSOR::Type Type, char* structure, 
						const TENSOR_LABELS& TensorLabels, const std::vector<double>& Coeffs, 
						std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorComponents)
	{
		UV_XZtoUZ_VX(Type, SU3xSU2::LABELS(structure[U]), SU3xSU2::LABELS(structure[V]), SU3xSU2::LABELS(structure[X]), SU3xSU2::LABELS(structure[Z]), TensorLabels, Coeffs, TensorComponents);
		std::swap(structure[V], structure[Z]); 
		std::swap(structure[X], structure[Z]);
	}


	void UV_X__ZtoZ__UV_X(char* structure, const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsZ__UV_X);

	void UV_XZtoUV_ZX(char* structure, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ);
	void UV_XZtoUV_ZX(char* structure, const SU3xSU2::LABELS& IR34, std::vector<double>& Coeffs);

	void UV_XZtoVU_XZ(char* structure, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ);
	void UV_XZtoVU_XZ(char* structure, const SU3xSU2::LABELS& IR12, std::vector<double>& Coeffs);

// this function facilitates the two transformations: UV_XZ -> UV_ZX and UV_XZ
// -> VU_XZ, where it is assumed that U, V, X, and Z are all composed by the
// single creation/annihilation operator
	void FlipTensors (const SU3::LABELS& irA, const SU3::LABELS& irB,  const SU3xSU2::LABELS& irAB, std::vector<double>& Coeffs);

//	(A.50)
	void UV_XZtoXZ_UV(const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsXZ_UV);
	void UV_XZtoXZ_UV(char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsXZ_UV);
	void UV_XZtoXZ_UV(char* structure, const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ,  std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsXZ_UV);


//	(A.53)
	void UV_XZtoU__V_XZ(	const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4, 
							const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU__V_XZ);
	void UV_XZtoU__V_XZ(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
							std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU__V_XZ);

//	(A.54)
	void UV_XZtoUV_X__Z(	const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4, 
							const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z);
	void UV_XZtoUV_X__Z(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
												std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z);
	void UV_XZtoUV_X__Z(	char* structure, 
							const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ, 
							std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_X__Z);


//	(A.55)
	void UV_XZtoU_XZ__V(	const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4, 
							const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
							std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU_XZ__V);
	void UV_XZtoU_XZ__V(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
							std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU_XZ__V);
	void UV_XZtoU_XZ__V(char* structure, const std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUV_XZ, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsU_XZ__V);


//	(A.56)
	void UV_XZtoUX_VZ(	const SU3xSU2::LABELS& ir1, const SU3xSU2::LABELS& ir2, const SU3xSU2::LABELS& ir3, const SU3xSU2::LABELS& ir4,
						const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs,
						std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUX_VZ);
	void UV_XZtoUX_VZ(	char* structure, const TENSOR_LABELS& TensorUV_XZ, const std::vector<double>& Coeffs, 
						std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TensorsUX_VZ);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
////////               Recoupling logic              ////////               
/////////////////////////////////////////////////////////////
	void RecoupleTwoShellTensor(	char* structure, const TENSOR_LABELS& TensorXY_ZX, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TransformedTensors);
	void RecoupleThreeShellTensor(	char* structure, const TENSOR_LABELS& TensorXY_ZX, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TransformedTensors);
	void RecoupleFourShellTensor(	char* structure, const TENSOR_LABELS& TensorXY_ZX, const std::vector<double>& Coeffs, std::vector<std::pair<TENSOR_LABELS, std::vector<double> > >& TransformedTensors);
};
#endif 
