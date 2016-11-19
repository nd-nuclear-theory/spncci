/****************************************************************
  shell_comparison.cpp
                                
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <cmath>
#include "cppformat/format.h"

#include "sp3rlib/u3coef.h"
#include "moshinsky/moshinsky_xform.h"
namespace u3shell
{
void
MoshinskyTransformTensor(
  const OperatorLabelsU3ST& operator_labels,
  int etap, int eta,
  const u3shell::TwoBodySubspaceU3ST& bra_subspace, 
  const u3shell::TwoBodySubspaceU3ST& ket_subspace,
  int rho0, 
  std::string normalization,
  double coef,
  u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_expansion
  )
{
      u3::SU3 x0(operator_labels.x0());
      u3::SU3 xp(bra_subspace.omega().SU3());
      u3::SU3 x(ket_subspace.omega().SU3());      
      HalfInt Sp(bra_subspace.S());
      HalfInt S(ket_subspace.S());
      HalfInt Tp(bra_subspace.T());
      HalfInt T(bra_subspace.T());

      int dimb=bra_subspace.size();
      int dimk=ket_subspace.size();
      // set up two-body sector for given operator,(omegap,Sp,Tp),(omega,S,T), rho0
      Eigen::MatrixXd sector(dimk,dimb);
      sector=MoshinskyTransform(x0, etap, eta, bra_subspace, ket_subspace, rho0, normalization);
      // Accumluating the moshinsky transformed coeffcients 
      for (int i=0; i<dimb; ++i)
      {
        const u3shell::TwoBodyStateU3ST bra_state(bra_subspace,i);
        int eta1p=bra_state.N1();
        int eta2p=bra_state.N2();
        u3shell::TwoBodyStateLabelsU3ST bra(eta1p, eta2p, xp, Sp, Tp);
        for (int j=0; j<dimk; ++j)
          {
            const u3shell::TwoBodyStateU3ST ket_state(ket_subspace,j);
            int eta1=ket_state.N1();
            int eta2=ket_state.N2();
            u3shell::TwoBodyStateLabelsU3ST ket(eta1, eta2, x, S, T);

            double rme=coef*sector(i,j);      
            TwoBodyUnitTensorLabelsU3ST tboperator(operator_labels,rho0,bra,ket);
            two_body_expansion[tboperator]+=rme;                
           }
        }
      // remove unit tensors with coefficient zero
      std::vector<TwoBodyUnitTensorLabelsU3ST> delete_list;
      for(auto key_value : two_body_expansion)
      {
        if(fabs(key_value.second)<10e-10)
          delete_list.push_back(key_value.first);
      }
      for(int i=0; i<delete_list.size(); ++i)
        {
          auto key=delete_list[i];
          two_body_expansion.erase(key);
        }
}

void 
  MoshinskyTransformUnitTensor(
    const u3shell::RelativeUnitTensorLabelsU3ST& tensor, 
    u3shell::RelativeCMUnitTensorCache& unit_relative_cm_expansion,
    u3shell::TwoBodySpaceU3ST& space,
    u3shell::TwoBodyUnitTensorCoefficientsU3ST& two_body_expansion,
    std::string normalization,
    double expansion_coef
  )
  {
    //extract operator labels 
    int etap=tensor.bra().eta();
    int eta=tensor.ket().eta();
    HalfInt S=tensor.ket().S();
    HalfInt T=tensor.ket().T();
    HalfInt Sp=tensor.bra().S();
    HalfInt Tp=tensor.bra().T();
    u3::SU3 x0(tensor.x0());
    HalfInt S0(tensor.S0());
    HalfInt T0(tensor.T0());

    u3shell::OperatorLabelsU3ST operator_labels(etap-eta,x0,tensor.S0(), tensor.T0(), tensor.g0());
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sectors are generated subject to the constraints that 
    // omega.N()+N0=omegap.N()
    // omega x x0 -> omegap
    // S x S0 ->Sp
    // T x T0 ->Tp
    u3shell::TwoBodySectorsU3ST sector_labels_list(space,tensor);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // For each two-body sector, moshinsky transform the relative-cm matrix element
    // and accumulate values in container two_body_expansion.
    for (int sector_index=0; sector_index<sector_labels_list.size(); sector_index++)
    {
      // for a given sector, extract bra and ket subspace information
      auto sector_labels=sector_labels_list.GetSector(sector_index);
      const u3shell::TwoBodySubspaceU3ST& bra_subspace(sector_labels.bra_subspace());
      const u3shell::TwoBodySubspaceU3ST& ket_subspace(sector_labels.ket_subspace());
      int rho0=sector_labels.multiplicity_index();

      if ((bra_subspace.S()!=Sp) || (bra_subspace.T()!=Tp) || (ket_subspace.S()!=S) || (ket_subspace.T()!=T))
        continue;

      u3::SU3 xp(bra_subspace.omega().SU3());
      u3::SU3 x(ket_subspace.omega().SU3());
      HalfInt Np(bra_subspace.omega().N());
      HalfInt N(ket_subspace.omega().N());

      int eta_cm=int(Np-etap);

      u3shell::RelativeCMStateLabelsU3ST bra(etap,eta_cm,xp,Sp,Tp);
      u3shell::RelativeCMStateLabelsU3ST ket(eta,eta_cm,x,S,T);
      u3shell::RelativeCMUnitTensorLabelsU3ST braket_u3st(x0,S0,T0,rho0,bra,ket);
      
      if(unit_relative_cm_expansion.count(braket_u3st)==0)
        continue;

      // Adding in center of mass
      expansion_coef=expansion_coef*unit_relative_cm_expansion[braket_u3st];
      if(fabs(expansion_coef)<10e-8)
        continue;
      ///////////////////////////////////////////////////////////////////////////////////////
      
      MoshinskyTransformTensor(operator_labels,etap, eta, bra_subspace, ket_subspace, rho0, 
        normalization, expansion_coef,two_body_expansion);
    }
  }

	void RelativeUnitTensorToTwobodyU3ST(int Nmax,  
	  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
	  u3shell::RelativeCMExpansion& unit_relative_cm_map,
	  TwoBodyExpansionMap& two_body_expansion_vector,
	  std::string normalization
	  )
	{
	  u3shell::TwoBodySpaceU3ST space(Nmax);
	  for(auto tensor : relative_unit_tensors)
	    {
	      u3shell::RelativeCMUnitTensorCache& unit_relative_cm_expansion=unit_relative_cm_map[tensor];
	      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_expansion;
	      // u3shell::MoshinskyTransformUnitTensor(tensor, unit_relative_cm_expansion,space, two_body_expansion,normalization);
	      two_body_expansion_vector[tensor]=two_body_expansion;
	    }
	}


}//namespace