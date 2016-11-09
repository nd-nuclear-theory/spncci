/****************************************************************
  relative_cm_xform.cpp
                                
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <fstream>
#include <iostream>
#include "cppformat/format.h"
#include "basis/lsjt_operator.h"

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/import_interaction.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/moshinsky_xform.h"
#include "u3shell/tensor_labels.h"
#include "moshinsky/relative_cm_xform.h"

namespace u3shell
{

  void
  BranchNLST(
    const u3shell::RelativeUnitTensorLabelsU3ST& relative_unit_tensor,
    u3::WCoefCache& w_cache,
    std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors
    )
  {
    int N0, etap,eta,eta_cm;
    HalfInt S0,T0,Sp,Tp,S,T;
    u3::SU3 x0;
    std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=relative_unit_tensor.FlatKey();
    MultiplicityTagged<int>::vector L0_set=u3::BranchingSO3(x0);
    for(auto Lk0 : L0_set)
      {
        int L0=Lk0.irrep;
        int kappa0_max=Lk0.tag;
        for(int kappa0=1; kappa0<=kappa0_max; ++kappa0)
          for(int Lp=etap%2; Lp<=etap; Lp+=2)
            for(int L=eta%2; L<=eta; L+=2)
              {
                if(not am::AllowedTriangle(L,L0,Lp))
                  continue;
                RelativeStateLabelsNLST bra_nlst(etap,Lp,Sp,Tp);
                RelativeStateLabelsNLST ket_nlst(eta,L,S,T);
                RelativeBraketNLST braket_nlst(L0,S0,T0,bra_nlst,ket_nlst);
                branched_rel_unit_tensors[braket_nlst]
                  +=u3::WCached(w_cache,u3::SU3(eta,0),1,L,x0,kappa0,L0,u3::SU3(etap,0),1,Lp,1);
              }
      }
  }


  void RelativeToCMLST(
    int Nmax, 
    const std::map<RelativeBraketNLST,double>& branched_rel_unit_tensors,
    std::map<u3shell::RelativeCMBraketNLST,double>& rel_cm_lst_map)
  // Coupling on center of mass as SO(3) level
  {
    for(auto it=branched_rel_unit_tensors.begin(); it!=branched_rel_unit_tensors.end(); ++it)
      {
        int L0,Lr,Lrp,N,Np;
        HalfInt S,T,Sp,Tp,S0,T0;
        RelativeStateLabelsNLST ket_rel,bra_rel;
        RelativeBraketNLST rel_tensor=it->first;
        double rme_rel=it->second;
        std::tie(L0,S0,T0,bra_rel,ket_rel)=rel_tensor;
        std::tie(Np,Lrp,Sp,Tp)=bra_rel;
        std::tie(N,Lr,S,T)=ket_rel;
        int Ncm_max=Nmax-std::max(Np,N);
        for(int Ncm=0; Ncm<=Ncm_max; ++Ncm)
          for(int Lcm=Ncm%2; Lcm<=Ncm; Lcm+=2)
            for(int Lp=abs(Lrp-Lcm); Lp<=(Lrp+Lcm); ++Lp)
              for(int L=abs(Lr-Lcm); L<=(Lr+Lcm); ++L)
                  {
                    u3shell::RelativeCMStateLabelsNLST bra_rel_cm(Np,Lrp,Ncm,Lcm,Lp,Sp,Tp);
                    u3shell::RelativeCMStateLabelsNLST ket_rel_cm(N,Lr,Ncm,Lcm,L,S,T);
                    u3shell::RelativeCMBraketNLST braket_rel_cm(L0,S0,T0,bra_rel_cm,ket_rel_cm);
                    rel_cm_lst_map[braket_rel_cm]
                      =am::Unitary9J(L0,0,L0,Lr,Lcm,L,Lrp,Lcm,Lp)
                      *parity(Lr+Lrp+L+Lp)*rme_rel;
                  }
      }
  }

  void RelativeToCMU3ST(int Nmax,  
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
    std::vector<RelativeCMUnitTensorCache>& unit_tensor_rel_cm_expansions
    )
  // Coupling on center of mass at U(3) level
  {
    int N0, Np,N, kappa0,L0;
    HalfInt S0,T0,Sp,Tp,S,T;
    u3::SU3 x0;
    u3shell::RelativeUnitTensorLabelsU3ST tensor;
    for(auto tensor:relative_unit_tensors)
      {
        RelativeCMUnitTensorCache relative_cm_u3st_map;
        std::tie(x0,S0,T0,Np,Sp,Tp,N,S,T)=tensor.FlatKey();
        u3::SU3 xr(N,0);
        u3::SU3 xrp(Np,0);
        for(int Ncm=0; Ncm<=Nmax; Ncm++)
          {
            u3::SU3 x_cm(Ncm,0);
            MultiplicityTagged<u3::SU3>::vector x_set=u3::KroneckerProduct(xr,x_cm);
            MultiplicityTagged<u3::SU3>::vector xp_set=u3::KroneckerProduct(xrp,x_cm);
            for(auto ip : xp_set)
              for(auto i: x_set)
                {
                  u3::SU3 x(i.irrep);
                  u3::SU3 xp(ip.irrep);
                  int rho0_max=u3::OuterMultiplicity(x,x0,xp);
                  u3shell::RelativeCMStateLabelsU3ST bra(Np,Ncm,xp,Sp,Tp);
                  u3shell::RelativeCMStateLabelsU3ST ket(N,Ncm,x,S,T);

                  for(int rho0=1; rho0<=rho0_max; ++rho0)
                  {
                    u3shell::RelativeCMUnitTensorLabelsU3ST tensor_cm(x0,S0,T0,rho0,bra,ket);
                    relative_cm_u3st_map[tensor_cm]
                      =//parity(x.lambda()+x.mu()+xp.lambda()+xp.mu()+Np+N)
                        u3::U(x0,xr,xp,x_cm,xrp,1,1,x,1,rho0);
                  }
                }
          }
        unit_tensor_rel_cm_expansions.push_back(relative_cm_u3st_map);
      }
  }

  void RelativeUnitTensorToRelativeCMUnitTensorU3ST(int Nmax,  
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensors,
    u3::WCoefCache& w_cache,
    RelativeCMExpansion& unit_relative_cm_map)
  // Coupling on center of mass by branching and re-upcoupling
  {
    for(int i=0; i<relative_unit_tensors.size(); ++i)
      {
        const u3shell::RelativeUnitTensorLabelsU3ST& tensor=relative_unit_tensors[i];

        std::map<RelativeBraketNLST,double> branched_rel_unit_tensors;
        BranchNLST(tensor,w_cache,branched_rel_unit_tensors);

        std::map<u3shell::RelativeCMBraketNLST,double> rel_cm_lst_map;
        RelativeToCMLST(Nmax, branched_rel_unit_tensors,rel_cm_lst_map);
        
        RelativeCMUnitTensorCache rel_cm_u3st_map;
        u3shell::UpcoupleCMU3ST(rel_cm_lst_map,w_cache,rel_cm_u3st_map);

        unit_relative_cm_map[tensor]=rel_cm_u3st_map;
      }
  }


}