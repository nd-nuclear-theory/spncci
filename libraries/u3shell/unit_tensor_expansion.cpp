/****************************************************************
  unit_tensor_expansion.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/unit_tensor_expansion.h"

#include <array>
#include <cmath>

#include "cppformat/format.h"
#include "am/am.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/moshinsky.h"

namespace u3shell
{
	//NOT CORRECT
  // void NumberOperatorUnitTensorExpansion(
  //       int Nmin, int Nmax, 
  //       u3shell::TwoBodyUnitTensorCoefficientsU3ST& N_operator_expansion, 
  //       int A
  //     )
  // {
  //   u3shell::OperatorLabelsU3ST NumberOperator(0,u3::SU3(0,0),0,0,0);
  //   u3shell::TwoBodySpaceU3ST space(Nmax);
  //   u3shell::TwoBodySectorsU3ST two_body_sectors(space,NumberOperator);

  //   for(int sector_index=0; sector_index<two_body_sectors.size(); ++sector_index)
  //     {
  //       const u3shell::TwoBodySectorsU3ST::SectorType& sector = two_body_sectors.GetSector(sector_index);
  //       const u3shell::TwoBodySectorsU3ST::SubspaceType& ket_subspace = sector.ket_subspace();

  //       // extract subspace labels, same for both bra and ket
  //       u3::U3 omega;
  //       HalfInt S, T;
  //       int g;
  //       std::tie(omega,S,T,g) =  ket_subspace.GetSubspaceLabels();
  //       int eta1,eta2;
  //       for(int index=0; index<ket_subspace.size(); ++index)
  //         {
  //           std::tie(eta1,eta2)=ket_subspace.GetStateLabels(index);

  //           u3shell::TwoBodyStateLabelsU3ST labels(eta1,eta2,omega.SU3(),S,T);
  //           u3shell::TwoBodyUnitTensorLabelsU3ST 
  //             unit_tensor_labels(NumberOperator,1,labels,labels);
  //           N_operator_expansion[unit_tensor_labels]+=(eta1+eta2);
  //         }
  //     }
  // }

  void BrelRelativeUnitTensorExpansion(int Nmin, int Nmax, u3shell::RelativeUnitTensorCoefficientsU3ST& Brel_operator, int A)
  {
    for(int N=Nmin; N<=Nmax; N+=2)
      for(int S=0; S<=1; ++S)
        for(int T=0; T<=1; ++T)
          if((N+S+T)%2==1)
            {
              int Np=N-2;
              int rho_max=u3::OuterMultiplicity(u3::SU3(N,0), u3::SU3(0,2), u3::SU3(Np,0));
              if(rho_max==0)
                continue;

              u3shell::RelativeStateLabelsU3ST bra(Np,S,T);
              u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
              u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(u3::SU3(0,2),0,0,bra,ket);
              double rme=u3shell::RelativeSp3rLoweringOperator(bra,ket);
              if (fabs(rme)>10e-10)
                Brel_operator[relative_unit_tensor]+=2.*rme/(A*(A-1));
            }
  }
  void NrelRelativeUnitTensorExpansion(int Nmin, int Nmax, u3shell::RelativeUnitTensorCoefficientsU3ST& Nrel_operator, int A)
  {
    for (int N=Nmin; N<=Nmax; N++)
      for (int S=0;S<=1; S++)
        for (int T=0;T<=1; T++)
          if ((N+S+T)%2==1)
            {
              u3shell::RelativeStateLabelsU3ST bra(N,S,T);
              u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
              u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(u3::SU3(0,0),0,0,bra,ket);
              double rme=u3shell::RelativeNumberOperator(bra,ket);
              if (fabs(rme)>10e-10)
                Nrel_operator[relative_unit_tensor]+=N;
            }
  }

  void IdentityRelativeUnitTensorExpansion(int Nmin, int Nmax, u3shell::RelativeUnitTensorCoefficientsU3ST& Nrel_operator, int A)
  {
    for (int N=Nmin; N<=Nmax; N++)
      for (int S=0;S<=1; S++)
        for (int T=0;T<=1; T++)
          if ((N+S+T)%2==1)
            {
              u3shell::RelativeStateLabelsU3ST bra(N,S,T);
              u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
              u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(u3::SU3(0,0),0,0,bra,ket);
              Nrel_operator[relative_unit_tensor]+=2./(A*(A-1));
            }
  }


}