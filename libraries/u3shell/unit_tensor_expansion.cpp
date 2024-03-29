/****************************************************************
  unit_tensor_expansion.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include "u3shell/unit_tensor_expansion.h"

#include <array>
#include <cmath>

#include "fmt/format.h"
#include "am/am.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"

// currently defined in two_body_operators.cpp
extern double zero_threshold;

namespace u3shell
{

  void BrelRelativeUnitTensorExpansion(int Nmin, int Nmax, 
        u3shell::RelativeUnitTensorCoefficientsU3ST& Brel_operator, int A)
  {
    for(int N=Nmin; N<=Nmax; N++)
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
              double rme=u3shell::Brel(bra,ket);
              if (fabs(rme)>zero_threshold)
                Brel_operator[relative_unit_tensor]+=2.*rme/A;
            }
  }

void ArelRelativeUnitTensorExpansion(int Nmin, int Nmax, 
      u3shell::RelativeUnitTensorCoefficientsU3ST& Arel_operator, int A)
  {
    for(int N=Nmin; N<=Nmax; N++)
      for(int S=0; S<=1; ++S)
        for(int T=0; T<=1; ++T)
          if((N+S+T)%2==1)
            {
              int Np=N+2;
              int rho_max=u3::OuterMultiplicity(u3::SU3(N,0), u3::SU3(2,0), u3::SU3(Np,0));
              if(rho_max==0)
                continue;

              u3shell::RelativeStateLabelsU3ST bra(Np,S,T);
              u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
              u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(u3::SU3(2,0),0,0,bra,ket);
              double rme=u3shell::Arel(bra,ket);
              if (fabs(rme)>zero_threshold)
                Arel_operator[relative_unit_tensor]+=2.*rme/A;
            }
  }

  void NrelRelativeUnitTensorExpansion(int Nmin, int Nmax, 
        u3shell::RelativeUnitTensorCoefficientsU3ST& Nrel_operator, int A)
  {
    for (int N=Nmin; N<=Nmax; N++)
      for (int S=0;S<=1; S++)
        for (int T=0;T<=1; T++)
          if ((N+S+T)%2==1)
            {
              u3shell::RelativeStateLabelsU3ST bra(N,S,T);
              u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
              u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(u3::SU3(0,0),0,0,bra,ket);
              double rme=u3shell::Nrel(bra,ket);
              if (fabs(rme)>zero_threshold)
                Nrel_operator[relative_unit_tensor]+=(2.*N)/A;
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

void TrelRelativeUnitTensorExpansion(int Nmin,int Nmax,u3shell::RelativeUnitTensorCoefficientsU3ST& Trel_operator,int A)
{
  u3shell::RelativeStateLabelsU3ST bra,ket;
  int Np, rho_max;
  u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
  for(int N=Nmin; N<=Nmax; N++)
  for(int S=0; S<=1; ++S)
    for(int T=0; T<=1; ++T)
      if((N+S+T)%2==1)
        {
          ket=u3shell::RelativeStateLabelsU3ST(N,S,T); 

          // Brel term
          Np=N-2;
           bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
          relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,2),0,0,bra,ket);
          double Brme=u3shell::Brel(bra,ket);
          if (fabs(Brme)>zero_threshold)
            // Division by 2 for comparison with Tomas Trel (convert from moshinsky to mechanics convention)
            Trel_operator[relative_unit_tensor]+=-sqrt(1.5)*Brme/2.;
          
          // Hrel term
          Np=N;
          bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
          relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket);
          double Nrme=N+3/2.;
          if (fabs(Nrme)>zero_threshold)
            // Division by 2 for comparison with Tomas Trel (convert from moshinsky to mechanics convention)
            Trel_operator[relative_unit_tensor]+=Nrme/2.;

          // Arel term
          Np=N+2;
          bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
          relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,0,bra,ket);
          double Arme=u3shell::Arel(bra,ket);
          if (fabs(Arme)>zero_threshold)
            // Division by 2 for comparison with Tomas Trel (convert from moshinsky to mechanics convention)
            Trel_operator[relative_unit_tensor]+=-sqrt(1.5)*Arme/2.;
        }

}




}
