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
              double rme=u3shell::RelativeSp3rLoweringOperator(bra,ket);
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
              double rme=u3shell::RelativeSp3rRaisingOperator(bra,ket);
              if (fabs(rme)>zero_threshold)
                Arel_operator[relative_unit_tensor]+=2.*rme/A;
            }
  }


void RaisingPolynomialRelativeUnitTensorExpansion(const u3::U3& n0, int Nmin, int Nmax, 
      u3shell::RelativeUnitTensorCoefficientsU3ST& Prel_operator, int A)
  {
    // Get B coefficients
    sp3r::BCoefCache bcoef_cache;
    sp3r::GenerateBCoefCache(bcoef_cache,Nmax);
    // std::vector<u3::U3> raising_polynomials=sp3r::RaisingPolynomialLabels(Nmax);
    // Get list of raising polynomials 
    // for(auto n0 : raising_polynomials)
    for(int S=0; S<=1; ++S)
      for(int T=0; T<=1; ++T)
        for(int N=Nmin; N<=Nmax; N++)
          {
            int Np=int(N+n0.N());
            // Check parity constraint on bra and ket
            if(((N+S+T)%2==0)||((Np+S+T)%2==0))
              continue;
            //Check outer multiplciity of SU(3) coupling
            int rho_max=u3::OuterMultiplicity(u3::SU3(N,0), n0.SU3(), u3::SU3(Np,0));
            if(rho_max==0)
              continue;
          
            u3shell::RelativeStateLabelsU3ST bra(Np,S,T);
            u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
            u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(n0.SU3(),0,0,bra,ket);
            double rme=u3shell::RelativeSp3rRaisingPolynomial(n0,bra,ket,bcoef_cache);
            std::cout<<fmt::format("{}  {} {} {}   {} {} {}  {}",n0.Str(),Np,S,T, N, S, T,rme)<<std::endl;
            if (fabs(rme)>zero_threshold)
              Prel_operator[relative_unit_tensor]+=rme;
          }
  }


  void NintrRelativeUnitTensorExpansion(int Nmin, int Nmax, 
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
              double rme=u3shell::RelativeNumberOperator(bra,ket);
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
          // rho_max=u3::OuterMultiplicity(u3::SU3(N,0), u3::SU3(0,2), u3::SU3(Np,0));
          // if(rho_max==0)
          //   continue;
          bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
          relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,2),0,0,bra,ket);
          double Brme=u3shell::RelativeSp3rLoweringOperator(bra,ket);
          if (fabs(Brme)>zero_threshold)
            // Division by 2 for comparison with Tomas Trel
            Trel_operator[relative_unit_tensor]+=-sqrt(1.5)*Brme/2.;
          // Hrel term
          Np=N;
          // rho_max=u3::OuterMultiplicity(u3::SU3(N,0), u3::SU3(0,0), u3::SU3(Np,0));
          // if(rho_max==0)
          //   continue;
          bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
          relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket);
          double Nrme=N+3/2.;
          if (fabs(Nrme)>zero_threshold)
            // Division by 2 for comparison with Tomas Trel
            Trel_operator[relative_unit_tensor]+=Nrme/2.;

          // Arel term
          Np=N+2;
          // rho_max=u3::OuterMultiplicity(u3::SU3(N,0), u3::SU3(2,0), u3::SU3(Np,0));
          // if(rho_max==0)
          //   continue;
          bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
          relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,0,bra,ket);
          double Arme=u3shell::RelativeSp3rRaisingOperator(bra,ket);
          if (fabs(Arme)>zero_threshold)
            // Division by 2 for comparison with Tomas Trel
            Trel_operator[relative_unit_tensor]+=-sqrt(1.5)*Arme/2.;
        }

}




}