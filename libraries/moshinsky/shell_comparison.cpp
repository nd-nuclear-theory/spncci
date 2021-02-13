/****************************************************************
  shell_comparison.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

****************************************************************/
#include <cmath>
#include "fmt/format.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "am/wigner_gsl.h"
#include "sp3rlib/u3coef.h"
#include "moshinsky/moshinsky_xform.h"
#include "u3shell/upcoupling.h"
extern double zero_threshold;

namespace u3shell
{
  void Q(int Nmax,u3shell::RelativeRMEsU3ST& Qintr, int A, double coef=1.0, bool moshinsky_convention=false)
  {
    u3shell::RelativeStateLabelsU3ST bra,ket;
    int Np, rho_max;
    u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
    std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key;
    int kappa0=1;
    int L0=2;

    double intrinsic_factor=2./A;
    for(int N=0; N<=Nmax; N++)
      for(int S=0; S<=1; ++S)
        for(int T=0; T<=1; ++T)
          if((N+S+T)%2==1)
          {
            ket=u3shell::RelativeStateLabelsU3ST(N,S,T);

            // Brel term
            Np=N-2;
            if(Np>=0)
            {
              bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
              relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,2),0,0,bra,ket);
              key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
              double Brme=u3shell::Brel(bra,ket);

              if (fabs(Brme)>zero_threshold)
                Qintr[key]+=sqrt(3)*Brme*intrinsic_factor*coef;
            }
            // Crel term
            Np=N;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(1,1),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Crme=std::sqrt(4.*(N*N+3*N)/3);

            if(fabs(Crme)>zero_threshold)
              Qintr[key]+=std::sqrt(3)*Crme*intrinsic_factor*coef;

            // Arel term
            Np=N+2;
            if(Np<=Nmax)
            {
              bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
              relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,0,bra,ket);
              key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
              double Arme=u3shell::Arel(bra,ket);
              if (fabs(Arme)>zero_threshold)
                Qintr[key]+=sqrt(3)*Arme*intrinsic_factor*coef;
            }
          }
  }

  typedef std::tuple<int,int,HalfInt,HalfInt,HalfInt>RelativeSubspaceLabelNLSJT;
  typedef std::tuple<HalfInt,HalfInt,RelativeSubspaceLabelNLSJT,RelativeSubspaceLabelNLSJT> RelativeBraketNLSJT;

  void BranchU3STtoNLSTJ(
    const u3shell::RelativeRMEsU3ST& rmes_u3st,
    HalfInt J0,
    int Nmax
    )
    {
      std::map<RelativeBraketNLSJT,double>rmes_lsjt;

      // Read in the interaction from file
      basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Nmax+2);
      std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
      std::array<basis::MatrixVector,3> isospin_component_matrices_lsjt;
      basis::OperatorLabelsJT operator_labels;

      for(auto it=rmes_u3st.begin(); it!=rmes_u3st.end(); ++it)
        {
          u3::SU3 x0;

          u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels;
          int kappa0,L0,Np,N;
          HalfInt Sp,Tp,S,T,S0,T0;
          std::tie(unit_tensor_labels,kappa0,L0)=it->first;
          std::tie(x0,S0,T0,Np,Sp,Tp,N,S,T)=unit_tensor_labels.FlatKey();
          double rme=it->second;
          for(int Lp=Np; Lp>=0; Lp=Lp-2)
            for(int L=N;L>=0; L=L-2)
              {
                int np((Np-Lp)/2);
                int n((N-L)/2);
                double rme_lst=u3::W(u3::SU3(N,0),1,L,x0,kappa0,L0,u3::SU3(Np,0),1,Lp,1)*rme;

                // Branching to J
                for(HalfInt Jp=abs(Lp-Sp); Jp<=Lp+Sp; ++Jp)
                  for(HalfInt J=abs(L-S); J<=L+S; ++J)
                    {
                      RelativeSubspaceLabelNLSJT bra(Np,Lp,Sp,Jp,Tp);
                      RelativeSubspaceLabelNLSJT ket(N,L,S,J,T);
                      RelativeBraketNLSJT braket(J0,T0,bra,ket);
                      rmes_lsjt[braket]+=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp)*rme_lst;
                    }
              }
          }
        for(auto it=rmes_lsjt.begin(); it!=rmes_lsjt.end(); ++it)
          {
            int L,Lp,L0,Np,N;
            HalfInt Sp,Tp,Jp,S,J,T,S0,T0;

            RelativeSubspaceLabelNLSJT bra,ket;
            std::tie(J0,T0,bra,ket)=it->first;
            std::tie( Np,Lp,Sp,Jp,Tp)=bra;
            std::tie( N,L,S,J,T)=ket;
            std::cout<<fmt::format("{} {}  {} {} {} {} {}   {} {} {} {} {}  {}",
              J0,T0,Np,Lp,Sp,Jp,Tp,N,L,S,J,T,it->second)
            <<std::endl;
          }
    }

}//namespace

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  u3shell::RelativeRMEsU3ST rmes_u3st;
  int J0=2;
  int Nmax=2;

  u3shell::Q(Nmax,rmes_u3st,2);
  u3shell::BranchU3STtoNLSTJ(rmes_u3st,J0,Nmax);

}
