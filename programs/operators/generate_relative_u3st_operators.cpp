/****************************************************************
  generate_relative_u3st_operators.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/14/17 (aem,mac): Created.
****************************************************************/
#include <fstream>
#include "u3shell/upcoupling.h"

double zero_threshold=10e-6;
namespace u3shell
{
  void Nintr(int Nmax, u3shell::RelativeRMEsU3ST& Nrel_operator,int A)
  {
    for (int N=0; N<=Nmax; N++)
      for (int S=0;S<=1; S++)
        for (int T=0;T<=1; T++)
          if ((N+S+T)%2==1)
            {
              u3shell::RelativeStateLabelsU3ST bra(N,S,T);
              u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
              u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(u3::SU3(0,0),0,0,bra,ket);
              int kappa0=1;
              int L0=1;
              std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(relative_unit_tensor,kappa0,L0);
              Nrel_operator[key]+=(2.*N)/A;
            }
  }

  void K2intr(int Nmax,u3shell::RelativeRMEsU3ST& K2intr, int A, bool moshinsky_convention=true)
  {
    u3shell::RelativeStateLabelsU3ST bra,ket;
    int Np, rho_max;
    u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
    std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key;
    int kappa0=1; 
    int L0=0;

    double intrinsic_factor=2.*(1.+KroneckerDelta(moshinsky_convention,false))/A;

    for(int N=0; N<=Nmax; N++)
    for(int S=0; S<=1; ++S)
      for(int T=0; T<=1; ++T)
        if((N+S+T)%2==1)
          {
            ket=u3shell::RelativeStateLabelsU3ST(N,S,T); 

            // Brel term
            Np=N-2;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,2),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Brme=u3shell::RelativeSp3rLoweringOperator(bra,ket);
            if (fabs(Brme)>zero_threshold)
              K2intr[key]+=-sqrt(1.5)*Brme*intrinsic_factor;

            // Hrel term
            Np=N;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Nrme=N+3/2.;
            K2intr[key]+=Nrme*intrinsic_factor;

            // Arel term
            Np=N+2;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Arme=u3shell::RelativeSp3rRaisingOperator(bra,ket);
            if (fabs(Arme)>zero_threshold)
              K2intr[key]+=-sqrt(1.5)*Arme*intrinsic_factor;
          }
  }
} // end namespace 
int main(int argc, char **argv)
{

  if(argc<5)
    std::cout<<"Syntax: Z  N  Nmax  N1B"<<std::endl;

  u3::U3CoefInit();

  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  int N1B=std::stoi(argv[4]);
  int A=N+Z;

  u3shell::RelativeRMEsU3ST Nintr;
  u3shell::Nintr(Nmax+2*N1B,Nintr, A);
  std::string nintr_filename="Nintr_u3st.dat";
  std::ofstream nintr_os(nintr_filename);
  u3shell::WriteRelativeOperatorU3ST(nintr_os, Nintr);  

  // Tintr=(hbar^2/(2m))*K2intr

  u3shell::RelativeRMEsU3ST Krel2;
  u3shell::K2intr(Nmax+2*N1B,Krel2, A);
  std::string krel2_filename="K2intr_u3st.dat";
  std::ofstream krel2_os(krel2_filename);
  u3shell::WriteRelativeOperatorU3ST(krel2_os, Krel2);  



}



