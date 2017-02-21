/****************************************************************
  generate_relative_u3st_operators.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/14/17 (aem,mac): Created.
****************************************************************/
#include <iostream>
#include <fstream>
#include "cppformat/format.h"
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

  void Lintr(int Nmax,u3shell::RelativeRMEsU3ST& Lintr, int A, bool moshinsky_convention=true)
    // Need to check intrinsic 
  {
    int kappa0=1; 
    int L0=1;

    double intrinsic_factor=2./A;

    for(int N=0; N<=Nmax; N++)
      for(int S=0; S<=1; ++S)
        for(int T=0; T<=1; ++T)
          if((N+S+T)%2==1)
          {
            u3shell::RelativeStateLabelsU3ST ket(N,S,T); 
            u3shell::RelativeStateLabelsU3ST bra(N,S,T);
            u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor(u3::SU3(1,1),0,0,bra,ket);
            std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(relative_unit_tensor,kappa0,L0);
          
            double Lrme=std::sqrt(2./3*(N*N+3*N));

            if(fabs(Lrme)>zero_threshold)
              Lintr[key]+=Lrme*intrinsic_factor;
          }
  }

  void Qintr(int Nmax,u3shell::RelativeRMEsU3ST& Qintr, int A, bool moshinsky_convention=true)
  {
    u3shell::RelativeStateLabelsU3ST bra,ket;
    int Np, rho_max;
    u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
    std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key;
    int kappa0=1; 
    int L0=2;

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
              Qintr[key]+=sqrt(3)*Brme*intrinsic_factor;

            // Hrel term
            Np=N;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Nrme=N+3/2.;
            Qintr[key]+=sqrt(3)*Nrme*intrinsic_factor;

            // Arel term
            Np=N+2;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Arme=u3shell::RelativeSp3rRaisingOperator(bra,ket);
            if (fabs(Arme)>zero_threshold)
              Qintr[key]+=sqrt(3)*Arme*intrinsic_factor;
          }

  }

  void Interaction(int Nmax, int Jmax, int J0, int T0, int g0, std::string& interaction_filename,
                    u3shell::RelativeRMEsU3ST& Interaction_u3st, int A)
    {

      basis::RelativeSectorsLSJT relative_lsjt_sectors;
      basis::RelativeSpaceLSJT relative_lsjt_space(Nmax, Jmax);
      basis::MatrixVector sector_vector;
      u3shell::GetInteractionMatrix(interaction_filename, relative_lsjt_space,relative_lsjt_sectors,sector_vector);

      //upcouple to LST
      std::map<u3shell::RelativeSectorNLST,Eigen::MatrixXd> Interaction_nlst;
      u3shell::UpcouplingNLST(relative_lsjt_space,relative_lsjt_sectors,sector_vector,J0,g0,T0,Nmax,Interaction_nlst);

      // Upcouple to U(3) level
      u3shell::UpcouplingU3ST(Interaction_nlst, T0, Nmax, Interaction_u3st);
    }

} // end namespace

int main(int argc, char **argv)
{

  if(argc<4)
    std::cout<<"Syntax: A Nmax N1B <load_filename> "<<std::endl;
  
  u3::U3CoefInit();

  int A=std::stoi(argv[1]);
  int Nmax=std::stoi(argv[2]);
  int N1B=std::stoi(argv[3]);
  std::string load_filename = argv[4];

  int Jmax, J0, T0, g0;
  std::string interaction_filename;

  std::string operator_type;
  std::ifstream is(load_filename);
  while(!is.eof())
  {
    is >> operator_type;

    if(operator_type=="INT")       
      is >> Jmax >> J0 >> T0 >> g0 >> interaction_filename;

    u3shell::RelativeRMEsU3ST Operator;

    if(operator_type=="Nintr" || operator_type=="ALL") 
      u3shell::Nintr(Nmax+2*N1B,Operator, A);
    
    else if(operator_type=="K2intr" || operator_type=="ALL")
      u3shell::K2intr(Nmax+2*N1B,Operator, A);

    else if(operator_type=="Lintr" || operator_type=="ALL")
      u3shell::Lintr(Nmax+2*N1B,Operator, A);

    else if(operator_type=="Qintr" || operator_type=="ALL")
      u3shell::Qintr(Nmax+2*N1B,Operator, A);

    else if(operator_type=="INT")
      Interaction(Nmax, Jmax, J0, T0, g0, interaction_filename,Operator,A);

    else
      {
        std::cout<<fmt::format("{} is not a valid operator type",operator_type)<<std::endl
               <<"The allowed operator types are:"<<std::endl
               <<"    Lintr, K2intr, Nintr, Qintr, INT, ALL"<<std::endl
               <<"ALL includes everything but INT"<<std::endl;

        exit;
      }

    // Writing operator to file 
    std::string operator_filename;
    if(operator_type=="INT")
      operator_filename=fmt::format("{}_u3st.dat",interaction_filename);
    else
      operator_filename=fmt::format("{}_u3st.dat",operator_type);
    
    std::ofstream os(operator_filename);
    u3shell::WriteRelativeOperatorU3ST(os, Operator);  

  }
    
  is.close();
}



