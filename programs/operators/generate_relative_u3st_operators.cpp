/****************************************************************
  generate_relative_u3st_operators.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/14/17 (aem,mac): Created.
  2/21/17 (aem,mac): Update input parsing.  Add parsing checks.
****************************************************************/
#include <iostream>
#include <fstream>
#include "boost/math/constants/constants.hpp"
#include "cppformat/format.h"
#include "mcutils/parsing.h"

#include "u3shell/upcoupling.h"

// Checked agains mfdn for R2intr, Tintr, Nintr


// TODO (mac): update function names to accurately reflect that they
// are imperatives, not fruitful, e.g., AppendOperatorNintr...

double zero_threshold=10e-6;
namespace u3shell
{

  void Id(int Nmax, u3shell::RelativeRMEsU3ST& Id_operator, int A, double coef=1.0, bool moshinsky_convention=false)
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
              int L0=0;
              std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(relative_unit_tensor,kappa0,L0);
              Id_operator[key]+=2./(A*(A-1))*coef;
            }
  }


  void Nintr(int Nmax, u3shell::RelativeRMEsU3ST& Nrel_operator, int A, double coef=1.0, bool moshinsky_convention=false)
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
              int L0=0;
              std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(relative_unit_tensor,kappa0,L0);
              Nrel_operator[key]+=(2.*N)/A*coef;
            }
  }

  void Spin(int Nmax, u3shell::RelativeRMEsU3ST& Spin_operator, int A, double coef=1.0, bool moshinsky_convention=false)
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
              int L0=0;
              std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key(relative_unit_tensor,kappa0,L0);
              Spin_operator[key]+=(2*S*(S+1.))/A/(A-1)*coef;
            }
  }

  void k2intr(int Nmax,u3shell::RelativeRMEsU3ST& K2intr, int A, double coef=1.0, bool moshinsky_convention=false)
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
              K2intr[key]+=-sqrt(1.5)*Brme*intrinsic_factor*coef;

            // Hrel term
            Np=N;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Nrme=N+3/2.;
            K2intr[key]+=Nrme*intrinsic_factor*coef;

            // Arel term
            Np=N+2;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Arme=u3shell::RelativeSp3rRaisingOperator(bra,ket);
            if (fabs(Arme)>zero_threshold)
              K2intr[key]+=-sqrt(1.5)*Arme*intrinsic_factor*coef;
          }
  }

  void Lintr(int Nmax,u3shell::RelativeRMEsU3ST& Lintr, int A, double coef=1.0, bool moshinsky_convention=false)
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
          
            double Lrme=std::sqrt(4./3*(N*N+3*N));

            if(fabs(Lrme)>zero_threshold)
              Lintr[key]+=Lrme*intrinsic_factor*coef;
          }
  }


  void Qintr(int Nmax,u3shell::RelativeRMEsU3ST& Qintr, int A, double coef=1.0, bool moshinsky_convention=false)
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
              double Brme=u3shell::RelativeSp3rLoweringOperator(bra,ket);

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
              double Arme=u3shell::RelativeSp3rRaisingOperator(bra,ket);
              if (fabs(Arme)>zero_threshold)
                Qintr[key]+=sqrt(3)*Arme*intrinsic_factor*coef;
            }
          }
  }

  void r2intr(int Nmax,u3shell::RelativeRMEsU3ST& R2intr, int A, double coef=1.0, bool moshinsky_convention=false)
  {
    u3shell::RelativeStateLabelsU3ST bra,ket;
    int Np, rho_max;
    u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
    std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key;
    int kappa0=1; 
    int L0=0;

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
              double Brme=u3shell::RelativeSp3rLoweringOperator(bra,ket);
              if (fabs(Brme)>zero_threshold)
               R2intr[key]+=sqrt(1.5)*Brme*intrinsic_factor*coef;
            }
            // Crel term
            Np=N;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,0),0,0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Nrme=N+3/2.;
            R2intr[key]+=Nrme*intrinsic_factor*coef;

            // Arel term
            Np=N+2;
            if(Np<=Nmax)
            {
              bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
              relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,0,bra,ket);
              key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
              double Arme=u3shell::RelativeSp3rRaisingOperator(bra,ket);
              if (fabs(Arme)>zero_threshold)
                R2intr[key]+=sqrt(1.5)*Arme*intrinsic_factor*coef;
            }
          }
  }

  void Tintr(int Nmax,u3shell::RelativeRMEsU3ST& Tintr, int A, double hbar_omega, double coef=1.0, bool moshinsky_convention=false)
  {
    coef*=hbar_omega/(4*(KroneckerDelta(moshinsky_convention,false)));
    k2intr(Nmax,Tintr, A, coef, moshinsky_convention);
  }

  void Interaction(int Nmax, int Jmax, int J0, int T0, int g0, std::string& interaction_filename,
                    u3shell::RelativeRMEsU3ST& interaction_u3st, int A, double coef=1.0)
    {
      
      // Read in the interaction from file
      basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);
      std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
      std::array<basis::MatrixVector,3> isospin_component_matrices_lsjt;
      basis::OperatorLabelsJT operator_labels;
      
      basis::ReadRelativeOperatorLSJT(
        interaction_filename,relative_space_lsjt,operator_labels,
        isospin_component_sectors_lsjt, isospin_component_matrices_lsjt, true
        );

      // Multiplying interaction by coefficient 
      for(int T0=0; T0<=2; ++T0)
        for(auto& matrix : isospin_component_matrices_lsjt[T0])
          matrix=coef*matrix;

      // upcouple interaction
      u3shell::Upcoupling(
        relative_space_lsjt,
        isospin_component_sectors_lsjt,
        isospin_component_matrices_lsjt,
        J0, g0, T0,Nmax, interaction_u3st);


    }

} // end namespace


// load file
//
// hbar_omega
// operator_type coef (if INT, also Jmax, J0, T0, g0)
// operator_type coef
// operator_type coef

int main(int argc, char **argv)
{
  // constants 
  double hbarc=197.327; //MeVfm
  double mc2=938.92; //MeV
  double alpha=1/137.036; 
  const double pi=boost::math::constants::pi<double>();

  if(argc<4)
    std::cout<<"Syntax: A Nmax N1B <operator_filename_base> "<<std::endl;
  
  u3::U3CoefInit();

  int A=std::stoi(argv[1]);
  int Nmax=std::stoi(argv[2]);
  int N1B=std::stoi(argv[3]);
  std::string operator_filename = argv[4];

  int Jmax, J0, T0, g0;
  std::string interaction_filename;

  std::string operator_type;
  std::ifstream is(fmt::format("{}.load",operator_filename));
  if(not is)
    std::cout<<"is not open"<<std::endl;
  assert(is);

  double hbar_omega;
  double b2;

  u3shell::RelativeRMEsU3ST Operator;
  double coef;

  std::string line;
  int line_count=0;
  while(std::getline(is,line))
  {
    ++line_count;
    std::istringstream line_stream(line);
    
    if(line_count==1)
      {
        line_stream >> hbar_omega;
        ParsingCheck(line_stream,line_count,line);
        b2=hbarc*hbarc/mc2/hbar_omega;
        std::cout<<"bsqr= "<<b2<<std::endl;
        std::cout<<"pi= "<<pi<<std::endl;
        continue;
      }

    line_stream >> operator_type >> coef;
    ParsingCheck(line_stream,line_count,line);

    std::cout<<operator_type<<"  "<<coef<<std::endl;

    if(operator_type=="ID") 
      u3shell::Id(Nmax+2*N1B,Operator, A, coef);
    else if(operator_type=="Nintr") 
      u3shell::Nintr(Nmax+2*N1B,Operator, A, coef);
    else if(operator_type=="Spin") 
      u3shell::Spin(Nmax+2*N1B,Operator, A, coef);
    else if(operator_type=="r2intr")
      // Factor of 1/A on coef is to compute
      // the mean square of the radius   
      u3shell::r2intr(Nmax+2*N1B,Operator, A, coef*b2/A);
    else if(operator_type=="k2intr")
      u3shell::k2intr(Nmax+2*N1B,Operator, A, coef/b2);
    else if(operator_type=="Lintr")
      u3shell::Lintr(Nmax+2*N1B,Operator, A, coef);
    else if(operator_type=="Qintr")
      u3shell::Qintr(Nmax+2*N1B,Operator, A, sqrt(5./(16*pi))*coef*b2);
    else if(operator_type=="Tintr")
      u3shell::Tintr(Nmax+2*N1B,Operator, A, hbar_omega, coef);
    else if(operator_type=="INT")
      {
        line_stream >> Jmax >> J0 >> T0 >> g0 >> interaction_filename;  
        ParsingCheck(line_stream,line_count,line);
        Interaction(Nmax+2*N1B, Jmax, J0, T0, g0, interaction_filename,Operator,A, coef);
      }
    else if(operator_type=="COUL")
      {
        line_stream >> Jmax >> J0 >> T0 >> g0 >> interaction_filename;  
        ParsingCheck(line_stream,line_count,line);
        std::cout<<alpha*sqrt(mc2*hbar_omega)<<std::endl;
        coef*=alpha*sqrt(mc2*hbar_omega);
        // coef*=alpha*sqrt(mc2*hbar_omega);
        Interaction(Nmax+2*N1B, Jmax, J0, T0, g0, interaction_filename,Operator,A, coef);
      }
    else
      {
        std::cout<<fmt::format("{} is not a valid operator type",operator_type)<<std::endl
               <<"The allowed operator types are:"<<std::endl
               <<"    Id, Lintr, k2intr, r2intr, Nintr, Qintr, Tintr, INT, COUL"<<std::endl;

        std::exit(EXIT_FAILURE);
      }
  }
  is.close();
  
  // Writing operator to file 
  std::ofstream os(fmt::format("{}_hw{:.1f}_Nmax{:02d}_u3st.dat",operator_filename,hbar_omega,Nmax));
  u3shell::WriteRelativeOperatorU3ST(os, Operator);  
}



