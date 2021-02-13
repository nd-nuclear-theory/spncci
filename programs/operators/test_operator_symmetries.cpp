/****************************************************************
  generate_relative_u3st_operators.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  2/14/17 (aem,mac): Created.
  2/21/17 (aem,mac): Update input parsing.  Add parsing checks.
  5/14/19 (aem): Updated   basis::OperatorLabelsJT to 
    basis::RelativeOperatorParametersLSJT for reading in operators
****************************************************************/
#include <iostream>
#include <fstream>
#include "boost/math/constants/constants.hpp"
#include "fmt/format.h"
#include "mcutils/parsing.h"

#include "u3shell/upcoupling.h"

// Checked agains mfdn for R2intr, Tintr, Nintr

double zero_threshold=10e-6;
namespace u3shell
{


  void Qintr(int Nmax,u3shell::RelativeRMEsU3ST& Qintr, int A, int T0,double coef=1.0, bool moshinsky_convention=false)
  {
    u3shell::RelativeStateLabelsU3ST bra,ket;
    int Np, rho_max;
    u3shell::RelativeUnitTensorLabelsU3ST relative_unit_tensor;
    std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> key;
    int kappa0=1; 
    int L0=2;

    std::cout<<"coefficient "<<coef<<std::endl;
    double intrinsic_factor=2./A;
    for(int N=0; N<=Nmax; N++)
      for(int S=0; S<=1; ++S)
        for(int T=0; T<=1; ++T)
          if((N+S+T)%2==1)
          {
            ket=u3shell::RelativeStateLabelsU3ST(N,S,T); 
            double isospin_coefficient=(T0==1)?2*sqrt(T*(T+1)):1;
            // Brel term
            Np=N-2;
            if(Np>=0)
            {
              bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
              relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(0,2),0,T0,bra,ket);
              key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
              double Brme=u3shell::Brel(bra,ket);


              if (fabs(Brme)>zero_threshold)
                Qintr[key]+=sqrt(3)*Brme*intrinsic_factor*coef*isospin_coefficient;
            }
            // Crel term
            Np=N;
            bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
            relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(1,1),0,T0,bra,ket);
            key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
            double Crme=std::sqrt(4.*(N*N+3*N)/3);
            
            if(fabs(Crme)>zero_threshold)
              Qintr[key]+=std::sqrt(3)*Crme*intrinsic_factor*coef*isospin_coefficient;

            // Arel term
            Np=N+2;
            if(Np<=Nmax)
            {
              bra=u3shell::RelativeStateLabelsU3ST(Np,S,T);
              relative_unit_tensor=u3shell::RelativeUnitTensorLabelsU3ST(u3::SU3(2,0),0,T0,bra,ket);
              key=std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>(relative_unit_tensor,kappa0,L0);
              double Arme=u3shell::Arel(bra,ket);
              if (fabs(Arme)>zero_threshold)
                Qintr[key]+=sqrt(3)*Arme*intrinsic_factor*coef*isospin_coefficient;
            }
          }
  }

  void Interaction(int Nmax, int Jmax, int J0, int T0, int g0, std::string& interaction_filename,
                    u3shell::RelativeRMEsU3ST& interaction_u3st, int A, double coef=1.0)
    {
      
      // Read in the interaction from file
      basis::RelativeSpaceLSJT relative_space_lsjt(Nmax, Jmax);
      std::array<basis::RelativeSectorsLSJT,3> isospin_component_sectors_lsjt;
      std::array<basis::MatrixVector,3> isospin_component_matrices_lsjt;
      basis::RelativeOperatorParametersLSJT operator_labels;
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


      void CheckInteraction(int A, int Nmax, int Jmax, int J0, int T0, int g0, std::string& interaction_filename)
        {

          double coef=1.0;
          u3shell::RelativeRMEsU3ST Operator;
            Interaction(Nmax, Jmax, J0, T0, g0, interaction_filename,Operator,A, coef);


          //Print out RMEs
          for(auto it=Operator.begin(); it!=Operator.end(); ++it)
            {

              u3shell::RelativeUnitTensorLabelsU3ST labels;
              int kappa0, L0;
              std::tie(labels,kappa0,L0)=it->first;
              double rme=it->second;

              u3::SU3 x0;
              HalfInt S0, T0, S, T, Sp,Tp;
              int Np,N;
              std::tie(x0,S0,T0,Np,Sp,Tp,N,S,T)=labels.FlatKey();

              u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
              u3shell::RelativeStateLabelsU3ST ket(N,S,T);
              u3shell::OperatorLabelsU3ST op_labels_conj(N-Np,u3::Conjugate(x0),S0,T0,g0);
              u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels_conj(op_labels_conj,ket,bra);
              std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> labels_conj(unit_tensor_labels_conj,kappa0,L0);
                   
              int phase=ParitySign(u3::ConjugationGrade(x0)+S0);

              if(x0.lambda()==x0.mu())
                phase=phase*ParitySign(L0);

              double rme_conj=std::sqrt(u3::dim(u3::SU3(N,0)))/std::sqrt(u3::dim(u3::SU3(Np,0)))
                              *Hat(S)*Hat(T)/Hat(Sp)/Hat(Tp)
                              *phase
                              *Operator[labels_conj];

              if(abs(rme-rme_conj)>1e-4)
              {
                std::cout<<"WRONG"<<std::endl;
                std::cout<<labels.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl;
                std::cout<<rme<<std::endl;
                std::cout<<rme_conj<<std::endl;
              }
            
              // std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>labels_conj()
            }
        }

} // end namespace

int main(int argc, char **argv)
{
  // constants 
  double hbarc=197.327; //MeVfm
  double mc2=938.92; //MeV
  double alpha=1/137.036; 
  const double pi=boost::math::constants::pi<double>();

  u3::U3CoefInit();
  int N=1;
  int Z=1;
  int A=N+Z;
  int N1B=1;
  int Nmax=6;
  double b2=1.0;  


  // Test 1: proton quadrupole operator
  if (true)
    {
      int J0=2;
      int g0=0;
      double coef=sqrt(5./(16*pi))*b2;
      
      u3shell::RelativeRMEsU3ST Operator;

      //Isoscalar term
      int T0=0;
      double coef_scalar=coef*(1-(Z-N)*1.0/A)/2;
      u3shell::Qintr(Nmax+2*N1B,Operator,A,T0,coef_scalar);
      
      //Isovector term
      T0=1;
      double coef_vector=coef/2;
      u3shell::Qintr(Nmax+2*N1B,Operator,A,T0,coef_vector);

      //Print out RMEs
      for(auto it=Operator.begin(); it!=Operator.end(); ++it)
        {
          u3shell::RelativeUnitTensorLabelsU3ST labels;
          int kappa0, L0;
          std::tie(labels,kappa0,L0)=it->first;
          double rme=it->second;
          // std::cout<<labels.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl;
          // std::cout<<rme<<std::endl;

          u3::SU3 x0;
          HalfInt S0, T0, S, T, Sp,Tp;
          int Np,N;
          std::tie(x0,S0,T0,Np,Sp,Tp,N,S,T)=labels.FlatKey();

          u3shell::RelativeStateLabelsU3ST bra(Np,Sp,Tp);
          u3shell::RelativeStateLabelsU3ST ket(N,S,T);
          u3shell::OperatorLabelsU3ST op_labels_conj(N-Np,u3::Conjugate(x0),S0,T0,g0);
          u3shell::RelativeUnitTensorLabelsU3ST unit_tensor_labels_conj(op_labels_conj,ket,bra);
          std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int> labels_conj(unit_tensor_labels_conj,kappa0,L0);
  
          int phase=ParitySign(u3::ConjugationGrade(x0)+S0);

          if(x0.lambda()==x0.mu())
            phase=phase*ParitySign(L0);

          double rme_conj=std::sqrt(u3::dim(u3::SU3(N,0)))/std::sqrt(u3::dim(u3::SU3(Np,0)))
                          *Hat(S)*Hat(T)/Hat(Sp)/Hat(Tp)
                          *phase
                          *Operator[labels_conj];

          
          // std::cout<<Operator[labels_conj]<<"  "<<ParitySign(u3::ConjugationGrade(x0)+S0+T0)<<std::endl;
          if(abs(rme-rme_conj)>1e-4)
          {
            std::cout<<"WRONG"<<std::endl;
            std::cout<<labels.Str()<<"  "<<kappa0<<"  "<<L0<<std::endl;
            std::cout<<rme<<std::endl;
            std::cout<<rme_conj<<std::endl;
          }
          // std::tuple<u3shell::RelativeUnitTensorLabelsU3ST,int,int>labels_conj()
          
        }

    }
//Test 2 interaction
  if(true)
    {
      int Jmax=4;
      int J0=0;
      int g0=0;
      int T0=-1;
      

      std::string daejeon_filename="/global/home/amccoy/scratch/data/spncci/interaction/rel/Daejeon16_Nmax40/Daejeon16_Nmax40_hw15.0_rel.dat";
      u3shell::CheckInteraction(A, Nmax, Jmax, J0, T0, g0, daejeon_filename);

      std::string coulomb_filename="/global/home/amccoy/scratch/data/spncci/interaction/rel/coulomb_Nmax20_steps500_rel.dat";
      u3shell::CheckInteraction(A, Nmax, Jmax, J0, T0, g0, coulomb_filename);

      std::string n3lo_filename="/global/home/amccoy/scratch/data/spncci/interaction/rel/N3LO/N3LO_srg2.15_Nmax30_hw20.0_rel.dat";
      u3shell::CheckInteraction(A, Nmax, Jmax, J0, T0, g0, n3lo_filename);

    }


  

}



