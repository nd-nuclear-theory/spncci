/****************************************************************
  generate_nintr.cpp
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  12/1/16 (aem): Created.
  12/2/16 (aem): Added bool for U(N)->U(3) restriction when A=2
****************************************************************/
#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"

int main(int argc, char **argv)
{
  if(argc<5)
    std::cout<<"Syntax: Z  N  Nmax  N1B"<<std::endl;

  u3::U3CoefInit();

  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  int N1B=std::stoi(argv[4]);
  int Nmin=Nmax%2;
  int A=N+Z;
  // assert(Nmax>=(N+Z));

  bool un_u3_restrict=false;
  if((N==0)||(Z==0))
    un_u3_restrict=true;
  
  lsu3shell::GenerateModelSpaceFile(Z, N, Nmax, Nmax%2);

  u3shell::RelativeUnitTensorCoefficientsU3ST identity_operator;
  u3shell::IdentityRelativeUnitTensorExpansion(0, Nmax+2*N1B, identity_operator, N+Z);
  std::string identity_file=fmt::format("Identity_{:02d}_Nmax{:02d}.recoupler",N+Z,Nmax);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, identity_operator, identity_file,un_u3_restrict);

  //Generate Nrel operator up to Nmax cutoff
  std::string nrel_file=fmt::format("Nrel_{:02d}_Nmax{:02d}.recoupler",N+Z,Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
  u3shell::NintrRelativeUnitTensorExpansion(0,Nmax+2*N1B, Nrel_operator,N+Z);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, Nrel_operator, nrel_file, un_u3_restrict);

  //Generate Nintr operator up to Nmax cutoff
  std::string brel_file=fmt::format("Brel_{:02d}_Nmax{:02d}.recoupler",N+Z,Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Brel_operator;
  u3shell::BrelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B, Brel_operator,A);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, Brel_operator, brel_file, un_u3_restrict);

  //Generate Nintr operator up to Nmax cutoff
  std::string arel_file=fmt::format("Arel_{:02d}_Nmax{:02d}.recoupler",N+Z,Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Arel_operator;
  u3shell::ArelRelativeUnitTensorExpansion(Nmin,Nmax+2*N1B, Arel_operator,A);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, Arel_operator, arel_file, un_u3_restrict);


  // u3shell::BrelRelativeUnitTensorExpansion(0,Nmax+2*N1B, Brel_operator, N+Z);
  // lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, Brel_operator, brel_file, un_u3_restrict);
}