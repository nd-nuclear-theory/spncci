/****************************************************************
  generate_nintr.cpp
  
	Generates kinetic energy for comparison with lsu3shell
	The kinetic energy is given by k^2 without any 2*mass
	factors.

  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT

  1/9/17 (aem): Created.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "fmt/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"


int main(int argc, char **argv)
{
  if(argc<5)
    std::cout<<"Syntax: Z  N  Nmax  N1B"<<std::endl;
  int max_lambda_plus_mu=39;
  u3::U3CoefInit(max_lambda_plus_mu);
  
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  int N1B=std::stoi(argv[4]);

  bool un_u3_restrict=false;
  if(((N==2)&&(Z==0))||((N==0)&&(Z==2)))
    un_u3_restrict=true;
  std::string trel_file_name_base=fmt::format("Krel2_{:02d}_Nmax{:02d}",N+Z,Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Trel_operator;
  u3shell::TrelRelativeUnitTensorExpansion(0,Nmax+2*N1B, Trel_operator,N+Z);
  // for(auto it=Trel_operator.begin(); it!=Trel_operator.end(); ++it)
  //   std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
  std::string trel_file_name=fmt::format("{}.recoupler",trel_file_name_base);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, Trel_operator, trel_file_name,un_u3_restrict);

}