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
  u3::U3CoefInit();

  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);

  assert(Nmax>=(N+Z));

  bool un_u3_restrict=false;
  if(((N==2)&&(Z==0))||((N==0)&&(Z==2)))
    un_u3_restrict=true;
  
  u3shell::RelativeUnitTensorCoefficientsU3ST identity_operator;
  u3shell::IdentityRelativeUnitTensorExpansion(0, Nmax, identity_operator, N+Z);

  std::string identity_file=fmt::format("Identity_{:02d}_Nmax{:02d}.recoupler",N+Z,Nmax);
  lsu3shell::GenerateLSU3ShellOperator(Nmax, identity_operator, identity_file,un_u3_restrict);

  //Generate Nintr operator up to Nmax cutoff
  std::string nintr_file=fmt::format("Nintr_{:02d}_Nmax{:02d}.recoupler",N+Z,Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
  u3shell::NintrRelativeUnitTensorExpansion(0,Nmax, Nrel_operator,N+Z);
  // for(auto it=Nrel_operator.begin(); it!=Nrel_operator.end(); ++it)
  // 	{
  // 		std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
  // 	}
  lsu3shell::GenerateLSU3ShellOperator(Nmax, Nrel_operator, nintr_file, un_u3_restrict);
}