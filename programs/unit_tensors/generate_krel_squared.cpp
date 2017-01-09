/****************************************************************
  generate_nintr.cpp
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

	Generates kinetic energy for comparison with lsu3shell
	The kinetic energy is given by k^2 without any 2*mass
	factors.

  1/9/17 (aem): Created.
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

  bool un_u3_restrict=false;
  if(((N==2)&&(Z==0))||((N==0)&&(Z==2)))
    un_u3_restrict=true;
  std::string trel_file_name_base=fmt::format("Trel_{:02d}_Nmax{:02d}",A,Nmax);
  u3shell::RelativeUnitTensorCoefficientsU3ST Trel_operator;
  u3shell::TrelRelativeUnitTensorExpansion(0,Nmax+2*N1B, Trel_operator,A);

  std::string trel_file_name=fmt::format("{}.recoupler",trel_file_name_base);
  lsu3shell::GenerateLSU3ShellOperator(Nmax+2*N1B, Trel_operator, trel_file_name,un_u3_restrict);

}