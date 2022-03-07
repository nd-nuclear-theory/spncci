/****************************************************************
  print_lsu3shell_basis.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  11/2/21 (aem): Created.
****************************************************************/
#include "u3ncsm/dimensions.h"


int main(int argc, char **argv)
{

  if(argc<5)
  {
    std::cout<<"syntax: Z N Nmax Nstep"<<std::endl;
    std::exit(EXIT_FAILURE);
  }
 // nuclide
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);

  // Basis parameters
  int Nmax=std::stoi(argv[3]);
  int Nstep=std::stoi(argv[4]);

  lsu3shell::print_lsu3shell_basis_info({Z,N},Nmax,Nstep);

}//end main
