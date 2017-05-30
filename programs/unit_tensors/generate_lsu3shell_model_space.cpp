/****************************************************************
  generate_lsu3shell_model_space.cpp

  Generates LSU3shell model space file for SU3RME and for basis lister.

  Command line parameters:
    Z : number of protons
    N : number of neutrons
    Nmax : max number of oscillator quant in relative basis
    Nstep : indicates if all or only same parity spaces are included
            in model space.  Can only take values of 1 or 2. 

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/29/17 (mac): Extracted from generate_lsu3shell_relative_tensors.cpp.

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

  // process arguments
  if(argc<4)
    {
      std::cout<<"Syntax: Protons Neutrons Nmax Nstep"<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  assert((Nstep==2)||(Nstep==1));
  int Nmin=Nmax%Nstep;

  // set up lsu3shell model space file for unit tensor calculations
  int parity=(Nstep==1)?-1:Nmin;
  // std::string model_space_filename = fmt::format("model_space_{:02d}_{:02d}_Nmax{:02d}.dat",Z,N,Nmax);
  std::string model_space_filename = "model_space.dat";
  lsu3shell::GenerateModelSpaceFile(model_space_filename,Z,N,Nmax,parity);

}
