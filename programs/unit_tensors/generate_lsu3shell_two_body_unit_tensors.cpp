/****************************************************************
  generate_lsu3shell_two_body_unit_tensors.cpp

  Generate control and operator files for two-body unit tensors.

  Output files:
    two_body_unit_*.recoupler -- operators in input format for Tomas recoupler
    model_space -- model space file
    two_body_operators.dat -- listing of operator file basenames

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/8/16 (aem,mac): Created.
  9/7/16 (mac): Remove header from operator list.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "spncci/sp_basis.h"
// #include "spncci/lgi_unit_tensor.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"
#include "u3shell/two_body_operator.h"

// Given a nucleus (protons, neutrons)
// Given Nsigma_0
// Given a Nmax truncation
// Generate a list of two_body unit tensors
int main(int argc, char **argv)
{
  u3::U3CoefInit();

  // process arguments
  if(argc<4)
    {
      std::cout<<"Syntax: Protons Neutrons Nmin Nmax "<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  assert(Nstep<=2);

  // set up unit tensor model space space
  int Nmin=Nmax%Nstep;
  // std::string model_space=fmt::format("model_space.{}_{}_Nmax{}",Z,N,Nmax);
  std::ofstream model_stream("model_space");
  model_stream<<Z<<"  "<<N<<"  "<<-1<<std::endl;
  for(int Nex=Nmin; Nex<=Nmax; Nex+=2)
    model_stream<<Nex<<"  "<<-1<<std::endl;
  model_stream.close();

  // generate all relative unit tensors up to Nmax cutoff
  std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST> two_body_unit_tensor_labels;
  u3shell::GenerateTwoBodyUnitTensorLabelsU3ST(Nmax, two_body_unit_tensor_labels);
  int num_unit=two_body_unit_tensor_labels.size();
  for(int i=0; i<num_unit; ++i)
    {
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
      two_body_unit_tensor_coefficients[two_body_unit_tensor_labels[i]]=1;
      lsu3shell::GenerateLSU3ShellOperator(Nmax, two_body_unit_tensor_coefficients,i);
    }

  // write operator listing file
  std::ofstream control_stream("two_body_operators.dat");
  for(int i=0; i<num_unit; ++i)
    control_stream<<fmt::format("two_body_unit_{:06d}",i)<<std::endl;
}
