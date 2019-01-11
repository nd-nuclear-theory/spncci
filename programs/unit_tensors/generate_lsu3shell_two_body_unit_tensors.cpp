/****************************************************************
  generate_lsu3shell_two_body_unit_tensors.cpp

  Generate control and operator files for two-body unit tensors.

  Output files:
    two_body_unit_*.recoupler -- operators in input format for Tomas recoupler
    model_space.dat -- model space file
    two_body_operators.dat -- listing of operator file basenames

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/8/16 (aem,mac): Created.
  9/7/16 (mac): Remove header from operator list.
  9/10/16 (mac): Update filenames.
  12/2/16 (aem): Added bool argument for restricting on U(N)->U(3)
                branching
****************************************************************/

#include <fstream>
#include <ostream>  

#include "fmt/format.h"

#include "lsu3shell/lsu3shell_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"
#include "u3shell/two_body_operator.h"

int main(int argc, char **argv)
// Given a nucleus (protons, neutrons)
// Given Nsigma_0
// Given a Nmax truncation
// Generate a list of two_body unit tensors
{
  u3::U3CoefInit();

  // process arguments
  if(argc<7)
    {
      std::cout<<"Syntax: Protons Neutrons Nmax Nstep start stop "<<std::endl;
      std::exit(EXIT_FAILURE);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  int start=std::stoi(argv[5]);
  int stop=std::stoi(argv[6]);
  assert((1<=Nstep)&&(Nstep<=2));

  // set up unit tensor model space space
  int Nmin=Nmax%Nstep;
  // std::string model_space=fmt::format("model_space.{}_{}_Nmax{}",Z,N,Nmax);
  std::ofstream model_stream("model_space.dat");
  model_stream<<Z<<"  "<<N<<"  "<<-1<<std::endl;
  for(int Nex=Nmin; Nex<=Nmax; Nex+=2)
    model_stream<<Nex<<"  "<<-1<<std::endl;
  model_stream.close();

  // generate all relative unit tensors up to Nmax cutoff
  std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST> two_body_unit_tensor_labels;
  u3shell::GenerateTwoBodyUnitTensorLabelsU3ST(Nmax, two_body_unit_tensor_labels);
  int num_unit=two_body_unit_tensor_labels.size();
  std::cout<<"number of two body "<<num_unit<<std::endl;
  // for(int i=0; i<num_unit; ++i)  
  for(int i=start; i<std::min(stop,num_unit); ++i)
    {
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
      two_body_unit_tensor_coefficients[two_body_unit_tensor_labels[i]]=1;
      lsu3shell::GenerateLSU3ShellOperator(Nmax, two_body_unit_tensor_coefficients,i,true);
    }

  // write operator listing file
  std::ofstream control_stream("two_body_operators.dat");
  for(int i=start; i<std::min(stop,num_unit); ++i)
    control_stream<<fmt::format("two_body_unit_{:06d}",i)<<std::endl;
}
