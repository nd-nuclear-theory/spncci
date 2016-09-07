/****************************************************************
  check_lsu3shell_twobody_tensors.cpp

  Read in RMEs for two body tensors and check against expected
  values.

  Input files:
    basis.dat -- lsu3shell tabular basis listing file
    twobody_unit_*.rme -- output of SU3RME for each operator

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/7/16 (aem,mac): Created, starting frrom
    generate_lsu3shell_twobody_tensors.cpp.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_rme.h"
#include "spncci/sp_basis.h"
#include "spncci/lgi_unit_tensor.h"
#include "u3shell/unit_tensor_expansion.h"
#include "u3shell/two_body_operator.h"

void ReadBasis(u3shell::SpaceU3SPN space

// Given a nucleus (protons, neutrons)
// Given Nsigma_0
// Given a Nmax truncation
// Generate a list of twobody unit tensors
int main(int argc, char **argv)
{
  u3::U3CoefInit();

  // process arguments
  if(argc<4)
    {
      std::cout<<"Syntax: Protons Neutrons Nmin Nmax "<<std::endl;
      std::exit(1);
    }
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // will be either 1 or 2; 
  int Nstep=std::stoi(argv[4]);
  assert(Nstep<=2);
  int Nmin=Nmax%Nstep;

  // generate list of relative unit tensors up to Nmax cutoff
  std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST> twobody_unit_tensor_labels;
  u3shell::GenerateTwoBodyUnitTensorLabelsU3ST(Nmax, twobody_unit_tensor_labels);

  
  int num_unit=twobody_unit_tensor_labels.size();
  for(int i=0; i<num_unit; ++i)
    {
    }

}
