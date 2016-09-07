/****************************************************************
  generate_lgi_unit_tensors.cpp

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/29 (aem,mac): Created.
****************************************************************/
#include "cppformat/format.h"

#include "spncci/lgi_unit_tensor.h"
#include "lsu3shell/lsu3shell_operator.h"


int main(int argc, char **argv)
{
  u3::U3CoefInit();
  int Nmax=2;
	std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;

  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,relative_unit_tensor_labels);
  
  lsu3shell::GenerateLSU3ShellOperator(Nmax, relative_unit_tensor_labels);
  // for(int i=0; i<relative_unit_tensor_labels.size(); ++i)
  // {
  // 	std::cout<<fmt::format("operator{:6d}  {}", i,relative_unit_tensor_labels[i].Str())<<std::endl;
  // }
}
