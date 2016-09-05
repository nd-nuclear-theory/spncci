/****************************************************************
  lgi_unit_tensor_test.cpp

  Unit tensor algorithms
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
****************************************************************/
#include "cppformat/format.h"

#include "spncci/lgi_unit_tensor.h"
#include "lsu3shell_io/lsu3shell_interface.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();

	std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,relative_unit_tensor_labels);

  lsu3shell::GenerateLSU3ShellOperators(Nmax, relative_unit_tensor_labels);


  std::string filename_lgi="test_lgi_expansion.dat";
  Eigen::MatrixXd bra;
  spncci::ReadLGISU3Expansion(filename, bra, "bra");

  std::string filename_lgi="test_lgi_expansion";
  Eigen::MatrixXd ket;
  spncci::ReadLGISU3Expansion(filename, ket, "ket");
  std::cout<<bra<<"  "<<ket<<std::endl; 

//  std::string filename_rme="test_lsu3shell_rme.dat"

}


