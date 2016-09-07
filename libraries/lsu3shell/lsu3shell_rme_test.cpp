/****************************************************************
  lsu3shell_rme_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/5/16 (aem,mac): Created.
****************************************************************/
#include "lsu3shell/lsu3shell_rme.h"

#include <fstream>
#include <ostream>  
#include "sp3rlib/u3coef.h"
#include "cppformat/format.h"



int main(int argc, char **argv)
{
  u3::U3CoefInit();
  std::string op_filename="hi";
  u3shell::OperatorLabelsU3S op_labels(-2,u3::SU3(0,2),0,0);
  bool scalar_op=true;
  // setup for test case
  int Nsigma_0=11;  // 11 for 6Li

  // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  std::string lsu3_filename("lsu3basis_table.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);

  u3shell::SectorsU3SPN sectors(space,op_labels,scalar_op);
  basis::MatrixVector matrix_vector;

  std::ifstream is(op_filename.c_str());
  lsu3shell::ReadLSU3ShellRMEs(is,op_labels,basis_table,space, sectors,matrix_vector);

}