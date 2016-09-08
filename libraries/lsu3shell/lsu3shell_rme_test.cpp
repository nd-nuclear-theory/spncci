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

// #include <Eigen/Eigenvalues>



int main(int argc, char **argv)
{
  u3::U3CoefInit();
  std::string op_filename="/Users/annamccoy/projects/spncci/data/lsu3shell_basis/Nrel2.dat";
  u3shell::OperatorLabelsU3S op_labels(0,u3::SU3(0,0),0,0);
  bool scalar_op=true;
  // setup for test case
  int Nsigma_0=3;  // 11 for 6Li, 3 for 2H

  // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  std::string lsu3_filename("/Users/annamccoy/projects/spncci/data/lsu3shell_basis/lsu3shell_basis_2H.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  std::cout<<"Read Basis complete"<<std::endl;
  u3shell::SectorsU3SPN sectors(space,op_labels,scalar_op);
  basis::MatrixVector matrix_vector;

  std::ifstream is(op_filename.c_str());
  if(!is)
    std::cout<<"file didn't open"<<std::endl;
  lsu3shell::ReadLSU3ShellRMEs(is,op_labels,basis_table,space, sectors,matrix_vector);
  is.close();

  for(int i=0; i<matrix_vector.size(); ++i)
  {
    // if(fabs(matrix_vector[i].sum())>10e-13)
      std::cout<<matrix_vector[i]<<std::endl;
      std::cout<<"  "<<std::endl;
      Eigen::VectorXcd eigs=matrix_vector[i].eigenvalues();
      std::cout<<eigs<<std::endl<<std::endl;

  }
}