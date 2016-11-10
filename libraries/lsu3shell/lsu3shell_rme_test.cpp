/****************************************************************
  lsu3shell_rme_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/5/16 (aem,mac): Created.
****************************************************************/
#include "lsu3shell/lsu3shell_rme.h"
#include <fstream>
#include <ostream>  
#include "cppformat/format.h"
#include "sp3rlib/u3coef.h"
#include "lgi/lgi_solver.h"
#include <eigen3/Eigen/LU>

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  std::string op_filename="../../data/lsu3shell/lsu3shell_rme_2H_Nmax02_Nrel.dat";
  u3shell::OperatorLabelsU3S op_labels(0,u3::SU3(0,0),0,0);
  bool scalar_op=true;
  // setup for test case
  int Nsigma_0=3;  // 11 for 6Li, 3 for 2H

  // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  std::string lsu3_filename("../../data/lsu3shell/lsu3shell_basis_2H_Nmax02.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  std::cout<<"Read Basis complete"<<std::endl;
  u3shell::SectorsU3SPN sectors(space,op_labels,scalar_op);
  basis::MatrixVector matrices;

  // reading in operator rme's obtained form SU3RME
  std::ifstream is(op_filename.c_str());
  if(!is)
    std::cout<<"file didn't open"<<std::endl;
  lsu3shell::ReadLSU3ShellRMEs(is,op_labels,basis_table,space, sectors,matrices);
  is.close();

  for(int i=0; i<matrices.size(); ++i)
  {
    // if(fabs(matrices[i].sum())>10e-13)
      // std::cout<<matrices[i]<<std::endl;
      std::cout<<"eigenvalues"<<std::endl;

      Eigen::VectorXcd eigs=matrices[i].eigenvalues();
      for(int j=0; j<eigs.size(); ++j)  
        std::cout<<eigs(j).real()<<std::endl;
  }
  basis::MatrixVector ncm_matrices;
  lgi::GenerateNcmMatrixVector(Nsigma_0,op_filename,basis_table,space,ncm_matrices);
  for(int i=0; i<ncm_matrices.size(); ++i)
  {
    // if(fabs(matrices[i].sum())>10e-13)
      // std::cout<<matrices[i]<<std::endl;
      std::cout<<"eigenvalues"<<std::endl;

      Eigen::VectorXcd eigs=ncm_matrices[i].eigenvalues();
      for(int j=0; j<eigs.size(); ++j)  
        std::cout<<eigs(j).real()<<std::endl;

      Eigen::MatrixXd null=ncm_matrices[i].fullPivLu().kernel();
      std::cout<<"null space"<<std::endl<<null<<std::endl<<std::endl;
  }

  // test matrix comparison
  //
  // // This is just a self-comparison, with one entry skewed...
  // basis::MatrixVector matrices_mod = matrices;
  // matrices_mod[0](0,0) += 3.14159;
  // double epsilon = 1e-8;
  // CompareLSU3ShellRMEs(
  //     std::cout,
  //     basis_provenance,
  //     space, 
  //     sectors,
  //     matrices,
  //     matrices_mod,
  //     epsilon,
  //     true  // verbose
  //   );


}
