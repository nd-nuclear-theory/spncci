/****************************************************************
  lgi_solver.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include <fstream>
#include <iostream>
#include <eigen3/Eigen/LU>

#include "cppformat/format.h"
#include "sp3rlib/u3coef.h"

#include "lsu3shell/lsu3shell_basis.h"
#include "lgi/lgi_solver.h"
#include "lsu3shell/lsu3shell_rme.h"



int main(int argc, char **argv)
{
	u3::U3CoefInit();

  // setup for test case
  int Nsigma_0=3;  // 11 for 6Li, 3 for 2H

  // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  std::string lsu3_filename("../../data/lsu3shell/lsu3shell_basis_2H_Nmax02.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  // Filling out basis_table, basis_provenance and space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  std::cout<<"Read Basis complete"<<std::endl;

 // Operator information
  std::string nrel_filename="../../data/lsu3shell/lsu3shell_rme_2H_Nmax02_Nrel.dat";
  std::string brel_filename="../../data/lsu3shell/lsu3shell_rme_2H_Nmax02_Brel.dat";
  // u3shell::OperatorLabelsU3S op_labels(0,u3::SU3(0,0),0,0);
  // bool scalar_op=true;
  // // Set up sector for operator
  // u3shell::SectorsU3SPN sectors(space,op_labels,scalar_op);
  
	basis::MatrixVector lgi_expansion_matrix_vector;
  lgi::GenerateLGIExpansion(Nsigma_0,basis_table,space,brel_filename, 
  													nrel_filename,lgi_expansion_matrix_vector);
}
