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

extern double zero_threshold;

int main(int argc, char **argv)
{
	u3::U3CoefInit();

  // setup for test case
  int Nsigma_0= 3; // 11 for 6Li, 3 for 2H

  // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  std::string lsu3_filename("../../data/lsu3shell/lsu3shell_basis_2H_Nmax02.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  // Filling out basis_table, basis_provenance and space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  std::cout<<"Read Basis complete"<<std::endl;

 // Operator information
  std::string nrel_filename="../../data/lsu3shell/lsu3shell_rme_2H_Nrel_Nmax02.rme";
  std::string brel_filename="../../data/lsu3shell/lsu3shell_rme_2H_Brel_Nmax02.rme";
	basis::MatrixVector lgi_expansion_matrix_vector(space.size());
  basis::MatrixVector nrel_matrix_vector;
  std::ifstream is_nrel(nrel_filename.c_str());
  std::ifstream is_brel(brel_filename.c_str());
  assert(is_nrel.is_open());
  assert(is_brel.is_open());
  // lgi::GenerateNcmMatrixVector(Nsigma_0, is_nrel,basis_table,space, nrel_matrix_vector);
  
  u3shell::OperatorLabelsU3S brel_labels(-2,u3::SU3(0,2),0,0);
  //generate sectors for brel.
  u3shell::SectorsU3SPN brel_sectors(space,brel_labels,true);
  basis::MatrixVector brel_matrix_vector(space.size());
  lsu3shell::ReadLSU3ShellRMEs(is_brel,brel_labels,basis_table,space, brel_sectors,brel_matrix_vector);

  for(int i=0; i<brel_matrix_vector.size(); ++i)
  {
      std::cout<<brel_sectors.GetSector(i).ket_subspace().labels().Str()<<std::endl;
      std::cout<<"Matrix"<<std::endl;
      std::cout<<brel_matrix_vector[i]<<std::endl;
      Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(brel_matrix_vector[i]);
      lu_decomp.setThreshold(1e-6);
      std::cout << "Rank: "<<lu_decomp.rank()<<std::endl;
      std::cout<<"Doing kernel..." << std::endl;
      Eigen::MatrixXd null=lu_decomp.kernel();
      std::cout<<"null space"<<std::endl<<null<<std::endl<<std::endl;
  }

  u3shell::OperatorLabelsU3S nrel_labels(0,u3::SU3(0,0),0,0);
  //generate sectors for brel.
  u3shell::SectorsU3SPN nrel_sectors(space,nrel_labels,true);
  nrel_matrix_vector.resize(space.size());
  lsu3shell::ReadLSU3ShellRMEs(is_nrel,nrel_labels,basis_table,space, nrel_sectors,
    nrel_matrix_vector);

  for(int i=0; i<nrel_matrix_vector.size(); ++i)
  {
      std::cout<<nrel_sectors.GetSector(i).ket_subspace().labels().Str()<<std::endl;
      std::cout<<"Matrix"<<std::endl;
      std::cout<<nrel_matrix_vector[i]<<std::endl;
      if(nrel_matrix_vector[i].cols()<2)
        continue;
      Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(nrel_matrix_vector[i]);
      lu_decomp.setThreshold(1e-6);
      std::cout << "Rank: "<<lu_decomp.rank()<<std::endl;
      std::cout<<"Doing kernel..." << std::endl;
      Eigen::MatrixXd null=lu_decomp.kernel();
      std::cout<<"null space"<<std::endl<<null<<std::endl<<std::endl;
  }



  // basis::MatrixVector BrelNcm_vector(space.size()); 
  // lgi::GenerateBrelNcmMatrices(Nsigma_0,is_brel,is_nrel,basis_table,space, BrelNcm_vector);


  // lgi::GenerateLGIExpansion(Nsigma_0,basis_table,space,is_brel,
  // 													is_nrel,lgi_expansion_matrix_vector);
  // std::cout<<"expansions"<<std::endl;
  // for(auto matrix : lgi_expansion_matrix_vector)
  //   std::cout<<matrix<<std::endl<<std::endl;

  // is_brel.close();
  // is_nrel.close();

 //   // reading in operator rme's obtained form SU3RME
 //  basis::MatrixVector matrices;
 //  std::ifstream is(op_filename.c_str());
 //  if(!is)
 //    std::cout<<"file didn't open"<<std::endl;
 //  lsu3shell::ReadLSU3ShellRMEs(is,op_labels,basis_table,space, sectors,matrices);
 //  is.close();

 //  for(int i=0; i<matrices.size(); ++i)
 //  {
 //    // if(fabs(matrices[i].sum())>zero_threshold)
 //      // std::cout<<matrices[i]<<std::endl;
 //      std::cout<<"eigenvalues"<<std::endl;

 //      Eigen::VectorXcd eigs=matrices[i].eigenvalues();
 //      for(int j=0; j<eigs.size(); ++j)  
 //        std::cout<<eigs(j).real()<<std::endl;
 //  }
 //  basis::MatrixVector ncm_matrices;
 //  lgi::GenerateNcmMatrixVector(Nsigma_0,op_filename,basis_table,space,ncm_matrices);
 //  for(int i=0; i<ncm_matrices.size(); ++i)
 //  {
 //      if (matrices[i].rows()<2)
 //      continue;
 //    // if(fabs(matrices[i].sum())>zero_threshold)
 //      std::cout<<fmt::format("Matrix {}",i) << std::endl;
 //      std::cout<<matrices[i]<<std::endl;
 //      std::cout<<"eigenvalues"<<std::endl;

 //      Eigen::VectorXcd eigs=ncm_matrices[i].eigenvalues();
 //      for(int j=0; j<eigs.size(); ++j)  
 //        std::cout<<eigs(j).real()<<std::endl;

 //      std::cout<<"Doing decomp..." << std::endl;
 //      Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(ncm_matrices[i]);
 //      lu_decomp.setThreshold(1e-6);
 //      std::cout << "Rank: "<<lu_decomp.rank()<<std::endl;
 //      std::cout<<"Doing kernel..." << std::endl;
 //      Eigen::MatrixXd null=lu_decomp.kernel();
 //      std::cout<<"null space"<<std::endl<<null<<std::endl<<std::endl;
 //  }

}
