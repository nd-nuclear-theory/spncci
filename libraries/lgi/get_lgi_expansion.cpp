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
	std::cout<<"Need overhall"<<std::endl;
  u3::U3CoefInit();
  int A=std::stoi(argv[1]);
  int Nsigma_0=std::stoi(argv[2]);
  std::string nrel_filename=argv[3];
  std::string brel_filename=argv[4];
  std::string basis_filename=argv[5];

  if(argc<5)
    std::cout<<"Syntax : A  Nsigma_0  <nrel>  <brel>  <basis>"<<std::endl;  

  // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  // e.g. "../../data/lsu3shell/lsu3shell_basis_2H_Nmax02.dat"
  std::string lsu3_filename(basis_filename.c_str());
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  // Filling out basis_table, basis_provenance and space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  std::cout<<"Read Basis complete"<<std::endl;

 // Operator information
  // std::string nrel_filename="../../data/lsu3shell/lsu3shell_rme_2H_Nrel_Nmax02.rme";
  // std::string brel_filename="../../data/lsu3shell/lsu3shell_rme_2H_Brel_Nmax02.rme";
  std::ifstream is_nrel(nrel_filename.c_str());
  std::ifstream is_brel(brel_filename.c_str());
  assert(is_nrel.is_open());
  assert(is_brel.is_open());

  basis::MatrixVector ncm_matrix_vector;
  lsu3shell::GenerateNcmMatrixVector(A, is_nrel,basis_table,space, ncm_matrix_vector);
  
  u3shell::OperatorLabelsU3ST brel_labels(-2,u3::SU3(0,2),0,0,0);
  //generate sectors for brel.
  u3shell::SectorsU3SPN brel_sectors(space,brel_labels,true);
  basis::MatrixVector brel_matrix_vector(space.size());
  // lsu3shell::ReadLSU3ShellRMEs(is_brel,brel_labels,basis_table,space, brel_sectors,brel_matrix_vector);

  basis::MatrixVector lgi_expansion_matrix_vector(space.size());

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

  // u3shell::OperatorLabelsU3ST nrel_labels(0,u3::SU3(0,0),0,0,0);
  // //generate sectors for brel.
  // u3shell::SectorsU3SPN nrel_sectors(space,nrel_labels,true);
  // basis::MatrixVector nrel_matrix_vector(space.size());
  // lsu3shell::ReadLSU3ShellRMEs(is_nrel,nrel_labels,basis_table,space, nrel_sectors,
  //   nrel_matrix_vector);

  // for(int i=0; i<nrel_matrix_vector.size(); ++i)
  // {
  //     std::cout<<nrel_sectors.GetSector(i).ket_subspace().labels().Str()<<std::endl;
  //     std::cout<<"Matrix"<<std::endl;
  //     std::cout<<nrel_matrix_vector[i]<<std::endl;
  //     if(nrel_matrix_vector[i].cols()<2)
  //       continue;
  //     Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(nrel_matrix_vector[i]);
  //     lu_decomp.setThreshold(1e-6);
  //     std::cout << "Rank: "<<lu_decomp.rank()<<std::endl;
  //     std::cout<<"Doing kernel..." << std::endl;
  //     Eigen::MatrixXd null=lu_decomp.kernel();
  //     std::cout<<"null space"<<std::endl<<null<<std::endl<<std::endl;
  // }



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
