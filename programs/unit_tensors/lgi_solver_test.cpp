/****************************************************************
  lgi_solver_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  11/7/16 (aem): Created.
****************************************************************/
#include <fstream>
#include "cppformat/format.h"

#include "sp3rlib/u3coef.h"
#include "lgi/lgi_solver.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "u3shell/unit_tensor_expansion.h"

void NullSpaceCheck(
  const u3shell::SpaceU3SPN& space, 
  basis::MatrixVector& ncm_matrices)
{  
  Eigen::MatrixXd null;
  for(int i=0; i< space.size(); ++i)
  // for(auto matrix : ncm_matrices)
    {
      Eigen::MatrixXd& matrix=ncm_matrices[i];
      // std::cout<<matrix<<std::endl;
      u3shell::U3SPN labels(space.GetSubspace(i).GetSubspaceLabels());          
      // std::cout<<"Doing decomp..." << std::endl;
      int size=space.GetSubspace(i).size();
      std::cout<<labels.Str()<<"  size: "<<size<<std::endl;
      // std::cout<<matrix<<std::endl;
      if((matrix.cols()<2)&&(matrix.rows()<2))
        {
          if(fabs(matrix(0,0))<10e-6)
          {
            null=Eigen::MatrixXd::Identity(1,1);
            std::cout << "Rank : "<<0<<std::endl;
            std::cout << "CMF : "<<1<<std::endl<<std::endl;
          }
          else
          { 
            null=Eigen::MatrixXd::Zero(1,1);
            std::cout << "Rank : "<<1<<std::endl;
            std::cout << "CMF : "<<0<<std::endl<<std::endl;
          }
        }
      else
        {
          Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(matrix);
          lu_decomp.setThreshold(10e-6);
          std::cout << "Rank : "<<lu_decomp.rank()<<std::endl;
          std::cout << "null : "<<size-lu_decomp.rank()<<std::endl<<std::endl;
          //  std::cout<<"Doing kernel..." << std::endl;
          null=lu_decomp.kernel();
        }
      Eigen::MatrixXd zeros=matrix*null;
      double sum=0;
      int count=0;
      double max=0;
      for(int m=0; m<zeros.cols(); ++m)
        for(int n=0; n<zeros.rows(); ++n)
          {
            sum+=fabs(zeros(n,m));
            count++;
            if(fabs(zeros(n,m))>fabs(max))
              max=zeros(n,m);
          }
      std::cout<<"Average error: "<<sum/count<<"  Max error: "<<max<<std::endl;
      // std::cout<<"null space"<<std::endl<<null<<std::endl<<std::endl;
    }
}


int main(int argc, char **argv)
{
  u3::U3CoefInit();

  if(argc<4)
    std::cout<<"Syntax : A twice_Nsigma_0 <basis> <nrel> <brel>"<<std::endl;
  int A=std::stoi(argv[1]);
  int twice_Nsigma_0=std::stoi(argv[2]);
  HalfInt Nsigma_0=HalfInt(twice_Nsigma_0,2);
  std::string lsu3_filename = argv[3];
  std::string nrel_filename = argv[4];
  std::string brel_filename = argv[5];
  
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  
  std::ifstream is_brel(brel_filename.c_str());
  std::ifstream is_nrel(nrel_filename.c_str());

  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  std::cout<<"LSU3Shell basis  "<<std::endl;

  basis::MatrixVector ncm_matrices;
  std::cout<<"Ncm_matrices"<<std::endl;
  lsu3shell::GenerateNcmMatrixVector(A,is_nrel,basis_table,space,ncm_matrices);
  NullSpaceCheck(space,ncm_matrices);

  std::cout<<"Brel_matrices"<<std::endl;
  basis::MatrixVector brel_matrices;
  u3shell::OperatorLabelsU3ST brel_labels(-2,u3::SU3(0,2),0,0,0);
  u3shell::SectorsU3SPN brel_sectors(space,brel_labels,true);
  lsu3shell::ReadLSU3ShellRMEs(is_brel,brel_labels,basis_table,space,brel_sectors,brel_matrices);
  NullSpaceCheck(space,brel_matrices);


  // basis::MatrixVector matrix_vector(space.size());
  basis::MatrixVector lgi_expansion_matrix_vector(space.size());
  // lgi::GenerateNcmMatrixVector(Nsigma_0,nrel_filename, basis_table, space, matrix_vector);
  std::cout<<"Generating LGI expansion"<<std::endl;
  lgi::GenerateLGIExpansion(A,basis_table,space, is_brel,
  	is_nrel,lgi_expansion_matrix_vector);
  // for(auto matrix: lgi_expansion_matrix_vector)
  // 	std::cout<<matrix<<std::endl<<std::endl;

  std::ofstream os("LGI_file.dat");
  lgi::LGIVector lgi_vector;
  lgi::GetLGILabels(Nsigma_0,space,lgi_expansion_matrix_vector, lgi_vector);
  lgi::WriteLGILabels(lgi_vector,os);
  os.close();

  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(2, relative_unit_tensor_labels);
  for(int i=0; i<relative_unit_tensor_labels.size(); ++i)
    {
      u3shell::OperatorLabelsU3ST operator_labels(relative_unit_tensor_labels[i]);
      u3shell::SectorsU3SPN sectors(space,operator_labels,false);
      std::ifstream is_operator(fmt::format("relative_unit_{:06d}",i));

      if(not is_operator)
        {
          std::cout<<fmt::format("relative_unit_{:06d} not found",i)<<std::endl;
          continue;
        }
      basis::MatrixVector lsu3shell_operator_matrices;
      basis::MatrixVector spncci_operator_matrices;
      lsu3shell::ReadLSU3ShellRMEs(is_operator,operator_labels, basis_table,space, sectors,lsu3shell_operator_matrices);
      lgi::TransformOperatorToSpBasis(sectors,lgi_expansion_matrix_vector,lsu3shell_operator_matrices,spncci_operator_matrices);
      for(auto matrix : spncci_operator_matrices)
        std::cout<<matrix<<std::endl<<std::endl;
    }

}