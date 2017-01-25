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
#include "lgi/null_solver.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "u3shell/unit_tensor_expansion.h"
#include "utilities/utilities.h"

int CountBasisIrreps(const u3shell::SpaceU3SPN& space)
{
  int count=0;
  for(int i=0; i<space.size(); ++i)
    {
      u3shell::SubspaceU3SPN subspace(space.GetSubspace(i));
      count+=subspace.size();
    }
  return count;
}

void NullSpaceCheck(
  const u3shell::SpaceU3SPN& space,
  u3shell::SectorsU3SPN& sectors, 
  basis::MatrixVector& operator_matrices,
  basis::MatrixVector& null_vectors_matrices)
{  
  int total_null=0;
  int total_rank=0;
  // double threshold=1e-6;
  // std::cout<<"Null solver threshold "<<threshold<<std::endl;
  for(int s=0; s<sectors.size(); ++s)
    {
      
      Eigen::MatrixXd& matrix=operator_matrices[s];
      int i=sectors.GetSector(s).ket_subspace_index();
      Eigen::MatrixXd& null_vectors=null_vectors_matrices[i];
      // std::cout<<"matrix"<<std::endl;
      // std::cout<<matrix<<std::endl;
      // std::cout<<"null vectors"<<std::endl;
      // std::cout<<null_vectors<<std::endl;
      u3shell::U3SPN labels(space.GetSubspace(i).GetSubspaceLabels());          
      // std::cout<<"Doing decomp..." << std::endl;
      int size=space.GetSubspace(i).size();
      // std::cout<<labels.Str()<<"subspace size: "<<size<<std::endl;
     
      // std::cout<<matrix<<std::endl;
      // bool verbose=false;
      // int dimension=matrix.rows();
      // lgi::FindNullSpaceSVD(matrix, null_vectors,threshold, verbose)
      int nullity=null_vectors.cols();
      int rank=size-nullity;
      total_null+=nullity;
      total_rank+=rank;
      std::cout << "rank : "<<rank<<std::endl;
      std::cout << "nullity : "<<nullity<<std::endl<<std::endl;
      // If there are no null vectors, continue to next sector
      if(nullity==0)
        continue;
      // Otherwise check is null_vectors are actually null vectors
      Eigen::MatrixXd zeros=matrix*null_vectors;
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
      std::cout<<"Average error : "<<sum/count<<"  Max error : "<<max<<std::endl;
      // std::cout<<"null space"<<std::endl<<null<<std::endl<<std::endl;
    }
  std::cout<<std::endl<<"Number of Sectors : "<<operator_matrices.size()<<"   Size of space : "<<space.size()<<std::endl;
  std::cout<<"Total rank : "<<total_rank<<"   Total nullilty : "<<total_null<<std::endl<<std::endl;
}


int main(int argc, char **argv)
{
  u3::U3CoefInit();
  double zero_threshold=1e-6;
  // parameters for null space solver
  double threshold=1e-4;
  bool verbose=false;


  if(argc<8)
    std::cout<<"Syntax : A  twice_Nsigma_0  Nmax  N1B  <basis>  <nrel>  <brel>"<<std::endl;
  int A=std::stoi(argv[1]);
  int twice_Nsigma_0=std::stoi(argv[2]);
  HalfInt Nsigma_0=HalfInt(twice_Nsigma_0,2);
  int Nmax=std::stoi(argv[3]);
  int N1B=std::stoi(argv[4]);
  std::string lsu3_filename = argv[5];
  std::string nrel_filename = argv[6];
  std::string brel_filename = argv[7];
  std::string arel_filename = argv[8];

  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;

  // Read in the basis
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  int num_states=CountBasisIrreps(space);
  std::cout<<"LSU3Shell basis  states : "<<num_states<<std::endl;
  
  // Read in Nrel and compute Ncm matrix sectors
  basis::MatrixVector ncm_matrices;
  std::cout<<"Ncm_matrices"<<std::endl;
  std::ifstream is_nrel(nrel_filename.c_str());
  lsu3shell::GenerateNcmMatrixVector(A,is_nrel,basis_table,space,ncm_matrices);
  is_nrel.close();
 
  u3shell::OperatorLabelsU3ST ncm_labels(0,u3::SU3(0,0),0,0,0);
  u3shell::SectorsU3SPN ncm_sectors(space,ncm_labels,true);

  // Checking the null vectors 
  basis::MatrixVector ncm_null_vectors;
  for(auto& ncm_matrix : ncm_matrices)
  {
    Eigen::MatrixXd null_vectors;
    lgi::FindNullSpaceSVD(ncm_matrix, null_vectors,threshold, verbose);
    ncm_null_vectors.push_back(null_vectors);
  }

  NullSpaceCheck(space,ncm_sectors,ncm_matrices,ncm_null_vectors);
  
  // Reading in Brel matrices 
  std::cout<<std::endl<<"Brel_matrices"<<std::endl;
  std::ifstream is_brel(brel_filename.c_str());
  u3shell::OperatorLabelsU3ST brel_labels(-2,u3::SU3(0,2),0,0,0);
  u3shell::SectorsU3SPN brel_sectors(space,brel_labels,true);
  basis::MatrixVector brel_matrices(brel_sectors.size());
  lsu3shell::ReadLSU3ShellRMEs(is_brel,brel_labels,basis_table,space,brel_sectors,brel_matrices);
  is_brel.close();

  // // Get Brel null vectors 
  // basis::MatrixVector brel_null_vectors;
  // for(auto& brel_matrix : brel_matrices)
  // {
  //   Eigen::MatrixXd null_vectors;
  //   lgi::FindNullSpaceSVD(brel_matrix, null_vectors,threshold, verbose);
  //   brel_null_vectors.push_back(null_vectors);
  // }
  // NullSpaceCheck(space, brel_matrices,brel_null_vectors);

  // Solving for lgi expansion from Brel+Ncm
  basis::MatrixVector lgi_expansion_matrix_vector;

  std::cout<<"Generating LGI expansion"<<std::endl;
  std::ifstream is_nrel2(nrel_filename.c_str());
  std::ifstream is_brel2(brel_filename.c_str());

  lgi::LGIVector lgi_vector;
  bool keep_empty_subspaces=true;
  lgi::GenerateLGIExpansion(A,Nsigma_0,basis_table,space, is_brel2,
  	is_nrel2,lgi_vector,lgi_expansion_matrix_vector, keep_empty_subspaces);
  
  is_nrel.close();
  is_brel.close();

  std::cout<<"Finished Generating expansion "<<std::endl;
  ZeroOutMatrix(lgi_expansion_matrix_vector,1e-5);
  for(int i=0; i< lgi_expansion_matrix_vector.size(); ++i)
  {
    auto& matrix=lgi_expansion_matrix_vector[i];
    std::cout<<lgi_vector[i].Str()<<std::endl;
  	std::cout<<matrix<<std::endl<<std::endl;
  }
  std::cout<<"finished printing "<<std::endl;

  // Writing lgi family labels to file
  std::ofstream os("LGI_file.dat");
  lgi::WriteLGILabels(lgi_vector,os);
  os.close();

  // Checking the null space of LGI
  NullSpaceCheck(space, brel_sectors, brel_matrices,lgi_expansion_matrix_vector);
  std::cout<<"checking that lgi's are cmf"<<std::endl;
  NullSpaceCheck(space,ncm_sectors,ncm_matrices,lgi_expansion_matrix_vector);
  

  std::cout<<std::endl<<"Arel_matrices"<<std::endl;
  std::ifstream is_arel(arel_filename.c_str());
  u3shell::OperatorLabelsU3ST arel_labels(2,u3::SU3(2,0),0,0,0);
  u3shell::SectorsU3SPN arel_sectors(space,arel_labels,true);
  basis::MatrixVector arel_matrices(arel_sectors.size());
  lsu3shell::ReadLSU3ShellRMEs(is_arel,arel_labels,basis_table,space,arel_sectors,arel_matrices);
  is_arel.close();

  for(int i=0; i<arel_sectors.size(); ++i)
  {
    auto& arel_sector=arel_sectors.GetSector(i);
    int lgi_index=arel_sector.ket_subspace_index();
    int w_index=arel_sector.bra_subspace_index();
    Eigen::MatrixXd& Arel=arel_matrices[i];
    Eigen::MatrixXd& LGIs=lgi_expansion_matrix_vector[lgi_index];
    // std::cout<<"Arel "<<std::endl;
    // std::cout<<Arel<<std::endl;
    // std::cout<<"LGI "<<std::endl;
    // std::cout<<LGIs<<std::endl;
    int j=ncm_sectors.LookUpSectorIndex(w_index,w_index,1);
    Eigen::MatrixXd& Ncm=ncm_matrices[j];
    // std::cout<<"Brel "<<std::endl;
    // std::cout<<Brel<<std::endl;
    Eigen::MatrixXd null_matrix=Ncm*Arel*LGIs;
    // std::cout<<"omegas "<<std::endl;
    // std::cout<<omega_matrix<<std::endl;
    // std::cout<<"Brel on omega_matrix "<<std::endl;
    basis::MatrixVector null_matrix_vector;
    null_matrix_vector.push_back(null_matrix);
    ZeroOutMatrix(null_matrix_vector,threshold);
    
    if(not CheckIfZeroMatrix(null_matrix_vector[0]))
      {
        std::cout<<null_matrix_vector[0]<<std::endl;
        std::cout<<"cmf matrix ?"<<std::endl;
      }
    // Check orthogonality of laddered states 
    // Eigen::MatrixXd w_matrix=Arel*LGIs;
    // Eigen::MatrixXd orthog_set=w_matrix.transpose()*w_matrix;
    // basis::MatrixVector orthog_vector;
    // orthog_vector.push_back(orthog_set);
    // ZeroOutMatrix(orthog_vector,threshold);
    // std::cout<<orthog_vector[0]<<std::endl<<std::endl;

  }


  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax+2*N1B, relative_unit_tensor_labels,-1,0,false);

  u3shell::OperatorLabelsU3ST id_labels(0,u3::SU3(0,0),0,0,0);
  u3shell::SectorsU3SPN id_sectors(space,id_labels,false);
  basis::MatrixVector id_lsu3shell(id_sectors.size()), id_spncci(id_sectors.size());
  basis::SetOperatorToZero(id_sectors,id_lsu3shell);
  basis::SetOperatorToZero(id_sectors,id_spncci);
  std::cout<<"number of relative unit tensors "<<relative_unit_tensor_labels.size()<<std::endl;
  for(int i=0; i<relative_unit_tensor_labels.size(); ++i)
    {
      // std::cout<<"tensor "<<i<<std::endl;
      u3shell::OperatorLabelsU3ST operator_labels(relative_unit_tensor_labels[i]);
      u3shell::SectorsU3SPN sectors(space,operator_labels,false);
      std::ifstream is_operator(fmt::format("relative_unit_{:06d}.rme",i));

      if(not is_operator)
        {
          std::cout<<fmt::format("relative_unit_{:06d}.rme not found",i)<<std::endl;
          continue;
        }
      basis::MatrixVector lsu3shell_operator_matrices;
      basis::MatrixVector spncci_operator_matrices;
      // std::cout<<"reading in"<<std::endl;
      lsu3shell::ReadLSU3ShellRMEs(is_operator,operator_labels, basis_table,space, sectors,lsu3shell_operator_matrices);
      // for(auto matrix : lsu3shell_operator_matrices)
      //   std::cout<<matrix<<std::endl<<std::endl;
      // std::cout<<"transforming"<<std::endl;
      lgi::TransformOperatorToSpBasis(sectors,lgi_expansion_matrix_vector,lsu3shell_operator_matrices,spncci_operator_matrices);
      // for(auto matrix : spncci_operator_matrices)
      //   std::cout<<matrix<<std::endl<<std::endl;

      // Checking identity before and after transformation 
      // std::cout<<"trans complete"<<std::endl;
      if(operator_labels==id_labels)
        {
          assert(sectors.size()==id_sectors.size());
          // std::cout<<" id sum"<<std::endl;
          for(int i=0; i<sectors.size(); ++i)
            {
              // std::cout<<id_lsu3shell[i]<<"  "<<lsu3shell_operator_matrices[i]<<std::endl;
              id_lsu3shell[i]+=lsu3shell_operator_matrices[i];
              // id_spncci[i]+=spncci_operator_matrices[i];
            }
        }
    }
    // std::cout<<"lsu3shell identity check"<<std::endl;
    // for(auto lsu3_matrix : id_lsu3shell)
    //   std::cout<<lsu3_matrix<<std::endl<<std::endl;
    std::cout<<"transforming identity "<<std::endl;
    lgi::TransformOperatorToSpBasis(id_sectors,lgi_expansion_matrix_vector,id_lsu3shell,id_spncci);

    ZeroOutMatrix(id_spncci,1e-5);

    std::cout<<"spncci"<<std::endl;
    for(auto sp_matrix : id_spncci)
      std::cout<<sp_matrix<<std::endl<<std::endl;


}