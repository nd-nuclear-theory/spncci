/****************************************************************
transform_seed_blocks.cpp

Transforms seed blocks by finding "Hamiltonian preferred" basis 
for Sp(3,R) subspace defined by all irreps with same Sp(3,R)SpSnS
labels.  

Transformations are obtained from eigenvectors of Hamiltonian matrix
in small spaces

7/19/18 (aem): Created. 
11/29/19 (aem): Overhall
****************************************************************/
#include <fstream>  

#include "fmt/format.h"
#include "lgi/lgi_unit_tensors.h"
#include "spncci/io_control.h"
#include "spncci/transform_basis.h"
#include "spncci/computation_control.h"

namespace spncci
{


void GetReorderingMatrix(const std::vector<int>& top_rows,
  const std::vector<int>& bottom_rows,
  spncci::OperatorBlock& reordering_matrix)
{
  int num_rows=top_rows.size()+bottom_rows.size();
  reordering_matrix=spncci::OperatorBlock::Zero(num_rows,num_rows);
  int i=0;
  for(int j : top_rows)
  {
    reordering_matrix(i,j)=1;
    ++i;
  }

  for(int j : bottom_rows)
    {
      reordering_matrix(i,j)=1;
      ++i;
    }
}


void ReorderBlock(
  spncci::OperatorBlock& block,
  spncci::OperatorBlock& reordering_matrix,
  int& num_above_threshold,
  double threshold
)
{
  int gamma_max=block.rows();
  auto norms=block.rowwise().squaredNorm();

  std::vector<int> above_threshold,below_threshold;

  // Determine which irreps are above the threshold 
  for(int g=0; g<gamma_max; ++g)
    {
      if(norms[g]>threshold)
        above_threshold.push_back(g);
      else
        below_threshold.push_back(g);
    }

  spncci::GetReorderingMatrix(above_threshold,below_threshold,reordering_matrix);
  
  //Reorder the block
  block=reordering_matrix*block;

  num_above_threshold=above_threshold.size();
}



 
}//namespace



int main(int argc, char **argv)
{

  if(argc<4)
  {
    std::cout<<"Syntax: <truncation_file_num> Nsmax Nmax "<<std::endl;
    exit(0);
  }

  //Opens file "spncci.dat" and reads in run parameters
  spncci::RunParameters run_parameters; 

  lgi::MultiplicityTaggedLGIVector lgi_families;
  spncci::SpNCCISpace spncci_space;
  spncci::BabySpNCCISpace baby_spncci_space;
  std::vector<spncci::SpaceSpBasis> spaces_spbasis;
  spncci::SigmaIrrepMap sigma_irrep_map; //Not used here
  spncci::SetUpSpNCCISpaces(
    run_parameters,lgi_families,
    spncci_space,sigma_irrep_map,
    baby_spncci_space,spaces_spbasis
    );

  // // Read in Jn sets for defining truncation 
  // // Read in eigenvectors from files
  // // Regroup eigenvectors 
 

  // spncci::RegroupIntoIrrepFamilies(
  //   spaces_spbasis,num_irrep_families,num_eigenvalues,
  //   eigenvectors,irrep_family_blocks
  // );

  // std::string test_filename=fmt::format("irrep_family_blocks_{}",hw);
  // spncci::WriteIrrepFamilyBlocks(
  //   run_parameters.J_values,  num_irrep_families,num_eigenvalues,
  //   lgi_full_space_index_lookup,irrep_family_blocks,test_filename
  // );



  // int truncation_file_num=std::stoi(argv[1]);
  // int Nsmax=std::stoi(argv[2]);
  // int Nmax=std::stoi(argv[3]);

  // std::array<int,2> nuclide;
  // std::vector<std::pair<int,int>> Jn_set;
  // int num_Jn;
  // double threshold;

  // // std::vector<std::pair<int,int>> Jn_set;
  // // Jn_set.emplace_back(1,0); // J,n where n is zero based
  // // Jn_set.emplace_back(3,0); // J,n where n is zero based

  // // double threshold=1e-6*Jn_set.size();

  // // std::pair<std::string,double> truncation_mode("None",threshold);
  // // std::pair<std::string,double> truncation_mode("Rank",threshold);
  // // std::pair<std::string,double> truncation_mode("Threshold",threshold);


  // std::ifstream is;
  // std::string truncation_file_name=fmt::format("truncation_{:03d}.dat",truncation_file_num);
  // is.open(truncation_file_name);

  // if(not is)
  // {
  //   std::cout<<"truncation file not found "<<std::endl;
  //   exit(0);
  // }

  // // Read in nuclide information
  // is>>nuclide[0]>>nuclide[1];

  // // Read in truncation information
  // std::string truncation_type;
  // is>>truncation_type>>threshold>>num_Jn;

  // if(truncation_type!="UI" && truncation_type!="IR" && truncation_type!="IF")
  //   { 
  //     std::cout<<truncation_type<<" "<< "is and invalid truncation type.  Possible truncation types are "<<std::endl
  //     <<"unitary_irrep : UI"<<std::endl<<"irrep : IR"<<std::endl<<"irrep_families : IF"<<std::endl;
  //     std::exit(EXIT_FAILURE);
  //   }

  // std::pair<std::string,double> truncation_mode(truncation_type,threshold);
  
  // std::cout<<truncation_type<<"  "<<threshold<<std::endl;
  // // Get list of J,n pairs which will contrubute to definition of transformation matrices 
  // std::string line;
  // int twice_J,n;
  // for(int j=0; j<num_Jn; j++)
  //   {
  //     is>>twice_J>>n;
  //     std::cout<<twice_J<<"  "<<n<<std::endl;
  //     Jn_set.emplace_back(twice_J,n);
  //   }
  // is.close();

  // std::cout<<"num Jn included "<<Jn_set.size()<<std::endl;
  // ///////////////////////////////////////////////////////////////////////////////////////////
  // ///////////////////////////////////////////////////////////////////////////////////////////


  // // Read in irrep family blocks 
  // std::cout<<"read in irrep family blocks"<<std::endl;
  // std::map<int,std::vector<spncci::OperatorBlocks>> irrep_family_blocks;
  // // std::vector<std::vector<spncci::OperatorBlocks>> irrep_family_blocks;
  // std::map<int,std::map<int,int>> J_index_lookup_table;
  // std::string irrep_family_blocks_file="irrep_family_blocks_20";
  // spncci::ReadIrrepFamilyBlocks(irrep_family_blocks,J_index_lookup_table,irrep_family_blocks_file);


  // // Get transformation matrices for each irrep family
  // std::cout<<"defined transformations "<<std::endl;
  // spncci::OperatorBlocks transformations;
  // // std::map<int,spncci::OperatorBlock> transformations;


  // spncci::DefineIrrepFamilyTransformations2(
  //   Jn_set,irrep_family_blocks,J_index_lookup_table,
  //   transformations,truncation_mode,threshold
  // );

  // // write truncated lgi families to file 
  // std::string truncated_lgi_filename
  //   =fmt::format("lgi_families_truncated_{:02d}_{:02d}_{:03d}.dat",Nsmax,Nmax,truncation_file_num);
  // spncci::WriteTruncatedLGIs(nuclide,transformations,truncated_lgi_filename);

  //  // write transformations to file
  // int Jn_file_num=Jn_set.size();
  // std::string transformations_file=fmt::format("transformations_{:02d}_{:02d}_{:03d}.dat",Nsmax,Nmax,truncation_file_num);
  // std::cout<<"write transformation matrices "<<std::endl;
  // spncci::WriteTransformationMatrices(transformations,transformations_file);

  // // //Read transformations from file 
  // // std::cout<<"Read transformation matrices "<<std::endl;
  // // spncci::OperatorBlocks transformations2;
  // // spncci::ReadTransformationMatrices(transformations_file,transformations2);

}