/****************************************************************
transform_seed_blocks.cpp

Code which transforms seed blocks based on different linear combinations
of irreps in the same irrep family 

7/19/18 (aem): Created. 
****************************************************************/
#include <fstream>  

#include "cppformat/format.h"
#include "lgi/lgi_unit_tensors.h"
#include "spncci/io_control.h"
#include "spncci/transform_basis.h"


int main(int argc, char **argv)
{

  if(argc<2)
  {
    std::cout<<"Syntax: <truncation file num> "<<std::endl;
    exit(0);
  }

  int truncation_file_num=std::stoi(argv[1]);




  std::array<int,2> nuclide;
  std::vector<std::pair<int,int>> Jn_set;
  int Nmax,Nsmax, num_Jn;
  double threshold;

  // std::vector<std::pair<int,int>> Jn_set;
  // Jn_set.emplace_back(1,0); // J,n where n is zero based
  // Jn_set.emplace_back(3,0); // J,n where n is zero based

  // double threshold=1e-6*Jn_set.size();

  // std::pair<std::string,double> truncation_mode("None",threshold);
  // std::pair<std::string,double> truncation_mode("Rank",threshold);
  // std::pair<std::string,double> truncation_mode("Threshold",threshold);


  std::ifstream is;
  std::string truncation_file_name=fmt::format("truncation_{:03d}.dat",truncation_file_num);
  is.open(truncation_file_name);

  if(not is)
  {
    std::cout<<"truncation file not found "<<std::endl;
    exit(0);
  }

  // Read in nuclide information
  is>>nuclide[0]>>nuclide[1];

  // Read in truncation information
  std::string truncation_type;
  is>>truncation_type>>threshold>>Nsmax>>Nmax>>num_Jn;
  std::pair<std::string,double> truncation_mode(truncation_type,threshold);
  
  std::cout<<truncation_type<<"  "<<threshold<<std::endl;
  // Get list of J,n pairs which will contrubute to definition of transformation matrices 
  std::string line;
  int twice_J,n;
  for(int j=0; j<num_Jn; j++)
    {
      // std::istringstream line_stream(line);
      // line_stream>>twice_J>>n;
      is>>twice_J>>n;
      std::cout<<twice_J<<"  "<<n<<std::endl;
      Jn_set.emplace_back(HalfInt(twice_J,2),n);
    }
  is.close();

  std::cout<<"num Jn included "<<Jn_set.size()<<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////


  // Read in irrep family blocks 
  std::cout<<"read in irrep family blocks"<<std::endl;
  std::map<int,std::vector<spncci::OperatorBlocks>> irrep_family_blocks;
  // std::vector<std::vector<spncci::OperatorBlocks>> irrep_family_blocks;
  std::map<int,std::map<int,int>> J_index_lookup_table;
  std::string irrep_family_blocks_file="irrep_family_blocks_20";
  spncci::ReadIrrepFamilyBlocks(irrep_family_blocks,J_index_lookup_table,irrep_family_blocks_file);


  // Get transformation matrices for each irrep family
  std::cout<<"defined transformations "<<std::endl;
  spncci::OperatorBlocks transformations;
  // std::map<int,spncci::OperatorBlock> transformations;


  spncci::DefineIrrepFamilyTransformations(
    Jn_set,irrep_family_blocks,J_index_lookup_table,
    transformations,truncation_mode
  );

  // write truncated lgi families to file 
  std::string truncated_lgi_filename
    =fmt::format("lgi_families_truncated_{:02d}_{:02d}_{:03d}.dat",Nsmax,Nmax,truncation_file_num);
  spncci::WriteTruncatedLGIs(nuclide,transformations,truncated_lgi_filename);

   // write transformations to file
  int Jn_file_num=Jn_set.size();
  std::string transformations_file=fmt::format("transformations_{:02d}_{:02d}_{:03d}.dat",Nmax,Nsmax,truncation_file_num);
  std::cout<<"write transformation matrices "<<std::endl;
  spncci::WriteTransformationMatrices(transformations,transformations_file);

  //Read transformations from file 
  std::cout<<"Read transformation matrices "<<std::endl;
  spncci::OperatorBlocks transformations2;
  spncci::ReadTransformationMatrices(transformations_file,transformations2);

}