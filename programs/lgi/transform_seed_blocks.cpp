/****************************************************************
transform_seed_blocks.cpp

Code which transforms seed blocks based on different linear combinations
of irreps in the same irrep family 

7/19/18 (aem): Created. 
****************************************************************/
#include <fstream>  

#include "cppformat/format.h"
#include "lgi/lgi_unit_tensors.h"
#include "spncci/transform_basis.h"


int main(int argc, char **argv)
{
  // lgi::binary_float_precision=8;
  std::vector<std::pair<int,int>> Jn_set;
  Jn_set.emplace_back(1,0); // J,n where n is zero based

  // Read in irrep family blocks 
  std::cout<<"read in irrep family blocks"<<std::endl;
  std::map<int,std::vector<spncci::OperatorBlocks>> irrep_family_blocks;
  std::map<int,std::map<int,int>> J_index_lookup_table;
  std::string irrep_family_blocks_file="temporary_test_file_20";
  spncci::ReadIrrepFamilyBlocks(irrep_family_blocks,J_index_lookup_table,irrep_family_blocks_file);


  // Get transformation matrices for each irrep family
  std::cout<<"defined transformations "<<std::endl;
  std::map<int,spncci::OperatorBlock> transformations;
  spncci::DefineIrrepFamilyTransformation(
    Jn_set,irrep_family_blocks,
    J_index_lookup_table,transformations
  );


  //Orthonormalize transformations....TODO

  // write transformations to file
  std::string transformations_file="transformations_out";
  std::cout<<"write transformation matrices "<<std::endl;
  spncci::WriteTransformationMatrices(transformations,transformations_file);

  // //Read transformations from file 
  // spncci::ReadTransformationMatrices(transformations_file,transformations);

  for(auto it_bra=transformations.begin(); it_bra!=transformations.end(); ++it_bra)
    for(auto it_ket=transformations.begin(); it_ket!=transformations.end(); ++it_ket)
      {
        int irrep_family_index_bra=it_bra->first;
        int irrep_family_index_ket=it_ket->first;

        std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensors;
        std::vector<int> rho0_values;


        std::string lgi_unit_tensor_filename
          =fmt::format("seeds/operators_{:06d}_{:06d}.dat",irrep_family_index_bra,irrep_family_index_bra);
        bool files_found=lgi::ReadUnitTensorLabels(lgi_unit_tensor_filename,lgi_unit_tensors,rho0_values);

        basis::OperatorBlocks<double> unit_tensor_seed_blocks;
        std::string seed_filename
          =fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",irrep_family_index_bra,irrep_family_index_ket);
        files_found&=lgi::ReadBlocks(seed_filename, lgi_unit_tensors.size(), unit_tensor_seed_blocks);

        if(not files_found)
          continue;

        spncci::OperatorBlock& bra_transformation=it_bra->second;
        spncci::OperatorBlock& ket_transformation=it_ket->second;

        basis::OperatorBlocks<double> transformed_seed_blocks(unit_tensor_seed_blocks.size());
        for(int b=0; b<transformed_seed_blocks.size(); ++b)
        {
            transformed_seed_blocks[b]
              =bra_transformation*unit_tensor_seed_blocks[b]*ket_transformation.transpose();
            // std::cout<<unit_tensor_seed_blocks[b]<<std::endl<<std::endl
            // <<bra_transformation<<std::endl<<std::endl
            // <<ket_transformation<<std::endl
            // <<std::endl<<transformed_seed_blocks[b]<<std::endl<<"***************************"<<std::endl<<std::endl;
        }



        lgi::WriteSeedsToFile(transformed_seed_blocks,irrep_family_index_bra, irrep_family_index_ket);

      }
}