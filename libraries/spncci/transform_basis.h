/****************************************************************
  transform_basis.h

  Code to generate spncci basis transformations
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  7/11/18 (aem): Created.
****************************************************************/

#ifndef SPNCCI_TRANSFORM_BASIS_H_
#define SPNCCI_TRANSFORM_BASIS_H_

#include <iostream>
#include <string>

#include "cppformat/format.h"
#include "spncci/branching.h"
#include "spncci/branching_u3s.h"
#include "spncci/parameters.h"

namespace spncci
{

void RegroupIntoIrrepFamilies(
    const spncci::SpaceSpJ& spj_space,
    int num_irrep_families,
    int num_eigenvalues,
    const std::vector<spncci::Matrix>& eigenvectors,
    std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks
  );
// Regroup eigenvectors into blocks organzied by irrep family
// Blocks are <blocks(dim,gamma_max)>

void WriteIrrepFamilyBlocks(  
  const spncci::SpaceSpJ& spj_space,
  int num_irrep_families,
  int num_eigenvalues,
  const std::vector<int>& lgi_full_space_index_lookup,
  const std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  const std::string& filename
);
// Write irrep family blocks to file
//
// In file:
//		<float precision, num_J_values, num_eigenvalues, num_irrep_families>
// 		
// 		For each irrep family,
// 			For each J value,
//				<2J, rows, cols>
//				For each n
//					<block(gamma_max,dim)>
//					

void ReadIrrepFamilyBlocks(
  std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
	// std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
	std::map<int,std::map<int,int>>& J_index_lookup_table,
  const std::string& filename
);


void  DefineIrrepFamilyTransformations(
  const std::vector<std::pair<int,int>>& Jn_set,
  std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  // const std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  std::map<int,std::map<int,int>>& J_index_lookup_table,
  spncci::OperatorBlocks& transformations,
  const std::pair<std::string,double>& truncation_mode
);


void WriteTransformationMatrices(
  const spncci::OperatorBlocks& transformations,
  const std::string& filename
);
// <float_precision>
// for each irrep_family 
// 		<irrep_family_index, gamma_max, transformation matrix>
// ...

void WriteTruncatedLGIs(
    // const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const std::array<int,2>& nuclide,
    // int Nsmax, int Nmax, int truncation_file_num,
    const spncci::OperatorBlocks& transformations,
    const std::string& truncated_lgi_filename
  );
// write truncated set of lgi familes to file
//    lgi_families_truncated_Nsigma,max_Nmax_index.dat


void ReadTransformationMatrices(  
	const std::string& filename,
  spncci::OperatorBlocks& transformations
);

void  RegroupBlocks(
  const std::vector<std::pair<int,int>>& Jn_set,
  const std::vector<spncci::OperatorBlocks>& blocks,
  std::map<int,int>& J_index_table,
  spncci::OperatorBlock& irrep_family_block
);

void TransformSeeds(
  int bra_index,int ket_index,
  spncci::OperatorBlocks& transformations,
  basis::OperatorBlocks<double>& unit_tensor_seed_blocks
  );



}

#endif