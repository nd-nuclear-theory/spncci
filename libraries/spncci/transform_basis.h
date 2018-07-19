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
	std::map<int,std::map<int,int>>& J_index_lookup_table,
  const std::string& filename
);


void  DefineIrrepFamilyTransformation(
  const std::vector<std::pair<int,int>>& Jn_set,
  const std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  std::map<int,std::map<int,int>>& J_index_lookup_table,
  std::map<int,spncci::OperatorBlock>& transformations
);


void WriteTransformationMatrices(  
	const std::map<int,spncci::OperatorBlock>& transformations,
  const std::string& filename
);

}

#endif