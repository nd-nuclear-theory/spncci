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
  const std::vector<spncci::Matrix>& eigenvectors,
  int num_irrep_families,
  std::vector<spncci::OperatorBlocks>& irrep_family_blocks
);
// Regroup eigenvectors into blocks organzied by irrep family

void DefineIrrepFamilyTransformation(
  const spncci::SpaceSpJ& spj_space,
  std::vector<spncci::OperatorBlocks>& irrep_family_blocks
);


}

#endif