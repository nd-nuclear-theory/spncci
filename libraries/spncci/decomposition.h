/****************************************************************
  decomposition.h

  Code to generate eigenfunction decompositions.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/22/17 (mac): Created.
****************************************************************/

#ifndef SPNCCI_SPNCCI_DECOMPOSITION_H_
#define SPNCCI_SPNCCI_DECOMPOSITION_H_

#include "spncci/branching.h"

namespace spncci
{

  void CalculateNexDecompositions(
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::MatrixType>& eigenvectors,
      std::vector<basis::MatrixVector> Nex_decompositions
    );

  void CalculateBabySpNCCIDecompositions(
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::MatrixType>& eigenvectors,
      std::vector<basis::MatrixVector> Nex_decompositions
    );


}  // namespace

#endif
