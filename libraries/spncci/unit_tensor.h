/****************************************************************
  unit_tensor.h

  Unit tensor recursive evaluation
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/15/16 (aem,mac): Created.
  4/14/16 (aem): Added Np,N sector look ups and changed iteration order
  4/20/16 (aem): Defined GenerateUnitTensorU3SectorLabels function and 
                 changed interation in UnitTensorMatrixGenerator
  12/7/16 (aem): Overhall of unit tensor rme calculation
  12/21/16 (aem): Factored out conjugate sector calculation
  2/2/17 (mac): Add typedef UnitTensorMatricesByIrrepFamily.
  2/20/17 (mac): Update argument types on GenerateUnitTensorMatrix.
  2/23/17 (aem): Add conjugate sector flag to GenerateUnitTensorU3SectorLabels.
  2/23/17 (aem): Removed redundent calculation of edge sectors in recurrence.
  5/5/17 (aem) : Update from souffle to baby spncci
  6/16/17 (aem) : Reimplement omp parallization
  7/1/17 (aem): Fix spin conjugation phase
  10/4/17 (aem): Adjusted recurrence to work for A<6
****************************************************************/

#ifndef SPNCCI_SPNCCI_UNIT_TENSOR_H_
#define SPNCCI_SPNCCI_UNIT_TENSOR_H_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash_fwd.hpp>
#include <eigen3/Eigen/Eigen>

#include "basis/operator.h"
#include "spncci/spncci_basis.h"
#include "spncci/vcs_cache.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/vcs.h"
#include "u3shell/tensor_labels.h"


namespace spncci
{
  // struct UnitTensorMatrixStatistics
  // // Structure for reporting diagnostic information on unit tensor matrices.
  // {
  //   UnitTensorMatrixStatistics()
  //   : num_irrep_family_index_pairs(0), num_Nn_pairs(0), num_unit_tensor_sectors(0), num_nonzero_unit_tensor_sectors(0), num_matrix_elements(0), num_nonzero_matrix_elements(0)
  //   {};

  //   int num_irrep_family_index_pairs;
  //   int num_Nn_pairs;
  //   int num_unit_tensor_sectors;
  //   int num_nonzero_unit_tensor_sectors;
  //   int num_matrix_elements;
  //   int num_nonzero_matrix_elements;
  // };

  // spncci::UnitTensorMatrixStatistics
  //   GenerateUnitTensorMatrixStatistics(
  //       const spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices
  //     );
  // // Count entries at each level of branching in unit tensor matrices container.
  // //
  // // Arguments:
  // //   unit_tensor_matrices (input): seed rmes (input) and computed rmes (output)
  // //
  // // Returns:
  // //   statistics



void
ComputeUnitTensorHyperblocks(
  int Nmax, int N1v,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::KMatrixCache& k_matrix_map,
  const spncci::KMatrixCache& kinv_matrix_map,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
  const std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
  bool restrict_sp3r_to_u3_branching,
  int number_of_threads=1

  );
} //namespace 

#endif
