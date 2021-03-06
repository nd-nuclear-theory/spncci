/****************************************************************
  explicit_construction.h

  Explicit construction of SpNCCI states in lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).
  2/21/17 (aem): Add explict construction of full spncci basis
  2/22/17 (aem): Add computation of unit tensor rmes in explicit
    basis
  4/20/17 (aem): Updated to use hyperblock structure.
  Validated using programs/unit_tensors/explicit for Nmax=4
  10/4/17 (aem): Adapted for all A (including A<6)
****************************************************************/

#ifndef SPNCCI_SPNCCI_EXPLICIT_CONSTRUCTION_H_
#define SPNCCI_SPNCCI_EXPLICIT_CONSTRUCTION_H_

#include "spncci/spncci_basis.h"
#include "spncci/unit_tensor.h"

namespace spncci
{

  void
  ConstructSpNCCIBasisExplicit(
      const u3shell::SpaceU3SPN& lsu3shell_space,
      const spncci::SpNCCISpace& sp_irrep_families,
      const basis::MatrixVector& lgi_expansions,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const spncci::KMatrixCache& kinv_matrix_cache,
      const u3shell::SectorsU3SPN& Arel_sectors,
      const basis::MatrixVector& Arel_matrices,
      basis::MatrixVector& spncci_expansions,
      bool restrict_sp3r_u3_branching=false
    );

  // void 
  // ComputeUnitTensorSectorsExplicit(
  //   const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor,
  //   const u3shell::SpaceU3SPN& lsu3shell_space,
  //   const u3shell::SectorsU3SPN& lsu3shell_operator_sectors,
  //   basis::MatrixVector& lsu3shell_operator_matrices,
  //   const spncci::BabySpNCCISpace& baby_spncci_space,
  //   const basis::MatrixVector& spncci_expansions,
  //   spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_sectors_explicit
  //   );


  void 
  ComputeUnitTensorSectorsExplicit(
    const u3::U3& sigmap, const u3::U3& sigma,
    const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const u3shell::SpaceU3SPN& lsu3shell_space,
    const u3shell::SectorsU3SPN& lsu3shell_operator_sectors,
    basis::MatrixVector& lsu3shell_operator_matrices,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const basis::MatrixVector& spncci_expansions,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
    );

}  // namespace

#endif
