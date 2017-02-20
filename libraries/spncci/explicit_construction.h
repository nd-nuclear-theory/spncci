/****************************************************************
  explicit_construction.h

  Explicit construction of SpNCCI states in lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).

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
      const basis::MatrixVector& lgi_expansions,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const u3shell::SectorsU3SPN& Arel_sectors,
      const basis::MatrixVector& Arel_matrices,
      basis::MatrixVector& spncci_expansions
    );
  // Generate expansions of SpNCCI basis states in terms of lsu3shell
  // basis states, broken up by "baby SpNCCI subspaces", i.e., U3
  // subspaces segregated by irrep family.
  //
  // Limitation: Presently only supports Nex=0 or 2.
  //
  // Arguments:
  //   ...

}  // namespace

#endif
