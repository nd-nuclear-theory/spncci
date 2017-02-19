/****************************************************************
  computation_control.h

  High-level control code for computations in SpNCCI calculation
  programs.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).

****************************************************************/

#ifndef SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_
#define SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_

#include "lgi/lgi.h"
#include "spncci/spncci_basis.h"
#include "spncci/unit_tensor.h"
#include "sp3rlib/u3.h"
#include "u3shell/relative_operator.h"
#include "u3shell/u3spn_scheme.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // processing seed RMEs
  ////////////////////////////////////////////////////////////////

  void
  TransformSeedUnitTensorRMEs(
      const basis::MatrixVector& lgi_expansions,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::vector<u3shell::SectorsU3SPN>& lgi_unit_tensor_sectors,
      const std::vector<basis::MatrixVector>& lgi_unit_tensor_lsu3shell_matrices,
      std::vector<basis::MatrixVector>& lgi_unit_tensor_spncci_matrices
    );

  void
  StoreSeedUnitTensorRMEs(
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::vector<u3shell::SectorsU3SPN>& lgi_unit_tensor_sectors,
      const std::vector<basis::MatrixVector>& lgi_unit_tensor_spncci_matrices,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices,
      double zero_threshold
    );

}  // namespace

#endif
