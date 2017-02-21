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

  void
  RecurseUnitTensors(
      int N1v, int Nmax,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::KMatrixCache k_matrix_cache,
      u3::UCoefCache& u_coef_cache,
      u3::PhiCoefCache& phi_coef_cache,
      const std::map<int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>> unit_tensor_labels,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices
    );
  // Recursively populate sectors for unit tensors.
  //
  // (Does recursing unit tensors implies unit tensors are already cursed?)
  //
  // FUTURE: Should be able to do recurrence without externally
  // imposing N1v, Nmax.  Rather, use truncation defined by the
  // spncci_space.  Requires corresponding change to
  // spncci::GenerateUnitTensorMatrix.
  //
  // Arguments:
  //    ...

}  // namespace

#endif
