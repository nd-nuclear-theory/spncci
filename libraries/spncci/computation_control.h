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
// #include "u3shell/upcoupling.h"
#include "u3shell/u3spn_scheme.h"
#include "spncci/branching_u3s.h"
#include "spncci/branching_u3lsj.h"

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
      HalfInt Nsigma_max,
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
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices,
      bool verbose = false
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


  void
  ConstructObservableU3S(
      int Nmax, int N1v,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceU3S& space_u3s,
      const spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices,
      const u3shell::RelativeRMEsU3ST& relative_rmes,
      std::vector<spncci::SectorLabelsU3S>& sectors_u3s,
      basis::MatrixVector& matrices_u3s
    );
  //   Nmax, N1v (input): why are these needed???
  //   baby_spncci_space (input): baby SpNCCI space
  //   space_u3s (input): U3S regrouped SpNCCI space
  //   unit_tensor_matrices (input): unit tensor rmes in SpNCCI space
  //   relative_rmes (input): observable upcoupled relative rmes
  //   sectors_u3s (output): operator sectors [TODO upgrade to a proper "sectors" container]
  //   matrices_u3s (output): operator matrices


  void 
  ConstructBranchedObservables(
    const spncci::SpaceU3S& space_u3s,
    const std::vector<std::vector<spncci::SectorLabelsU3S>>& observable_sectors_u3s,
    const std::vector<basis::MatrixVector>& observable_matrices_u3s,
    std::map<HalfInt,spncci::SpaceLS>& spaces_lsj,
    int num_observables,
    const std::vector<HalfInt>& J_values,
    int J0,
    std::vector<std::map<HalfInt,Eigen::MatrixXd>>& observable_matrices
    );
  // Construct J branched observable matrices 


}  // namespace

#endif
