/****************************************************************
  unit_tensor.h

  Unit tensor recursive evaluation
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/16/18 (aem): Created based on unit_tensor.h
****************************************************************/

#ifndef SPNCCI_SPNCCI_UNIT_TENSOR_H_
#define SPNCCI_SPNCCI_UNIT_TENSOR_H_
#include "spncci/spncci_basis.h"
#include "spncci/vcs_cache.h"
#include "u3shell/tensor_labels.h"


namespace spncci
{

  void 
  PopulateHypersectorsWithSeeds(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors_Nn0,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
    const std::vector<int>& rho0_values,
    basis::MatrixVector& unit_tensor_seed_blocks,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  );
  // Adds seeds to hyperblocks for both Nn=0 sectors and all other sectors


  void 
  AddNn0BlocksToHyperblocks(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors_Nn0,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  );

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
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
    );
} //namespace 

#endif
