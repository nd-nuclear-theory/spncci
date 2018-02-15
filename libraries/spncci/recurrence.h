/****************************************************************
  unit_tensor.h

  Unit tensor recursive evaluation
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/16/18 (aem): Created based on unit_tensor.h
  1/30/18 (aem): Extracted all seed and recurrence functions from
    computational control 
****************************************************************/

#ifndef SPNCCI_SPNCCI_UNIT_TENSOR_H_
#define SPNCCI_SPNCCI_UNIT_TENSOR_H_
#include "spncci/spncci_basis.h"
#include "spncci/vcs_cache.h"
#include "u3shell/tensor_labels.h"


namespace spncci
{
  typedef std::pair<int,int> LGIPair;
  typedef std::unordered_map<spncci::LGIPair,spncci::ObservableBabySpNCCIHypersectors,boost::hash<LGIPair>>
            BabySpNCCIObservableHypersectorTable;

  void GenerateRecurrenceUnitTensors(
    int Nmax, int N1v,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    std::map<spncci::NnPair,std::set<int>>& operator_subsets_NnpNn
  );
  // Identify unit tensor subspaces for recurrence from list of unit tensor labels 
  // with non-zero rmes between a given lgi pair
  //
  // Unit tensor subspaces whose hyperblocks will be computed by recurrence are
  // stored in operator_subsets_NnpNn by (Nnp,Nn) for Nnp>=Nn

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
  // Transfers seed rmes from unit_tensor_seed_blocks to unit_tensor_hyperblocks
  // and adds conjugate rmes to unit_tensor_hyperblocks_Nn0 for use in recurrence

  void 
  AddNn0BlocksToHyperblocks(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors_Nn0,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_Nn0,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  );
  // Conjugates (0,Nn) hyperblocks computed recursively to 
  // unit_tensor_hyperblocks i.e., add (Nnp,0) to recurrence hyperblocks
  // container

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
  // Recursively computes unit tensor hyperblocks for Nnp>=Nn from 
  // recurrence relation derived in McCoy2018

  bool
  GenerateUnitTensorHyperblocks(
    const spncci::LGIPair& lgi_pair,
    int Nmax, int N1v,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const std::vector<int>& lgi_full_space_index_lookup,
    const spncci::SpNCCISpace& spncci_space,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::KMatrixCache& k_matrix_cache,
    const spncci::KMatrixCache& kinv_matrix_cache,
    u3::UCoefCache& u_coef_cache,
    u3::PhiCoefCache& phi_coef_cache,
    spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
    );

} //namespace 

#endif
