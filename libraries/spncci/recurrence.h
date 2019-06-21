/****************************************************************
  recurrence.h

  Unit tensor recursive evaluation
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/16/18 (aem): Created based on unit_tensor.h
  1/30/18 (aem): Extracted all seed and recurrence functions from
    computational control.
  6/21/19 (aem): Factored seed set up and recurrence into separate 
    functions.
****************************************************************/

#ifndef SPNCCI_SPNCCI_RECURRENCE_
#define SPNCCI_SPNCCI_RECURRENCE_
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

void GetLGIPairsForRecurrence(
      const lgi::MultiplicityTaggedLGIVector& lgi_families,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      std::vector<spncci::LGIPair>& lgi_pairs
    );
  // Not currently used.
  //
  // Generate vector of lgi pairs with ket<=bra and ordered based on size of the sum of
  // the dimension of the bra and ket Sp3RSpaces.  

void GetLGIPairsForRecurrence(
      const lgi::MultiplicityTaggedLGIVector& lgi_families,
      const std::vector<int>& lgi_full_space_index_lookup,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      int Nmax,
      std::vector<spncci::LGIPair>& lgi_pairs
    );
  // Alternative sorting of LGI pairs with ordering based on number of steps in recurrence (Nn,max)
  // Checks that file exits before adding it to the set of lgi pairs.
  // Only considers pairs with ket<=bra


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
    basis::OperatorBlocks<double>& unit_tensor_seed_blocks,
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

 // bool
 //    GenerateUnitTensorHyperblocks(
 //      const spncci::LGIPair& lgi_pair,
 //      int Nmax, int N1v,
 //      const lgi::MultiplicityTaggedLGIVector& lgi_families,
 //      const std::vector<int>& lgi_full_space_index_lookup,
 //      const spncci::SpNCCISpace& spncci_space,
 //      const spncci::BabySpNCCISpace& baby_spncci_space,
 //      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
 //      const spncci::KMatrixCache& k_matrix_cache,
 //      const spncci::KMatrixCache& kinv_matrix_cache,
 //      spncci::OperatorBlocks& lgi_transformations,
 //      bool transform_lgi_families,
 //      u3::UCoefCache& u_coef_cache,
 //      u3::PhiCoefCache& phi_coef_cache,
 //      spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
 //      basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
 //      );
 //    //DEPRECATED
 //    // Compute baby spncci hyperblocks for unit tensors from seeds
 //    //
 //    // Output:
 //    //  baby_spncci_hypersectors 
 //    //  unit_tensor_hyperblocks

void DoRecurrenceInitialization(
  int Nmax, int N1v,
  const spncci::LGIPair& lgi_pair,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  spncci::OperatorBlocks& lgi_transformations,
  bool transform_lgi_families,
  std::map<spncci::NnPair,std::set<int>>& unit_tensor_subspace_subsets,
  spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds,
  spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds_conj,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds,
  basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds_conj
  );
  // Setting up seeds for use in recurrence.
  //
  // + Reads in seed labels and RMEs from file
  // + Applies lgi_transformations to do change of basis if transform_lgi_families is true
  // + Seed hypersectors and their conjugates stored in baby_spncci_hypersector_seeds and 
  //   baby_spncci_hypersector_seeds_conj, respectively  
  // + RMEs and conjugate of RMEs stored in unit_tensor_hyperblocks_seeds and 
  //   unit_tensor_hyperblocks_seeds_conj, respectively
  // + Based on list of unit tensors with non-zero matrix elements obtained from file,
  //   creates map containing list of unit_tensor_indices corresponding to unit tensors 
  //   which will have non-zero matrix elements in the recurrence for each pair <Nnp,Nn> 


bool GenerateUnitTensorHyperblocks(
      const spncci::LGIPair& lgi_pair,
      int Nmax, int N1v,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const spncci::KMatrixCache& kinv_matrix_cache,
      std::map<spncci::NnPair,std::set<int>>& unit_tensor_subspace_subsets,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersector_seeds_conj,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks_seeds_conj,
      u3::UCoefCache& u_coef_cache,
      u3::PhiCoefCache& phi_coef_cache,
      spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
    );
    // Compute unit tensor RMEs between many-body basis states with Nnp>=Nn in irreps
    // determined by lgi_pair using recurrence relation derived in McCoy2018.  
    //
    // Output:
    //  baby_spncci_hypersectors 
    //  unit_tensor_hyperblocks


} //namespace 

#endif
