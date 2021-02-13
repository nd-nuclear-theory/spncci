
/****************************************************************
  hyperblocks_u3s.h

  U(3)xS layer of SpNCCI basis branching.
                                  
  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  1/2/19 (aem): Created based on branching_u3s
  6/21/19 (aem): Updated ComputeManyBodyRMEs to use new seed set
    up and recurrence functions. 
****************************************************************/

#ifndef SPNCCI_SPNCCI_HYPERBLOCKS_U3S_H_
#define SPNCCI_SPNCCI_HYPERBLOCKS_U3S_H_

#include "spncci/spncci_basis.h"
#include "spncci/recurrence.h"
#include "spncci/parameters.h"
namespace spncci
{

  typedef std::vector<spncci::ObservableBabySpNCCIHypersectors> ObservableHypersectorsTable;
  typedef std::tuple<int,int,int,int> ObservableHypersectorLabels;
  typedef std::vector< std::vector<basis::OperatorHyperblocks<double>>> ObservableHyperblocksTable;


	void ContractBabySpNCCISymmetricHypersectors(
	  const spncci::LGIPair& lgi_pair,
	  int num_observables, int num_hw_values,
	  const spncci::BabySpNCCISpace& baby_spncci_space,
	  const std::vector<u3shell::ObservableSpaceU3S>& observable_spaces,
	  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
	  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors1,
	  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors2,
	  const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks1,
	  const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks2,
	  const std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>>& observables_relative_rmes,
	  spncci::ObservableHypersectorsTable& observable_hypersectors_table,
	  spncci::ObservableHyperblocksTable& observable_hyperblocks_table
  	);

void GetBabySpNCCIHyperBlocks(
  const int observable_index,
  const int hw_index,
  const spncci::LGIPair& lgi_pair,
  std::vector<spncci::ObservableHypersectorLabels>& list_baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks
  );

void GetOperatorTile(
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::ObservableSpaceU3S& observable_space,
  const spncci::SubspaceSpBasis& spbasis_subspace_bra,
  const spncci::SubspaceSpBasis& spbasis_subspace_ket,
  const std::vector<int>& offsets_bra_subspace,
  const std::vector<int>& offsets_ket_subspace,
  const HalfInt& J0, const HalfInt& Jp, const HalfInt& J,
  const int hw_index,
  const int observable_index,
  const spncci::LGIPair& lgi_pair,
  u3::WCoefCache& w_cache,
  const std::vector<spncci::ObservableHypersectorLabels>& list_baby_spncci_hypersectors,
  basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks,
  spncci::OperatorBlock& tile
  );

  void ConstructSymmetricOperatorMatrix(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    const HalfInt& J0,
    const spncci::SpaceSpBasis& spbasis_bra, //For a given J
    const spncci::SpaceSpBasis& spbasis_ket, //For a given J
    const std::vector<spncci::LGIPair>& lgi_pairs,
    int observable_index, int hw_index,
    spncci::OperatorBlock& operator_matrix
  );

void ComputeManyBodyRMEs(
  const spncci::RunParameters& run_parameters,
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::SpNCCISpace& spncci_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const std::vector<u3shell::ObservableSpaceU3S>& observable_spaces,
  const std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>>& observables_relative_rmes,
  const spncci::KMatrixCache& k_matrix_cache,
  const spncci::KMatrixCache& kinv_matrix_cache,
  spncci::OperatorBlocks& lgi_transformations,
  u3::UCoefCache& u_coef_cache,
  u3::PhiCoefCache& phi_coef_cache,
  const spncci::LGIPair& lgi_pair
  );
  //Calculates many-body matrix elements of operators defined by observables_relative_rmes.  
  //
  // Calculates unit tensor hyperblocks for given lgi pair and it's conjugate
  // Contracts unit tensor hyperblocks with relative rmes of operators to get operator rme hyperblocks
  // Writes hyperblocks to file

}//namespace
#endif