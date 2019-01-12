
/****************************************************************
  hyperblocks_u3s.h

  U(3)xS layer of SpNCCI basis branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/2/19 (aem): Created based on branching_u3s 
****************************************************************/

#ifndef SPNCCI_SPNCCI_HYPERBLOCKS_U3S_H_
#define SPNCCI_SPNCCI_HYPERBLOCKS_U3S_H_

#include "spncci/spncci_basis.h"
#include "spncci/unit_tensor.h"

#include "spncci/branching2.h" //SpaceSpBasis
#include "spncci/branching_u3s.h" //ObservableHypersectorsTable
namespace spncci
{

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

}//namespace
#endif