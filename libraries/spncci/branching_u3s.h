/****************************************************************
  branching_u3s.h

  U(3)xS layer of SpNCCI basis branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/31/17 (mac): Extracted from sp_basis (as spncci_branching_u3s).
  2/17/17 (mac): Extract BabySpNCCISubspace to spncci_basis.
  2/19/17 (mac): Rename to branching_u3s.
  5/25/17 (mac): Overhaul implementation of U3S subspaces and store
    parent irrep info.
  7/1/17 (aem): Fix incorrect spin conjugation phase
  1/16/18 (aem) : Add U3S block initialization and printing. Extracted
    from spncci.cpp
  2/1/18 (aem) : Add ObservableHypersectorsU3S
  2/5/18 (aem): Removed obselete U3SSectors
****************************************************************/

#ifndef SPNCCI_SPNCCI_BRANCHING_U3S_H_
#define SPNCCI_SPNCCI_BRANCHING_U3S_H_

#include <unordered_map>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "spncci/spncci_common.h"
#include "spncci/unit_tensor.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"  
#include "u3shell/upcoupling.h"
#include "lgi/lgi.h"

namespace spncci
{
  typedef std::pair<int,int> LGIPair;

  typedef std::vector<std::vector<std::unordered_map<spncci::LGIPair,basis::OperatorHyperblocks<double>,boost::hash<spncci::LGIPair>>>>
      ObservableHyperblocksByLGIPairTable; 

  typedef std::vector<std::unordered_map<spncci::LGIPair,spncci::ObservableBabySpNCCIHypersectors,boost::hash<spncci::LGIPair>>>
      ObservableHypersectorsByLGIPairTable;
  ////////////////////////////////////////////////////////////////
  // basis indexing in U3S scheme for spncci basis branching
  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (omega,S) = U3S
  //
  // state labels within subspace: (baby_spncci_subspace_index)
  //
  //   baby_spncci_subspace_index (int): index of BabySpNCCI subspace
  //   from which state is drawn
  //
  // ----------------
  //
  // The idea is that states are grouped in the hierarchy
  //
  //   subspace: (omega,S)
  //     state: (sigma,Sp,Sn) => actually baby_spncci_subspace_index stored
  //       substates: (gamma,upsilon)
  //
  // Except this has all been made rather opaque, as the state labels
  // are hidden behind an index into the BabySpNCCI space...  Can't we
  // just label a state by its labels, like shown, and provide an
  // auxiliary mechanism for looking up and/or looking up by the
  // corresponding BabySpNCCI subspace???
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Each "state" actually represents a multiplicity of states
  // governed by (gamma_max,upsilon_max).
  //
  // The lgi family index and dimension (upsilon_max*gamma_max) as
  // well as gamma_max and upsilon_max can be extracted from baby
  // spncci.
  //
  // Associated with each subspace is a look-up table which can look up
  // the starting index of the particular "state" in the U3S sector.  
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (omega,S).
  // 
  // NOT TRUE?
  //
  // Aren't they just ordered by first appearance of (omega,S) in the
  // labels for a BabySpNCCI subspace?
  //
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class StateU3S;  // forward declaration (to permit use as "friend" of SubspaceU3S)

  class SubspaceU3S
    : public basis::BaseSubspace<u3::U3S,std::tuple<int>>
    // SubspaceLabelsType (u3shell::U3S): (omega,S)
    // StateLabelsType (int): index into BabySpNCCI space
    {
      public:

      // constructors

      SubspaceU3S() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      SubspaceU3S(const u3::U3S& omegaS, const BabySpNCCISpace& baby_spncci_space);

      // accessors
      u3::U3S omegaS() const {return labels_;}
      u3::U3 omega() const {return omegaS().U3();}
      u3::SU3 x() const {return omegaS().SU3();}
      HalfInt N() const {return omegaS().U3().N();}
      HalfInt S() const {return omegaS().S();}

      int full_dimension() const {return full_dimension_;}
      // int sector_dim() const {return full_dimension();} // DEPRECATED in favor of full_dimension

      // diagnostic output
      std::string Str() const;

      // state indexing lookup -- DEPRECATED
      int sector_index(int state_index) const
      // DEPRECATED -- in favor of StateU3S::starting_index()
      {
        return state_substate_offset_.at(state_index);
      }

      private:

      // state auxiliary data
      //
      // Note: These parallel arrays are getting cumbersome.  Perhaps
      // this information should be bundled into a struct.

      friend class StateU3S;
      std::vector<int> state_substate_offset_;  // starting index, counting (gamma,upsilon) multiplicity
      std::vector<int> state_dimension_;  // number of substates, counting (gamma,upsilon) multiplicity
      std::vector<int> state_gamma_max_;  // gamma_max
      std::vector<u3shell::U3SPN> state_sigmaSPN_;  // Sp irrep symmetry labels

      // dimension, counting (gamma,upsilon) multiplicity of subspace
      int full_dimension_;
      
    };

  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateU3S
    : public basis::BaseState<SubspaceU3S>
  {
    
    public:

    // pass-through constructors

    StateU3S(const SubspaceType& subspace, int& index)
      // Construct state by index.
      : basis::BaseState<SubspaceU3S>(subspace, index) {}

    StateU3S(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
      // Construct state by reverse lookup on labels.
      : basis::BaseState<SubspaceU3S> (subspace, state_labels) 
      {}

    // pass-through accessors
    u3::U3S omegaS() const {return subspace().omegaS();}
    u3::U3 omega() const {return subspace().omega();}
    HalfInt S() const {return subspace().S();}
    HalfInt N() const {return subspace().N();}

    // supplemental data accessors
    int substate_offset() const
    // Provide offset of first substate into fully expanded listing of
    // substates in subspace.
    {
      return subspace().state_substate_offset_[index()];
    }
    int dimension() const
    // Provide number of substates of this composite state.
    {
      return subspace().state_dimension_[index()];
    }
    u3shell::U3SPN sigmaSPN() const
      // Provide full symmetry labels (sigma,Sp,Sn,S) of Sp irrep.
      {
        return subspace().state_sigmaSPN_[index()];
      }
    int gamma_max() const
    // Provide gamma multiplicity of Sp irrep.
    {
      return subspace().state_gamma_max_[index()];
    }

    private:
 
  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class SpaceU3S
    : public basis::BaseSpace<SubspaceU3S>
  {
    
    public:

    // constructor
    SpaceU3S() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceU3S(const BabySpNCCISpace& baby_spncci_space);

    // diagnostic output
    std::string Str() const;

  };

  ////////////////////////////////////////////////////////////////
  // Hypersector
  // omegaS hypersectors
  ////////////////////////////////////////////////////////////////

  class ObservableHypersectorsU3S
    : public basis::BaseHypersectors<spncci::SpaceU3S,u3shell::ObservableSpaceU3S>
  {

  public:
      // constructor

    ObservableHypersectorsU3S() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes

    ObservableHypersectorsU3S(
        const spncci::SpaceU3S& space,
        const u3shell::ObservableSpaceU3S& operator_space
      );
      // Enumerate sector pairs connected by U3S subspaces of 
      // observables
  };


 void InitializeU3SObservableBlocks(
      const spncci::SpaceU3S& space_u3s,
      int num_observables, int num_hw_values,
      const std::vector<spncci::ObservableHypersectorsU3S>& observable_hypersectors_by_observable,
      std::vector<std::vector<spncci::OperatorBlocks>>& observables_blocks_array
    );


  void 
  ContractBabySpNCCIU3S(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      basis::OperatorHyperblocks<double>& observable_hyperblocks
    );

void 
  ContractAndRegroupU3Sv2(
      int num_lgi,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::SpaceU3S& u3s_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableHypersectorsU3S& observable_sectors,
      spncci::OperatorBlocks& observable_blocks
    );

  void 
  ContractAndRegroupU3S(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::SpaceU3S& target_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableHypersectorsU3S& target_sectors_u3s,
      spncci::OperatorBlocks& observables_blocks
    );

  void
  ContractBabySpNCCIU3S(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      basis::OperatorHyperblocks<double>& observable_hyperblocks
    );

  void 
  RegroupU3S(
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::SpaceU3S& u3s_space,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      const basis::OperatorHyperblocks<double>& observable_hyperblocks,
      const spncci::ObservableHypersectorsU3S& observable_sectors,
      spncci::OperatorBlocks& observable_blocks
    );

  void 
  ContractAndRegroupU3S(
      const std::pair<int,int>& lgi_pair,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::SpaceU3S& u3s_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableHypersectorsU3S& observable_sectors,
      spncci::OperatorBlocks& observable_blocks
    );

void 
  ContractBabySpNCCIU3S2(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      basis::OperatorHyperblocks<double>& observable_hyperblocks
    );
}  // namespace

#endif
