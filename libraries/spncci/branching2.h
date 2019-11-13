/****************************************************************
  branching2.h

  Basis definitions for U(3)xS, LxS, and J branchings of SpNCCI basis.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/18/18 (aem): Created based on branching.h
****************************************************************/

#ifndef SPNCCI_BRANCHING2_H_
#define SPNCCI_BRANCHING2_H_

#include "basis/degenerate.h"
#include "sp3rlib/u3.h"
#include "spncci/spncci_basis.h"
#include "u3shell/u3spn_scheme.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // convenience typedef for (L,S)
  ////////////////////////////////////////////////////////////////

  typedef std::pair<u3::U3,int> omegaLLabels;


  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to J level
  ////////////////////////////////////////////////////////////////
  //
  //   subspace: (sigma,Sp,Sn,S) irrep_family_index
  //     state: (omega,L)
  //       substates: (kappa,gamma,upsilon)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (sigma,Sp,Sn,S) => u3shell::U3SPN
  //
  // state labels within subspace: (omega,L) => omegaLLabels
  //
  // substate labels (implied): (kappa,gamma,upsilon)
  //
  //   (See BabySpNCCI docstring in spncci_basis for definitions of
  //   these basis labels.)
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, states are ordered by first appearance of
  // omega in irrep_family (sigma,Sp,Sn,S) then by L in order of
  // appearence in u3::BranchingSO3(omega).
  //
  // Degeneracy is kappa_max*gamma_max*upsilon_max.
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by appearance in
  // babyspncci
  //
  // ALTERNATE CODE AVAILABLE BUT NOT SELECTED:
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (omega,S).

  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class SubspaceSpBasis
    : public basis::BaseDegenerateSubspace<int,omegaLLabels>
    // SubspaceLabelsType (int) : irrep_family_index
    //     Formerly   SubspaceLabelsType (u3::U3SPN): (sigma,Sp,Sn,S)
    // StateLabelsType (omegaLLabels): (omega,L)
    {
      public:

      // constructors

      SubspaceSpBasis() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      SubspaceSpBasis(
        const HalfInt& J,
        const u3shell::U3SPN& sigmaSPN,
        int irrep_family_index,
        const BabySpNCCISpace& baby_spncci_space
      );

      // subspace label accessors
      // u3shell::U3SPN sigmaSPN() const {return labels_;}
      u3shell::U3SPN sigmaSPN() const {return sigmaSPN_;}

      HalfInt N() const {return sigmaSPN().U3().N();}
      u3::U3 sigma() const {return sigmaSPN().U3();}
      HalfInt S() const {return sigmaSPN().S();}

      // int irrep_family_index() const {return irrep_family_index_;}
      int irrep_family_index() const {return labels_;}
      int gamma_max() const {return gamma_max_;}
      int dimension() const {return dimension_;}

      // state auxiliary data accessors
      const std::vector<int>& state_kappa_max() const {return state_kappa_max_;}
      const std::vector<int>& state_upsilon_max() const {return state_upsilon_max_;}
      const std::vector<int>& state_baby_spncci_subspace_index() const {return state_baby_spncci_subspace_index_;}

      // diagnostic output
      std::string LabelStr() const;
      std::string DebugStr() const;

      private:

      int gamma_max_;
      int irrep_family_index_;
      int dimension_;
      u3shell::U3SPN sigmaSPN_;

      // state auxiliary data
      std::vector<int> state_kappa_max_;
      std::vector<int> state_upsilon_max_;
      std::vector<int> state_baby_spncci_subspace_index_;
    };

  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateSpBasis
    : public basis::BaseDegenerateState<SubspaceSpBasis>
  {

    public:

    // pass-through constructors

    StateSpBasis(const SubspaceType& subspace, int& index)
      // Construct state by index.
      : basis::BaseDegenerateState<SubspaceSpBasis>(subspace,index) {}

    StateSpBasis(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
      // Construct state by reverse lookup on labels.
      : basis::BaseDegenerateState<SubspaceSpBasis>(subspace,state_labels)
      {}

    // pass-through accessors for subspace labels
    u3shell::U3SPN sigmaSPN() const {return subspace().sigmaSPN();}
    u3::U3 sigma() const {return subspace().sigma();}
    HalfInt S() const {return subspace().S();}
    HalfInt N() const {return subspace().N();}


    // state label accessors
    omegaLLabels omegaL() const {return labels();}
    u3::U3 omega() const {return labels().first;}
    int L() const {return labels().second;}
    int Nn() const
    {
      return int(omega().N()-N());
    }

    // diagnostic output
    // std::string LabelStr() const;

    // state auxiliary data accessors
    int upsilon_max() const
    {
      return subspace().state_upsilon_max()[index()];
    }
    int kappa_max() const
    {
      return subspace().state_kappa_max()[index()];
    }

    int baby_spncci_subspace_index() const
    {
      return subspace().state_baby_spncci_subspace_index()[index()];
    }

    private:

  };

  // ////////////////////////////////////////////////////////////////
  // // space
  // ////////////////////////////////////////////////////////////////

  class SpaceSpBasis
    : public basis::BaseDegenerateSpace<SubspaceSpBasis>
  {

    public:

    // constructor
    SpaceSpBasis() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J);

    //Alternate constructor that constructs a subspace of the full space base on irrep_family_subset
    SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J, std::set<int>irrep_family_subset);

    HalfInt J() const {return J_;}

    // diagnostic output
    std::string DebugStr(bool show_subspaces=false) const;

    private:
      HalfInt J_;

  };


  ////////////////////////////////////////////////////////////////
  // sectors
  ////////////////////////////////////////////////////////////////

  class SectorsSpBasis
    : public basis::BaseSectors<SpaceSpBasis>
  {

    public:

    // constructor

    SectorsSpBasis() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SectorsSpBasis(
        const SpaceSpBasis& space,
        HalfInt J0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial character.

  };


  void GetSpBasisOffsets(
    const spncci::SpaceSpBasis& spbasis,
    std::vector<std::vector<int>>& offsets
    );

  void ConstructOperatorMatrix(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    const HalfInt& J0,
    // u3::WCoefCache& w_cache,
    const spncci::SpaceSpBasis& spbasis_bra, //For a given J
    const spncci::SpaceSpBasis& spbasis_ket, //For a given J
    std::vector<int>& nums_lgi_pairs,
    int observable_index, int hw_index,
    spncci::OperatorBlock& operator_matrix
  );






}  // namespace

#endif
