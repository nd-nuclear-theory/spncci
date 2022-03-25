/****************************************************************
  sp3r.h

  Sp(3,R) labeling and branching.

  Anna E. McCoy [1,2,3] and Mark A. Caprio[1]
  [1] University of Notre Dame
  [2] TRIUMF
  [3] Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/9/16 (aem,mac): Created based on prototype spstates.py, sp3r.py,
    and coefficients.py.
  2/1/17 (mac): Rename DebugString to DebugStr.
  7/1/17 (aem): Add modified branching rule for sp3r->u3 for A<6
  9/27/17 (aem):
    + Updated modified branching rule
    + Added upsilon_max accessor for U3Subspace giving upsilon max
      U3Subspace size still gives number of (n,rho) pairs in
      corresponding U3boson space
    + Broke off sp3r coefficients into separate module
  11/4/21 (aem): Add new functions checking if Sp(3,R)->U(3) branching
      must be restricted.
  3/11/22 (aem): 
    + Fixed modified branching rule for Sp3R->U3 for A<6
    + Switch U3Subspace to be a BaseDegenerateSubspace with degeneracy rho_max
    + Add option for labels only construction of Sp3RSpace
    + Store K matrices with U3Subspace when constructing full space
  3/24/22 (aem):
    + Changes state U3Subspace to SO3States.
    + Raising polynomial information now stored by shared ptr to
      u3boson::U3Subspace which is accessed by nonorthogonal_basis()
****************************************************************/

#ifndef SP3R_H_
#define SP3R_H_

#include <map>
#include <string>

#include "basis/basis.h"
#include "basis/operator.h"
#include "basis/degenerate.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"

namespace sp3r
{

  // TODO: Move into vcs?  Or separate sp3r_utils file?
  // Used in vcs.cpp but, vcs K matrix functions used in sp3r.cpp
  bool IsUnitary(const u3::U3& sigma);
  // Check if sigma is label of unitary Sp(3,R) irrep
  // based on the criteria given in jpa-18-1985-939-Rowe.

  bool ModifySp3RBranching(const u3::U3& sigma);
  // Returns true if Sp(3,R)->U(3) branching obtained by coupling
  // Sp(3,R) raising polynomials onto sigma must be modified.


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Class for carrier space of Sp(3,R)>U(3)>SO(3) irrep
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // sp3r::Sp3RSpace() 
  //  ->sp3r::U3Subspace()
  //    ->sp3r::SO3State()
  //
  //  sp3r::U3Subspace also contains:
  //    + K and Kinv matrices for basis orthogonalization.  Only stored if 
  //      constructor flag subspace_labels_only = false. 
  //      Elements of K_matrix are (sigma upsilon omega||K||sigma n rho omega) and
  //      elements of Kinv_matrix are (sigma n rho omega||Kinv||sigma upsilon omega).
  //
  //    + shared pointer to u3boson::U3Subspace containing u3boson subspaces
  //      which corresponds to the non-orthogonal Sp(3,R) basis
  //      Within u3boson subspace, "states" are raising polynomial labels n
  //      with degeneracy rho_max.
  //
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // space labels: sigma (u3::U3)
  // space truncation: Nn,max (integer)
  //
  // subspace labels: omega (u3::U3)
  // subspace degeneracy: upsilon_max
  // 
  // state labels: L (unsigned int)
  // state degeneracy: kappa_max
  //
  // state labels within subspace: L (unsigned int) with degeneracy kappa_max
  //    state labels and multiplicities only stored if space constructor flag
  //    subspace_labels_only = false.
  //
  // Within a space, the subspaces are ordered by:
  //   -- "canonically" increasing omega
  //      which is defined for us as lexicographical by N(lambda,mu)
  //
  // Within a subspace, the states are ordered by:
  //   -- "canonically" increasing L

  class U3Subspace;
  class Sp3RSpace;
  class SO3State;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class U3Subspace
    : public
      basis::BaseDegenerateSubspace<
        U3Subspace,
        std::tuple<u3::U3>,
        SO3State,
        std::tuple<unsigned int>
      >
  {
  public:

    U3Subspace() = default;
    // constructor

    inline U3Subspace(const u3::U3& omega, unsigned int upsilon_max)
      : BaseDegenerateSubspace{omega}, upsilon_max_{upsilon_max}
     {}
    // This is a lightweight constructor which only stores the subspace
    // labels without populating the subspace with states and storing the
    // Kmatrices and raising polynomial quantum numbers.

    template<typename K1, typename K2>
    U3Subspace(
      const u3::U3& omega,
      int upsilon_max,
      std::shared_ptr<const u3boson::U3Subspace> u3boson_ptr,
      K1&& K_matrix__,
      K2&& Kinv_matrix__
    );
    // Full subspace constructor

    ////////////////////////////////////////////////////////////////////////
    // accessors which can always be used
    const u3::U3 omega() const {return std::get<0>(labels());}
    const u3::U3 U3() const {return omega(); } //Deprecate?
    inline unsigned int upsilon_max() const {return upsilon_max_;}
    inline std::size_t dimension() const {return upsilon_max();}
    std::string LabelStr() const {return omega().Str();}

    // Accessors which should only be used when full subspace constructed.
    inline const basis::OperatorBlock<double>& K_matrix() const {return K_matrix_;}
    inline const basis::OperatorBlock<double>& Kinv_matrix() const {return Kinv_matrix_;}
    const u3boson::U3Subspace& nonorthogonal_basis() const {return *u3boson_ptr_;}
    inline std::size_t nonorthogonal_basis_dimension() const {return K_matrix().cols();}
    std::string DebugStr() const;

  private:
    unsigned int upsilon_max_;
    basis::OperatorBlock<double> K_matrix_, Kinv_matrix_;
    std::shared_ptr<const u3boson::U3Subspace> u3boson_ptr_;
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class SO3State
      : public basis::BaseDegenerateState<U3Subspace>
    {
     public:

      SO3State() = default;
      // pass-through constructors

      SO3State(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
          : basis::BaseDegenerateState<U3Subspace>(subspace, index)
      {}

      SO3State(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
      // Construct state by reverse lookup on labels.
          : basis::BaseDegenerateState<U3Subspace>(subspace, state_labels)
      {}

      // pass-through accessors for subspace labels
      unsigned int L() const { return std::get<0>(labels()); }
      unsigned int kappa_max() const { return subspace().GetStateDegeneracy(index());}
      // private:
    };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class Sp3RSpace
      : public basis::BaseSpace<Sp3RSpace, U3Subspace>
  // : public basis::BaseSpace<Sp3RSpace, U3Subspace, u3::U3>
  {
   public:
    Sp3RSpace() = default;
    Sp3RSpace(const u3::U3& sigma, int Nn_max, const bool subspace_labels_only = false);

    // accessors
    u3::U3 sigma() const {return sigma_;}

    // u3::U3 sigma() const {return std::get<0>(labels());}
    unsigned int Nn_max() const {return Nn_max_;}

    // diagnostic output
    std::string DebugStr() const;

  private:
    unsigned int Nn_max_;
    u3::U3 sigma_;
  };


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Sectors: U3Subspaces connected by operator w0
  class Sp3RSectors
    : public basis::BaseSectors<Sp3RSpace>
  {
    public:

    // Default constructor
    Sp3RSectors() = default;

    // Constructor
    Sp3RSectors(
        const Sp3RSpace& space,
        const u3::U3& omega0,
        const bool& su3_generator=false
      );

  private:
    u3::U3 omega0_;
  };

}  // namespace

#endif
