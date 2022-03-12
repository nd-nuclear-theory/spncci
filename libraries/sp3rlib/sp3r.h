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

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep (space) -> U(3) irrep (subspace)
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // space labels: sigma (u3::U3)
  // space truncation: Nn,max (integer)
  // subspace labels: omega (u3::U3)
  // state labels within subspace: n (u3::U3) with degeneracy rho_max (integer)
  //
  // Within a space, the subspaces are ordered by:
  //   -- "canonically" increasing omega
  //      which is defined for us as lexicographical by N(lambda,mu)
  //
  // Within a subspace, the states are ordered by:
  //   -- "canonically" increasing n
  //      which is defined for us as lexicographical by N(lambda,mu)
  //   -- numerically increasing rho

  class U3Subspace;
  class Sp3RSpace;
  class U3State;

  ////////////////////////////////////////////////////////////////
  // U(3) subspace
  ////////////////////////////////////////////////////////////////
  
class U3Subspace
    : public basis::
          BaseDegenerateSubspace<U3Subspace, std::tuple<u3::U3>, U3State, std::tuple<u3::U3>>
  {
  public:

    U3Subspace() = default;
    // constructor
    U3Subspace(const u3::U3& omega, int upsilon_max)
      : BaseDegenerateSubspace{omega}, upsilon_max_{upsilon_max}
     {}
    // U3Subspace(const u3::U3& omega, int upsilon_max);
    // Construct U(3) subspace.
    //
    // This is a lightweight constructor which only stores the labels,
    // without populating the subspace with states.
    //
    // Arguments:
    //   omega (u3::U3) : labels for subspace


    // U3Subspace(
    //     const u3::U3& omega,
    //     int upsilon_max,
    //     const MultiplicityTagged<u3::U3>::vector& state_set
    //   )
    //   : U3Subspace(omega, upsilon_max)
    // {Init(state_set);}

    // void Init(const MultiplicityTagged<u3::U3>::vector& state_set);
    // Alternative constructor from list of (n,rho) states

    ////////////////////////////////////////////////////////////////////////
    template<typename K1, typename K2>
    U3Subspace(
      const u3::U3& omega,
      int upsilon_max,
      const u3boson::U3Subspace& u3boson_subpace,
      K1&& K_matrix__,
      K2&& Kinv_matrix__
    )
    : BaseDegenerateSubspace{omega},
      upsilon_max_{upsilon_max},
      K_matrix_{std::forward<K1>(K_matrix__)},
      Kinv_matrix_{std::forward<K2>(Kinv_matrix__)}
    {
      assert(K_matrix_.rows()==Kinv_matrix_.cols());
      assert(K_matrix_.cols()==Kinv_matrix_.rows());
      Init(u3boson_subpace);
    }

    void Init(const u3boson::U3Subspace& u3boson_subpace);

    ////////////////////////////////////////////////////////////////////////

    // accessors
    const u3::U3 U3() const {return std::get<0>(labels());}
    const u3::U3 omega() const { return U3(); }
    std::string LabelStr() const {return omega().Str();}
    int upsilon_max() const {return upsilon_max_;}
    std::string DebugStr() const;

    inline std::size_t nonorthogonal_basis_size() const
    {
      return BaseDegenerateSubspace::dimension();
    }

    // Currently K_matrix_ and Kinv_matrix_ only stored when
    // basis constructed using U3BosonSpace. 
    inline const basis::OperatorBlock<double>& K_matrix() const
      {
        return K_matrix_;
      }
    inline const basis::OperatorBlock<double>& Kinv_matrix() const
      {
        return Kinv_matrix_;
      }

  private:
    int upsilon_max_;
    basis::OperatorBlock<double> K_matrix_, Kinv_matrix_;

  };

  ////////////////////////////////////////////////////////////////
  // U3State [n] -> rho_max
  ////////////////////////////////////////////////////////////////

  class U3State
      : public basis::BaseDegenerateState<U3Subspace>
    {
     public:

      U3State() = default;
      // pass-through constructors

      U3State(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
          : basis::BaseDegenerateState<U3Subspace>(subspace, index)
      {}

      U3State(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
      // Construct state by reverse lookup on labels.
          : basis::BaseDegenerateState<U3Subspace>(subspace, state_labels)
      {}

      // pass-through accessors for subspace labels
      u3::U3 n() const { return std::get<0>(labels()); }
      int rho_max() const { return subspace().GetStateDegeneracy(index()); }
      MultiplicityTagged<u3::U3> n_multiplicity_tagged() const { return {n(),rho_max()}; }

    };

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) space
  ////////////////////////////////////////////////////////////////

class Sp3RSpace
    : public basis::BaseSpace<Sp3RSpace, U3Subspace>
// : public basis::BaseSpace<Sp3RSpace, U3Subspace, u3::U3>
  {
   public:
    Sp3RSpace() = default;
    Sp3RSpace(const u3::U3& sigma, int Nn_max, const bool subspace_labels_only = false);

    // diagnostic output
    std::string DebugStr() const;

    // accessors
    u3::U3 sigma() const {return sigma_;}
    // u3::U3 sigma() const {return std::get<0>(labels());}
    int Nn_max() const {return Nn_max_;}

  private:
    int Nn_max_;
    u3::U3 sigma_;
  };



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
