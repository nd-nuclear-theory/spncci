/****************************************************************
recurrence_indexing_spatial.h

    Indexing for SpNCCI recurrence

    SpNCCIIrrep
    ->  sigma
        -> omega (states (n,rho) -> upsilon)

    spncci::spatial::Space() []
    ->spncci::spatial::LGISpace() [sigma]
        ->spncci::spatial::U3Subspace() [omega]
          ->spncci::spatial::U3State() [n,rho] n=Nn(lambda_n,mu_n)/(nx,ny,nz)

    spatial::RecurrenceSpace() []
    ->spatial::RecurrenceLGISpace() [sigma,sigma',parity_bar]
      ->spatial::RecurrenceNnsumSpace() [Nsum]
        ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
          ->spatial::RecurrenceOperatorSubspace() [x0] ->rho0_max
            ->spatial::RecurrenceOperatorState() [Nbar,Nbar']

  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT

  + 06/24/21 (aem): Created.
  + 08/05/21 (pjf): Split into separate headers for spin and spatial.
****************************************************************/

#ifndef RECURRENCE_INDEXING_SPATIAL_H_
#define RECURRENCE_INDEXING_SPATIAL_H_

#include <array>
#include <functional>  // for std::hash
#include <unordered_map>


// #include "basis/hypersector.h"
#include "am/halfint.h"
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "lgi/lgi.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"
#include "spncci/recurrence_indexing_spin.h"
// #include "spncci/spncci_common.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"
#include "u3shell/unit_tensor_space_u3s.h"
// #include "u3shell/upcoupling.h"

namespace spncci::spatial
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// spncci::spatial::Space() []
// ->spncci::spatial::LGISpace() [sigma]
//     ->spncci::spatial::U3Subspace() [omega]
//       ->spncci::spatial::U3State() [n,rho] n=Nn(lambda_n,mu_n)/(nx,ny,nz)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class U3Subspace;
class U3State;
class LGISpace;
class Space;

class U3Subspace
    : public basis::BaseDegenerateSubspace<U3Subspace, std::tuple<u3::U3>, U3State, std::tuple<u3::U3>>
{
 public:
  U3Subspace() = default;

  U3Subspace(const u3::U3& omega, const MultiplicityTagged<u3::U3>::vector& nrho_vector);

  U3Subspace(const sp3r::U3Subspace& u3subspace);

  u3::U3 omega() const { return std::get<0>(labels()); }
  inline int upsilon_max() const { return dimension(); }
};

class U3State
    : public basis::BaseDegenerateState<U3Subspace>
{
 public:
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

  // private:
};

class LGISpace
    : public basis::BaseSpace<LGISpace, U3Subspace, std::tuple<u3::U3>>
{
 public:
  LGISpace() = default;
  LGISpace(const u3::U3& sigma, const int Nn_max);

  u3::U3 sigma() const { return std::get<0>(labels()); }
};

class Space
    : public basis::BaseSpace<Space, LGISpace>
{
 public:
  Space() = default;

  // Construct from list of sigma
  Space(const std::vector<u3::U3>& sigma_vector, const HalfInt& Nsigma0, const int Nmax);

  // Construct from spin::Space
  template<typename tLGIType>
  Space(const spin::Space<tLGIType>& spin_space, const HalfInt& Nsigma0, const int& Nmax)
  {
    // Nsigma0_=spin_space.Nsigma0();
    for (int i = 0; i < spin_space.size(); ++i)
    {
      const u3::U3& sigma = spin_space.GetSubspace(i).sigma();
      int Nn_max = Nmax - int(sigma.N() - Nsigma0);
      PushSubspace(LGISpace(sigma, Nn_max));
    }
  }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// spatial::RecurrenceSpace() []
// ->spatial::RecurrenceLGISpace() [sigma,sigma']
//   ->spatial::RecurrenceNnsumSpace() [Nsum]
//     ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
//       ->spatial::RecurrenceOperatorSubspace() [x0] ->rho0_max
//         ->spatial::RecurrenceOperatorState() [Nbar,Nbar']
/////////////////////////////////////////////////////////////////////////////////////////////////////////
class RecurrenceSpace;
class RecurrenceLGISpace;
class RecurrenceNnsumSpace;
class RecurrenceU3Space;
class RecurrenceOperatorSubspace;
class RecurrenceOperatorState;

struct UnitTensorConstraintParameters
{
  UnitTensorConstraintParameters(int N1v_, HalfInt Nsigma0_, uint8_t parity_bar_)
      : N1v{N1v_}, Nsigma0{Nsigma0_}, parity_bar{parity_bar_}
  {}

  int N1v;
  HalfInt Nsigma0;
  uint8_t parity_bar;
};


class RecurrenceOperatorSubspace
    : public basis::BaseSubspace<
          RecurrenceOperatorSubspace,
          std::tuple<u3::SU3>,
          RecurrenceOperatorState,
          std::tuple<int, int>
        >
{
 public:
  RecurrenceOperatorSubspace() = default;

  RecurrenceOperatorSubspace(
      const u3::SU3& x0, const std::vector<std::tuple<int, int>>& Nbar_pairs
    );

  u3::SU3 x0() const { return std::get<0>(labels()); }
};

class RecurrenceOperatorState
    : public basis::BaseState<RecurrenceOperatorSubspace>
{
 public:
  RecurrenceOperatorState(const SubspaceType& subspace, std::size_t index)
  // Construct state by index.
      : BaseState{subspace, index}
  {}

  RecurrenceOperatorState(
      const SubspaceType& subspace,
      const typename SubspaceType::StateLabelsType& state_labels
    )
  // Construct state by reverse lookup on labels.
      : BaseState{subspace, state_labels}
  {}

  // // pass-through accessors for subspace labels
  int Nbar() const { return std::get<0>(labels()); }
  int Nbar_p() const { return std::get<1>(labels()); }
};

// spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
class RecurrenceU3Space
    : public basis::BaseDegenerateSpace<RecurrenceU3Space, RecurrenceOperatorSubspace, std::tuple<u3::U3, u3::U3>>
{
 public:
  RecurrenceU3Space() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceU3Space(
      const std::tuple<u3::U3, u3::U3>& omega_pair,
      const UnitTensorConstraintParameters& unit_tensor_parameters
    );

  u3::U3 omega_ket() const { return std::get<0>(labels()); }
  u3::U3 omega_bra() const { return std::get<1>(labels()); }
};


// spatial::RecurrenceNnsumSpace() [Nsum]
class RecurrenceNnsumSpace
    : public basis::BaseDegenerateSpace<RecurrenceNnsumSpace, RecurrenceU3Space, std::tuple<int>>
{
 public:
  RecurrenceNnsumSpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceNnsumSpace(
      int Nnsum,
      const std::vector<std::tuple<int, int>> u3subspace_index_pairs,
      const LGISpace& lgi_space_ket,
      const LGISpace& lgi_space_bra,
      const UnitTensorConstraintParameters& unit_tensor_parameters
    );

  int Nnsum() const { return std::get<0>(labels()); }

  int upsilon_max_ket(const int i) const
  {
    return std::get<0>(upsilon_pairs_[i]);
  }
  int upsilon_max_bra(const int i) const
  {
    return std::get<1>(upsilon_pairs_[i]);
  }
  uint8_t parity_bar() const { return parity_bar_; }

 private:
  uint8_t parity_bar_;
  std::vector<std::tuple<int, int>> upsilon_pairs_;
};


// spatial::RecurrenceLGISpace() [sigma,sigma',parity_bar]
class RecurrenceLGISpace
    : public basis::BaseSpace<RecurrenceLGISpace, RecurrenceNnsumSpace, std::tuple<u3::U3, u3::U3, uint8_t>>
{
 public:
  RecurrenceLGISpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceLGISpace(
      const LGISpace& lgi_space_ket,
      const LGISpace& lgi_space_bra,
      const UnitTensorConstraintParameters& unit_tensor_constraints
    );

  u3::U3 sigma_ket() const { return std::get<0>(labels()); }
  u3::U3 sigma_bra() const { return std::get<1>(labels()); }
  uint8_t parity_bar() const { return std::get<2>(labels()); }
};

// spatial::RecurrenceSpace() []
class RecurrenceSpace
    : public basis::BaseSpace<RecurrenceSpace, RecurrenceLGISpace>
{
 public:
  RecurrenceSpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceSpace(
      const spatial::Space& space_ket,
      const spatial::Space& space_bra,
      const int& N1v,
      const HalfInt& Nsigma0
    );

  // private:
};

}  // namespace spncci::spatial

#endif  // RECURRENCE_INDEXING_H_