/****************************************************************
recurrence_indexing_spatial.h

    Indexing for SpNCCI recurrence

    SpNCCIIrrep
    ->  sigma
        -> omega (states (n,rho) -> upsilon)

    spncci::spatial::Space() []
    ->spncci::spatial::Sp3RSpace() [sigma]
        ->spncci::spatial::U3Subspace() [omega]
          ->spncci::spatial::U3State() [n,rho] n=Nn(lambda_n,mu_n)/(nx,ny,nz)

    spatial::RecurrenceSpace() []
    ->spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
      ->spatial::RecurrenceNnsumSpace() [Nnsum]
        ->spatial::RecurrenceU3Space() [omega,omega'] (upsilon,upsilon')
          ->spatial::RecurrenceOperatorSubspace() [x0] (rho0)
            ->spatial::RecurrenceOperatorState() [Nbar,Nbar']

    spatial::ContractionSpace() [J0]
    ->spatial::ContractionSp3RSpace() [sigma,sigma',parity_bar]
      ->spatial::ContractionU3Space() [omega,omega'] (upsilon,upsilon')
        ->spatial::ContractionOperatorSubspace() [L0] (kappa0) <-- J0-2<=L0<=J0+2
          ->spatial::ContractionOperatorState() [x0] (rho0)

    spatial::BranchingSpace() [J,J',J0]
    ->spatial::BranchingSp3RSpace() [sigma,sigma',parity_bar]
      ->spatial::BranchingU3Subspace() [omega,omega'] (upsilon,upsilon')
        ->spatial::BranchingState() [L,L'] (kappa,kappa')
    TODO: Find efficient way to actually store in (L,kappa,L',kappa') order

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
#include <memory>      // for std::shared_ptr
#include <unordered_map>
#include <utility>  // for std::forward

#include "am/halfint.h"
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "fmt/format.h"
#include "lgi/lgi.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"
#include "spncci/recurrence_indexing_spin.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"
#include "u3shell/unit_tensor_space_u3s.h"

namespace spncci::spatial
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// spncci::spatial::Space() []
// ->spncci::spatial::Sp3RSpace() [sigma]
//     ->spncci::spatial::U3Subspace() [omega]
//       ->spncci::spatial::U3State() [n,rho] n=Nn(lambda_n,mu_n)/(nx,ny,nz)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

class U3Subspace;
class U3State;
class Sp3RSpace;
class Space;

class U3Subspace
    : public basis::
          BaseDegenerateSubspace<U3Subspace, std::tuple<u3::U3>, U3State, std::tuple<u3::U3>>
{
 public:
  U3Subspace() = default;

  template<typename K1, typename K2>
  U3Subspace(
      const u3::U3& omega,
      const MultiplicityTagged<u3::U3>::vector& nrho_vector,
      K1&& K_matrix__,
      K2&& Kinv_matrix__
    )
      : BaseDegenerateSubspace{omega},
        K_matrix_{std::forward<K1>(K_matrix__)},
        Kinv_matrix_{std::forward<K2>(Kinv_matrix__)}
  {
    for (const auto& [n, rho_max] : nrho_vector) PushStateLabels(n, rho_max);
    assert(K_matrix().rows() == Kinv_matrix().cols());
    assert(
        (nonorthogonal_basis_size() == K_matrix().cols())
        && (nonorthogonal_basis_size() == Kinv_matrix().rows())
      );
  }

  // U3Subspace(const sp3r::U3Subspace& u3subspace);

  u3::U3 omega() const { return std::get<0>(labels()); }
  inline int upsilon_max() const { return dimension(); }

  inline const basis::OperatorBlock<double>& K_matrix() const
  {
    return K_matrix_;
  }
  inline const basis::OperatorBlock<double>& Kinv_matrix() const
  {
    return Kinv_matrix_;
  }

  inline std::size_t dimension() const { return K_matrix().rows(); }
  inline std::size_t nonorthogonal_basis_size() const
  {
    return BaseDegenerateSubspace::dimension();
  }

 private:
  basis::OperatorBlock<double> K_matrix_, Kinv_matrix_;
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

class Sp3RSpace
    : public basis::BaseSpace<Sp3RSpace, U3Subspace, std::tuple<u3::U3>>
{
 public:
  Sp3RSpace() = default;
  Sp3RSpace(const u3::U3& sigma, const int Nn_max);

  u3::U3 sigma() const { return std::get<0>(labels()); }
};

class Space
    : public basis::BaseSpace<Space, Sp3RSpace>
{
 public:
  Space() = default;

  // Construct from list of sigma
  Space(
      const std::vector<u3::U3>& sigma_vector, const HalfInt& Nsigma0, const int Nmax
    );

  // Construct from spin::Space
  template<typename tLGIType>
  Space(const spin::Space<tLGIType>& spin_space, const HalfInt& Nsigma0, const int& Nmax)
  {
    // Nsigma0_=spin_space.Nsigma0();
    for (int i = 0; i < spin_space.size(); ++i)
    {
      const u3::U3& sigma = spin_space.GetSubspace(i).sigma();
      int Nn_max = Nmax - int(sigma.N() - Nsigma0);
      PushSubspace(Sp3RSpace(sigma, Nn_max));
    }
  }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// spatial::RecurrenceSpace() []
// ->spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
//   ->spatial::RecurrenceNnsumSpace() [Nsum]
//     ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
//       ->spatial::RecurrenceOperatorSubspace() [x0] -> rho0
//         ->spatial::RecurrenceOperatorState() [Nbar,Nbar']
/////////////////////////////////////////////////////////////////////////////////////////////////////////
class RecurrenceSpace;
class RecurrenceSp3RSpace;
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
  std::string LabelStr() const { return x0().Str(); }
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
  int Nbarp() const { return std::get<1>(labels()); }
  std::string LabelStr() const { return fmt::format("{} {}", Nbar(), Nbarp()); }
};

// spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
class RecurrenceU3Space
    : public basis::BaseDegenerateSpace<
          RecurrenceU3Space,
          RecurrenceOperatorSubspace,
          std::tuple<u3::U3, u3::U3>
        >
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
  std::string LabelStr() const
  {
    return fmt::format("{} {}", omega_ket().Str(), omega_bra().Str());
  }
};


// spatial::RecurrenceNnsumSpace() [Nsum]
class RecurrenceNnsumSpace
    : public basis::
          BaseDegenerateSpace<RecurrenceNnsumSpace, RecurrenceU3Space, std::tuple<int>>
{
 private:
  using BaseDegenerateSpaceType =
      basis::BaseDegenerateSpace<RecurrenceNnsumSpace, RecurrenceU3Space, std::tuple<int>>;

 public:
  RecurrenceNnsumSpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceNnsumSpace(
      int Nnsum,
      const std::vector<std::tuple<int, int>> u3subspace_index_pairs,
      const Sp3RSpace& sp3r_space_ket,
      const Sp3RSpace& sp3r_space_bra,
      const UnitTensorConstraintParameters& unit_tensor_parameters
    );

  int Nnsum() const { return std::get<0>(labels()); }

  int upsilon_max_ket(const int i) const
  {
    return std::get<0>(upsilon_max_pairs_[i]);
  }
  int upsilon_max_bra(const int i) const
  {
    return std::get<1>(upsilon_max_pairs_[i]);
  }

  std::size_t GetSubspaceOffset(
      std::size_t i, int upsilon_ket = 1, int upsilon_bra = 1
    ) const
  {
    const std::size_t degeneracy_index =
        (upsilon_ket - 1) * upsilon_max_bra(i) + (upsilon_bra - 1);
    return BaseDegenerateSpaceType::GetSubspaceOffset(i, degeneracy_index);
  }

  uint8_t parity_bar() const { return parity_bar_; }

  inline std::string LabelStr() const { return fmt::format("{}", Nnsum()); }

 private:
  uint8_t parity_bar_;
  std::vector<std::tuple<int, int>> upsilon_max_pairs_;
};


// spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
class RecurrenceSp3RSpace
    : public basis::BaseSpace<
          RecurrenceSp3RSpace,
          RecurrenceNnsumSpace,
          std::tuple<u3::U3, u3::U3, uint8_t>
        >
{
 public:
  RecurrenceSp3RSpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceSp3RSpace(
      std::shared_ptr<const Sp3RSpace> sp3r_space_ket_ptr,
      std::shared_ptr<const Sp3RSpace> sp3r_space_bra_ptr,
      const UnitTensorConstraintParameters& unit_tensor_constraints
    );

  u3::U3 sigma_ket() const { return std::get<0>(labels()); }
  u3::U3 sigma_bra() const { return std::get<1>(labels()); }
  uint8_t parity_bar() const { return std::get<2>(labels()); }
  const Sp3RSpace& ket_space() const { return *sp3r_space_ket_ptr_; }
  const Sp3RSpace& bra_space() const { return *sp3r_space_bra_ptr_; }
  const std::shared_ptr<const Sp3RSpace> ket_space_ptr() const
  {
    return sp3r_space_ket_ptr_;
  }
  const std::shared_ptr<const Sp3RSpace> bra_space_ptr() const
  {
    return sp3r_space_bra_ptr_;
  }

  inline std::string LabelStr() const
  {
    return fmt::format(
        "{} {}  {}", sigma_ket().Str(), sigma_bra().Str(), parity_bar()
      );
  }

 private:
  std::shared_ptr<const Sp3RSpace> sp3r_space_ket_ptr_, sp3r_space_bra_ptr_;
};

// spatial::RecurrenceSpace() []
class RecurrenceSpace
    : public basis::BaseSpace<RecurrenceSpace, RecurrenceSp3RSpace>
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


////////////////////////////////////////////////////////////////
// recurrence sectors
////////////////////////////////////////////////////////////////

class RecurrenceU3Sector
    : public basis::BaseSector<RecurrenceU3Space>
{
 public:
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  using BaseSector::BaseSector;

  std::size_t source_subspace_index() const
  {
    return BaseSector::ket_subspace_index();
  }
  std::size_t target_subspace_index() const
  {
    return BaseSector::bra_subspace_index();
  }
  const RecurrenceU3Space& source_subspace() const
  {
    return BaseSector::ket_subspace();
  }
  const RecurrenceU3Space& target_subspace() const
  {
    return BaseSector::bra_subspace();
  }
};

class RecurrenceU3Sectors
    : public basis::BaseSectors<RecurrenceNnsumSpace, RecurrenceNnsumSpace, RecurrenceU3Sector>
{
 public:
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  RecurrenceU3Sectors() = default;

  RecurrenceU3Sectors(
      const RecurrenceSp3RSpace& sp3r_space, int target_Nnsum, int source_Nnsum
    );
};

}  // namespace spncci::spatial

#endif  // RECURRENCE_INDEXING_H_
