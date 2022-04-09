/****************************************************************
recurrence_indexing_spin.h

    Indexing for SpNCCI recurrence

    SpNCCIIrrep
    ->  sigma
        -> omega (states (n,rho) -> upsilon)

    spin::Space() []
    ->spin::LGISpace() [sigma]
      ->spin::SpinSubspace() [S]
        ->spin::[Iso]SpinState() [SpSn]/[T] -> gamma_max

    spin::RecurrenceSpace() []
      spin::RecurrenceLGISpace() [sigma,sigma',exchange_symm_bar]
      ->spin::RecurrenceSpinSpace() [S,S']
        -> spin::RecurrenceSpinSubspace() [Sp,Sn,Sp',Sn']/[T,T']->(gamma,gamma')
          ->spin::RecurrenceOperatorState() [operator_index]

    operator_index -> (S0,T0,Sbar,Sbar',Tbar,Tbar') bit representation

  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT

  06/24/21 (aem) : Created.
****************************************************************/

#ifndef RECURRENCE_INDEXING_SPIN_H_
#define RECURRENCE_INDEXING_SPIN_H_

#include <array>
#include <functional>  // for std::hash
#include <map>
#include <tuple>
#include <utility>

#include "am/am.h"
#include "am/halfint.h"
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "fmt/format.h"
// #include "lgi/lgi.h"
#include "sp3rlib/multiplicity_tagged.h"
#include "sp3rlib/u3.h"

namespace spncci
{
namespace spin
{
///////////////////////////////////////////////////////////////////////////////////////////
// spin::Space() []
// ->spin::LGISpace() [sigma]
//   ->spin::SpinSubspace() [S]
//     ->spin::[Iso]SpinState() [SpSn]/[T] -> gamma_max
///////////////////////////////////////////////////////////////////////////////////////////

// declarations
template<typename tLGIType> class Space;
template<typename tLGIType> class LGISpace;
template<typename tLGIType> class SpinSubspace;
template<typename tLGIType> class SpinState;

template<typename tLGIType>
class SpinSubspace
    : public basis::BaseDegenerateSubspace<
          SpinSubspace<tLGIType>,               /* derived type */
          std::tuple<HalfInt>,                  /* subspace labels */
          SpinState<tLGIType>,                  /* state type */
          typename tLGIType::UpstreamLabelsType /* state labels */
        >
{
 private:
  using BaseDegenerateSubspaceType = basis::BaseDegenerateSubspace<
      SpinSubspace<tLGIType>,
      std::tuple<HalfInt>,
      SpinState<tLGIType>,
      typename tLGIType::UpstreamLabelsType
    >;
  using UpstreamLabelsType = typename tLGIType::UpstreamLabelsType;

 public:
  using LGIType = tLGIType;
  using SubspaceLabelsType =
      typename BaseDegenerateSubspaceType::SubspaceLabelsType;
  using StateLabelsType = typename BaseDegenerateSubspaceType::StateLabelsType;

  SpinSubspace() = default;

  SpinSubspace(
      const SubspaceLabelsType& S,
      const typename MultiplicityTagged<StateLabelsType>::vector& upstreams_vector
    )
      : BaseDegenerateSubspaceType{S}
  {
    for (const auto& [upstreams, gamma_max] : upstreams_vector)
      BaseDegenerateSubspaceType::PushStateLabels(upstreams, gamma_max);
  }

  HalfInt S() const{return std::get<0>(BaseDegenerateSubspaceType::labels());}

  std::string LabelStr() const {return S().Str();}

 private:
};

template<typename tLGIType>
class SpinState
    : public basis::BaseDegenerateState<SpinSubspace<tLGIType>>
{
 private:
  using BaseDegenerateStateType =
      basis::BaseDegenerateState<SpinSubspace<tLGIType>>;

 public:
  using LGIType = tLGIType;
  using SubspaceType = typename BaseDegenerateStateType::SubspaceType;
  using StateLabelsType = typename BaseDegenerateStateType::StateLabelsType;

  // pass-through constructors

  SpinState(const SubspaceType& subspace, std::size_t index)
  // Construct state by index.
      : BaseDegenerateStateType{subspace, index}
  {}

  SpinState(const SubspaceType& subspace, const StateLabelsType& state_labels)
  // Construct state by reverse lookup on labels.
      : BaseDegenerateStateType{subspace, state_labels}
  {}

  // pass-through accessors for subspace labels
  HalfInt S() const { return BaseDegenerateStateType::subspace().S(); }

  int gamma_max() const
  {
    return BaseDegenerateStateType::subspace().degeneracy(index);
  }
};

template<typename tLGIType>
class LGISpace
    : public basis::BaseSpace<
          /* derived type */ LGISpace<tLGIType>,
          /* subspace type */ SpinSubspace<tLGIType>,
          /* space labels */ std::tuple<u3::U3>
        >
{
 private:
  using BaseSpaceType =
      basis::BaseSpace<LGISpace<tLGIType>, SpinSubspace<tLGIType>, std::tuple<u3::U3>>;
  using UpstreamLabelsType = typename tLGIType::UpstreamLabelsType;

 public:
  using LGIType = tLGIType;
  using SubspaceType = typename BaseSpaceType::SubspaceType;
  using SpaceLabelsType = typename BaseSpaceType::SpaceLabelsType;

  LGISpace() = default;

  LGISpace(
      const u3::U3& sigma,
      const std::map<
          typename SubspaceType::SubspaceLabelsType,
          typename MultiplicityTagged<UpstreamLabelsType>::vector
        >& spin_map
    )
      : BaseSpaceType{sigma}
  {
    for (const auto& [S, spin_vector] : spin_map)
      BaseSpaceType::EmplaceSubspace(S, spin_vector);
  }

  u3::U3 sigma() const { return std::get<0>(BaseSpaceType::labels()); }

  std::string LabelStr() const {return sigma().Str();}

};

template<typename tLGIType>
class Space
    : public basis::BaseSpace<
          /* derived type */ Space<tLGIType>,
          /* subspace type */ LGISpace<tLGIType>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<Space<tLGIType>, LGISpace<tLGIType>>;
  using UpstreamLabelsType = typename tLGIType::UpstreamLabelsType;

 public:
  using LGIType = tLGIType;
  using SubspaceType = typename BaseSpaceType::SubspaceType;

 public:
  Space() = default;

  Space(const typename MultiplicityTagged<tLGIType>::vector& lgi_vector, int Nmax)
  {
    std::map<u3::U3, std::map<std::tuple<HalfInt>, typename MultiplicityTagged<UpstreamLabelsType>::vector>>
        sigma_upstreams_map;
    for (const auto& [lgi, gamma_max] : lgi_vector)
    {
      const auto& upstream_labels = lgi.upstream_labels();
      if (lgi.Nex() > Nmax)
        continue;

      if (lgi.Nex() == 0)
        Nsigma0_ = lgi.sigma().N();

      sigma_upstreams_map[lgi.sigma()][lgi.S()].emplace_back(
          upstream_labels, gamma_max
        );
    }

    for (const auto& [sigma, spin_map] : sigma_upstreams_map)
      BaseSpaceType::EmplaceSubspace(sigma, spin_map);
  }

  HalfInt Nsigma0() const { return Nsigma0_; }

  HalfInt Nsigma0_;
};

////////////////////////////////////////////////////////////////////////////////
// spin::RecurrenceSpace() []
//   spin::RecurrenceLGISpace() [sigma,sigma',exchange_symm_bar]
//   ->spin::RecurrenceSpinSpace() [S,S']
//     -> spin::RecurrenceSpinSubspace() [Sp,Sn,Sp',Sn']/[T,T']->(gamma,gamma')
//       ->spin::RecurrenceOperatorState() [S0,T0,Sbar,Sbar',Tbar,Tbar']
////////////////////////////////////////////////////////////////////////////////

template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceSpace;
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceLGISpace;
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceSpinSpace;
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceSpinSubspace;
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceOperatorState;

// spin::RecurrenceS[PN/T]Subspace() [Sp,Sn,Sp',Sn']/[T,T']->(gamma,gamma')
template<typename tLGIType, typename tUnitTensorLabelsType>

class RecurrenceSpinSubspace
    : public basis::BaseSubspace<
          RecurrenceSpinSubspace<tLGIType, tUnitTensorLabelsType>, /* derived type */
          std::tuple<typename tLGIType::UpstreamLabelsType, typename tLGIType::UpstreamLabelsType>, /* subspace labels */
          RecurrenceOperatorState<tLGIType, tUnitTensorLabelsType>, /* state type */
          tUnitTensorLabelsType /* state labels */
        >
{
 private:
  using BaseSubspaceType = basis::BaseSubspace<
      RecurrenceSpinSubspace<tLGIType, tUnitTensorLabelsType>,
      std::tuple<typename tLGIType::UpstreamLabelsType, typename tLGIType::UpstreamLabelsType>,
      RecurrenceOperatorState<tLGIType, tUnitTensorLabelsType>,
      tUnitTensorLabelsType
    >;
  using UpstreamLabelsType = typename tLGIType::UpstreamLabelsType;

 public:
  using SubspaceLabelsType = typename BaseSubspaceType::SubspaceLabelsType;
  using StateLabelsType = typename BaseSubspaceType::StateLabelsType;
  using UnitTensorLabelsType = tUnitTensorLabelsType;

  RecurrenceSpinSubspace() = default;

  RecurrenceSpinSubspace(
      const SpinState<tLGIType>& spin_state_ket,
      HalfInt S_ket,
      const SpinState<tLGIType>& spin_state_bra,
      HalfInt S_bra,
      uint8_t exchange_symm_bar,
      uint8_t exchange_symm_barp
    );

  UpstreamLabelsType ket_upstream_labels() const
  {
    return std::get<0>(BaseSubspaceType::labels());
  }
  UpstreamLabelsType bra_upstream_labels() const
  {
    return std::get<1>(BaseSubspaceType::labels());
  }
};

template<typename tLGIType, typename tUnitTensorLabelsType>
RecurrenceSpinSubspace<tLGIType, tUnitTensorLabelsType>::RecurrenceSpinSubspace(
    const SpinState<tLGIType>& spin_state_ket,
    HalfInt S_ket,
    const SpinState<tLGIType>& spin_state_bra,
    HalfInt S_bra,
    uint8_t exchange_symm_bar,
    uint8_t exchange_symm_barp
  )
    : BaseSubspaceType{{spin_state_ket.labels(), spin_state_bra.labels()}}
{
  for (const auto tensor_labels : UnitTensorLabelsType::ALLOWED_LABELS)
  {
    if (!am::AllowedTriangle(S_ket, tensor_labels.S0(), S_bra))
      continue;
    if ((tensor_labels.exchange_symm_bar() != exchange_symm_bar) || (tensor_labels.exchange_symm_barp() != exchange_symm_barp))
      continue;
    if (!UnitTensorAllowed(
            spin_state_ket.labels(), tensor_labels, spin_state_bra.labels()
          ))
      continue;
    BaseSubspaceType::PushStateLabels(tensor_labels);
  }
}


// spin::RecurrenceOperatorState() [S0,T0,Sbar,Sbar',Tbar,Tbar']
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceOperatorState
    : public basis::BaseState<RecurrenceSpinSubspace<tLGIType, tUnitTensorLabelsType>>
{
 private:
  using BaseStateType =
      basis::BaseState<RecurrenceSpinSubspace<tLGIType, tUnitTensorLabelsType>>;

 public:
  using SubspaceType = typename BaseStateType::SubspaceType;
  using StateLabelsType = typename BaseStateType::StateLabelsType;

  // pass-through constructors

  RecurrenceOperatorState(const SubspaceType& subspace, std::size_t index)
  // Construct state by index.
      : BaseStateType{subspace, index}
  {}

  RecurrenceOperatorState(
      const SubspaceType& subspace, const StateLabelsType& state_labels
    )
  // Construct state by reverse lookup on labels.
      : BaseStateType{subspace, state_labels}
  {}

  // pass-through accessors
  //   Unit tensor properties should be accessed through the labels() method:
  //
  //     RecurrenceOperatorState state = subspace.GetState(i);
  //     std::cout << state.labels.S0() << "\n";
};


// spin::RecurrenceSpinSpace() [S,S']
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceSpinSpace
    : public basis::BaseDegenerateSpace<
          RecurrenceSpinSpace<tLGIType, tUnitTensorLabelsType>,
          RecurrenceSpinSubspace<tLGIType, tUnitTensorLabelsType>,
          std::tuple<HalfInt, HalfInt>
        >
{
 private:
  using BaseDegenerateSpaceType = basis::BaseDegenerateSpace<
      RecurrenceSpinSpace<tLGIType, tUnitTensorLabelsType>,
      RecurrenceSpinSubspace<tLGIType, tUnitTensorLabelsType>,
      std::tuple<HalfInt, HalfInt>
    >;
  using UpstreamLabelsType = typename tLGIType::UpstreamLabelsType;

 public:
  using SpaceLabelsType = typename BaseDegenerateSpaceType::SpaceLabelsType;

  RecurrenceSpinSpace() = default;

  RecurrenceSpinSpace(
      const SpinSubspace<tLGIType>& spin_subspace_ket,
      const SpinSubspace<tLGIType>& spin_subspace_bra,
      uint8_t exchange_symm_bar,
      uint8_t exchange_symm_barp
    );

  HalfInt S_ket() const
  {
    return std::get<0>(BaseDegenerateSpaceType::labels());
  }
  HalfInt S_bra() const
  {
    return std::get<1>(BaseDegenerateSpaceType::labels());
  }

  int GetKetSubspaceDegeneracy(std::size_t i) const
  {
    return std::get<0>(degeneracy_tuples_[i]);
  }
  int GetBraSubspaceDegeneracy(std::size_t i) const
  {
    return std::get<1>(degeneracy_tuples_[i]);
  }

  std::size_t GetSubspaceOffset(std::size_t i, int gamma_ket, int gamma_bra) const
  {
    const std::size_t degeneracy_index =
        (gamma_ket - 1) * GetBraSubspaceDegeneracy(i) + gamma_bra;
    return BaseDegenerateSpaceType::GetSubspaceOffset(i, degeneracy_index);
  }

  std::string LabelStr() const {return fmt::format("{} {}",S_ket(),S_bra());}

 private:
  std::vector<std::tuple<int, int>> degeneracy_tuples_;
};

template<typename tLGIType, typename tUnitTensorLabelsType>
RecurrenceSpinSpace<tLGIType, tUnitTensorLabelsType>::RecurrenceSpinSpace(
    const SpinSubspace<tLGIType>& spin_subspace_ket,
    const SpinSubspace<tLGIType>& spin_subspace_bra,
    uint8_t exchange_symm_bar,
    uint8_t exchange_symm_barp
  )
    : BaseDegenerateSpaceType{{spin_subspace_ket.S(), spin_subspace_bra.S()}}
{
  for (std::size_t i_ket = 0; i_ket < spin_subspace_ket.size(); ++i_ket)
  {
    for (std::size_t i_bra = 0; i_bra < spin_subspace_bra.size(); ++i_bra)
    {
      const auto& spin_state_ket = spin_subspace_ket.GetState(i_ket);
      const auto& spin_state_bra = spin_subspace_bra.GetState(i_bra);
      if (!UpstreamLabelsAllowed(spin_state_ket.labels(), spin_state_bra.labels()))
        continue;
      int degeneracy = spin_state_ket.degeneracy() * spin_state_bra.degeneracy();
      typename BaseDegenerateSpaceType::SubspaceType subspace(
          spin_state_ket,
          spin_subspace_ket.S(),
          spin_state_bra,
          spin_subspace_bra.S(),
          exchange_symm_bar,
          exchange_symm_barp
        );
      if (subspace.dimension() > 0)
      {
        BaseDegenerateSpaceType::PushSubspace(std::move(subspace), degeneracy);
        degeneracy_tuples_.push_back(
            {spin_state_ket.degeneracy(), spin_state_bra.degeneracy()}
          );
      }
    }
  }
}


// spin::RecurrenceLGISpace() [sigma,sigma',exchange_symm_bar]
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceLGISpace
    : public basis::BaseSpace<
          RecurrenceLGISpace<tLGIType, tUnitTensorLabelsType>,
          RecurrenceSpinSpace<tLGIType, tUnitTensorLabelsType>,
          std::tuple<u3::U3, u3::U3, uint8_t>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<
      RecurrenceLGISpace<tLGIType, tUnitTensorLabelsType>,
      RecurrenceSpinSpace<tLGIType, tUnitTensorLabelsType>,
      std::tuple<u3::U3, u3::U3, uint8_t>
    >;
  using UpstreamLabelsType = typename tLGIType::UpstreamLabelsType;

 public:
  using SpaceLabelsType = typename BaseSpaceType::SpaceLabelsType;

  RecurrenceLGISpace() = default;

  RecurrenceLGISpace(
      const LGISpace<tLGIType>& lgi_space_ket,
      const LGISpace<tLGIType>& lgi_space_bra,
      uint8_t exchange_symm_bar
    );

  u3::U3 sigma_ket() const { return std::get<0>(BaseSpaceType::labels()); }
  u3::U3 sigma_bra() const { return std::get<1>(BaseSpaceType::labels()); }
  uint8_t exchange_symm_bar() const { return std::get<2>(BaseSpaceType::labels()); }
  uint8_t exchange_symm_barp() const
  {
    return uint8_t(exchange_symm_bar() + abs(sigma_bra().N() - sigma_ket().N())) % 2;
  }

  std::string LabelStr() const 
  {
    return fmt::format("{} {} {}",sigma_ket().Str(),sigma_bra().Str(),exchange_symm_bar());
  }
};

template<typename tLGIType, typename tUnitTensorLabelsType>
RecurrenceLGISpace<tLGIType, tUnitTensorLabelsType>::RecurrenceLGISpace(
    const LGISpace<tLGIType>& lgi_space_ket,
    const LGISpace<tLGIType>& lgi_space_bra,
    uint8_t exchange_symm_bar
  )
    : BaseSpaceType{{lgi_space_ket.sigma(), lgi_space_bra.sigma(), exchange_symm_bar}}
{
  for (std::size_t i_ket = 0; i_ket < lgi_space_ket.size(); ++i_ket)
  {
    for (std::size_t i_bra = 0; i_bra < lgi_space_bra.size(); ++i_bra)
    {
      const auto& spin_subspace_ket = lgi_space_ket.GetSubspace(i_ket);
      const auto& spin_subspace_bra = lgi_space_bra.GetSubspace(i_bra);

      // two-body operator cannot change spin by more than 2; short circuit
      if (abs(spin_subspace_ket.S() - spin_subspace_bra.S()) > 2)
        continue;

      typename BaseSpaceType::SubspaceType subspace(
          spin_subspace_ket, spin_subspace_bra, exchange_symm_bar, exchange_symm_barp()
        );

      if (subspace.dimension() > 0)
        BaseSpaceType::PushSubspace(std::move(subspace));
    }
  }
}

// spin::RecurrenceSpace() []
template<typename tLGIType, typename tUnitTensorLabelsType>
class RecurrenceSpace
    : public basis::BaseSpace<
          RecurrenceSpace<tLGIType, tUnitTensorLabelsType>,
          RecurrenceLGISpace<tLGIType, tUnitTensorLabelsType>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<
      RecurrenceSpace<tLGIType, tUnitTensorLabelsType>,
      RecurrenceLGISpace<tLGIType, tUnitTensorLabelsType>
    >;
  using UpstreamLabelsType = typename tLGIType::UpstreamLabelsType;

 public:
  RecurrenceSpace() = default;

  RecurrenceSpace(const Space<tLGIType>& space_ket, const Space<tLGIType>& space_bra);
};

template<typename tLGIType, typename tUnitTensorLabelsType>
RecurrenceSpace<tLGIType, tUnitTensorLabelsType>::RecurrenceSpace(
    const Space<tLGIType>& space_ket, const Space<tLGIType>& space_bra
  )
    : BaseSpaceType{}
{
  for (std::size_t i_ket = 0; i_ket < space_ket.size(); ++i_ket)
    for (std::size_t i_bra = 0; i_bra < space_bra.size(); ++i_bra)
    {
      const auto& lgi_space_ket = space_ket.GetSubspace(i_ket);
      const auto& lgi_space_bra = space_bra.GetSubspace(i_bra);
      for (uint8_t exchange_symm_bar : {0, 1})
      {
        typename BaseSpaceType::SubspaceType subspace(
            lgi_space_ket, lgi_space_bra, exchange_symm_bar
          );

        if (subspace.dimension() > 0)
          BaseSpaceType::PushSubspace(std::move(subspace));
      }
    }
}

}  // namespace spin
}  // namespace spncci

#endif  // RECURRENCE_INDEXING_SPIN_H_
