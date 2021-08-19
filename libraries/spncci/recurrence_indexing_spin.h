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
      spin::RecurrenceLGISpace() [sigma,sigma',parity_bar]
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
#include "lgi/lgi.h"
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

  HalfInt S() const
  {
    return std::get<0>(BaseDegenerateSubspaceType::labels());
  }

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
//   spin::RecurrenceLGISpace() [sigma,sigma']
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
      uint8_t parity_bar,
      uint8_t parity_barp
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
    uint8_t parity_bar,
    uint8_t parity_barp
  )
    : BaseSubspaceType{{spin_state_ket.labels(), spin_state_bra.labels()}}
{
  for (const auto tensor_labels : UnitTensorLabelsType::ALLOWED_LABELS)
  {
    if (!am::AllowedTriangle(S_ket, tensor_labels.S0(), S_bra))
      continue;
    if ((tensor_labels.parity_bar() != parity_bar) || (tensor_labels.parity_barp() != parity_barp))
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

  RecurrenceOperatorState(const SubspaceType& subspace, const StateLabelsType& state_labels)
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
      uint8_t parity_bar,
      uint8_t parity_barp
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

 private:
  std::vector<std::tuple<int, int>> degeneracy_tuples_;
};

template<typename tLGIType, typename tUnitTensorLabelsType>
RecurrenceSpinSpace<tLGIType, tUnitTensorLabelsType>::RecurrenceSpinSpace(
    const SpinSubspace<tLGIType>& spin_subspace_ket,
    const SpinSubspace<tLGIType>& spin_subspace_bra,
    uint8_t parity_bar,
    uint8_t parity_barp
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
          parity_bar,
          parity_barp
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


// spin::RecurrenceLGISpace() [sigma,sigma',parity_bar]
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
      uint8_t parity_bar
    );

  u3::U3 sigma_ket() const { return std::get<0>(BaseSpaceType::labels()); }
  u3::U3 sigma_bra() const { return std::get<1>(BaseSpaceType::labels()); }
  uint8_t parity_bar() const { return std::get<2>(BaseSpaceType::labels()); }
  uint8_t parity_barp() const
  {
    return uint8_t(parity_bar() + abs(sigma_bra().N() - sigma_ket().N())) % 2;
  }
};

template<typename tLGIType, typename tUnitTensorLabelsType>
RecurrenceLGISpace<tLGIType, tUnitTensorLabelsType>::RecurrenceLGISpace(
    const LGISpace<tLGIType>& lgi_space_ket,
    const LGISpace<tLGIType>& lgi_space_bra,
    uint8_t parity_bar
  )
    : BaseSpaceType{{lgi_space_ket.sigma(), lgi_space_bra.sigma(), parity_bar}}
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
          spin_subspace_ket, spin_subspace_bra, parity_bar, parity_barp()
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
      for (uint8_t parity_bar : {0, 1})
      {
        typename BaseSpaceType::SubspaceType subspace(
            lgi_space_ket, lgi_space_bra, parity_bar
          );

        if (subspace.dimension() > 0)
          BaseSpaceType::PushSubspace(std::move(subspace));
      }
    }
}

}  // namespace spin
}  // namespace spncci

////////////////////////////////////////////////////////////////
// labels for good-ST unit tensors
////////////////////////////////////////////////////////////////
namespace spncci::spin
{
class UnitTensorLabelsST;
}
namespace std
{
template<> struct hash<spncci::spin::UnitTensorLabelsST>;
}
namespace spncci::spin
{
class UnitTensorLabelsST
{
 private:
  // conversions from std::bitset and integer
  constexpr explicit UnitTensorLabelsST(std::bitset<8> bitrep)
      : bitrep_{bitrep}
  {}

  constexpr explicit UnitTensorLabelsST(uint8_t i)
      : bitrep_{i}
  {}

  constexpr operator std::bitset<8>() const { return bitrep_; }
  friend struct std::hash<UnitTensorLabelsST>;

 public:
  unsigned long id() const { return bitrep_.to_ulong(); }

  constexpr UnitTensorLabelsST(int S0, int T0, int Sbar, int Sbarp, int Tbar, int Tbarp)
      : bitrep_{static_cast<uint>(
          64 * S0 + 16 * T0 + 8 * Sbar + 4 * Sbarp + 2 * Tbar + Tbarp
        )}
  {}

  static std::array<UnitTensorLabelsST, 36> const ALLOWED_LABELS;

  inline bool Allowed() const
  {
    if (Sbar() < 0 || Sbar() > 1)
      return false;
    if (Sbarp() < 0 || Sbarp() > 1)
      return false;
    if (!am::AllowedTriangle(Sbar(), S0(), Sbarp()))
      return false;
    if (Tbar() < 0 || Tbar() > 1)
      return false;
    if (Tbarp() < 0 || Tbarp() > 1)
      return false;
    if (!am::AllowedTriangle(Tbar(), T0(), Tbarp()))
      return false;
    return true;
  }

  inline unsigned int S0() const
  {
    // S0 is stored in bits 7 and 6.
    return ((bitrep_ & std::bitset<8>{0b11000000}) >> 6).to_ulong();
  }
  inline unsigned int T0() const
  {
    // T0 is stored in bits 5 and 4.
    return ((bitrep_ & std::bitset<8>{0b00110000}) >> 4).to_ulong();
  }
  inline unsigned int Sbar() const
  {
    // Sbar is stored in bit 3.
    return ((bitrep_ & std::bitset<8>{0b00001000}) >> 3).to_ulong();
  }
  inline unsigned int Sbarp() const
  {
    // Sbarp is stored in bit 2.
    return ((bitrep_ & std::bitset<8>{0b00000100}) >> 2).to_ulong();
  }
  inline unsigned int Tbar() const
  {
    // Tbar is stored in bit 1.
    return ((bitrep_ & std::bitset<8>{0b00000010}) >> 1).to_ulong();
  }
  inline unsigned int Tbarp() const
  {
    // Tbarp is stored in bit 0.
    return ((bitrep_ & std::bitset<8>{0b00000001}).to_ulong());
  }

  inline uint8_t parity_bar() const { return (Sbar() + Tbar()) % 2; }
  inline uint8_t parity_barp() const { return (Sbarp() + Tbarp()) % 2; }

  inline bool operator==(const UnitTensorLabelsST& rhs) const
  {
    return std::bitset<8>(*this) == std::bitset<8>(rhs);
  }

 private:
  std::bitset<8> bitrep_;
};

inline constexpr std::array<UnitTensorLabelsST, 36> const UnitTensorLabelsST::ALLOWED_LABELS{{
    // (S0, T0, Sbar, Sbarp, Tbar, Tbarp)
    UnitTensorLabelsST{0b00000000u},  // (0, 0, 0, 0, 0, 0)
    UnitTensorLabelsST{0b00000011u},  // (0, 0, 0, 0, 1, 1)
    UnitTensorLabelsST{0b00001100u},  // (0, 0, 1, 1, 0, 0)
    UnitTensorLabelsST{0b00001111u},  // (0, 0, 1, 1, 1, 1)
    UnitTensorLabelsST{0b00010001u},  // (0, 1, 0, 0, 0, 1)
    UnitTensorLabelsST{0b00010010u},  // (0, 1, 0, 0, 1, 0)
    UnitTensorLabelsST{0b00010011u},  // (0, 1, 0, 0, 1, 1)
    UnitTensorLabelsST{0b00011101u},  // (0, 1, 1, 1, 0, 1)
    UnitTensorLabelsST{0b00011110u},  // (0, 1, 1, 1, 1, 0)
    UnitTensorLabelsST{0b00011111u},  // (0, 1, 1, 1, 1, 1)
    UnitTensorLabelsST{0b00100011u},  // (0, 2, 0, 0, 1, 1)
    UnitTensorLabelsST{0b00101111u},  // (0, 2, 1, 1, 1, 1)
    UnitTensorLabelsST{0b01000100u},  // (1, 0, 0, 1, 0, 0)
    UnitTensorLabelsST{0b01000111u},  // (1, 0, 0, 1, 1, 1)
    UnitTensorLabelsST{0b01001000u},  // (1, 0, 1, 0, 0, 0)
    UnitTensorLabelsST{0b01001011u},  // (1, 0, 1, 0, 1, 1)
    UnitTensorLabelsST{0b01001100u},  // (1, 0, 1, 1, 0, 0)
    UnitTensorLabelsST{0b01001111u},  // (1, 0, 1, 1, 1, 1)
    UnitTensorLabelsST{0b01010101u},  // (1, 1, 0, 1, 0, 1)
    UnitTensorLabelsST{0b01010110u},  // (1, 1, 0, 1, 1, 0)
    UnitTensorLabelsST{0b01010111u},  // (1, 1, 0, 1, 1, 1)
    UnitTensorLabelsST{0b01011001u},  // (1, 1, 1, 0, 0, 1)
    UnitTensorLabelsST{0b01011010u},  // (1, 1, 1, 0, 1, 0)
    UnitTensorLabelsST{0b01011011u},  // (1, 1, 1, 0, 1, 1)
    UnitTensorLabelsST{0b01011101u},  // (1, 1, 1, 1, 0, 1)
    UnitTensorLabelsST{0b01011110u},  // (1, 1, 1, 1, 1, 0)
    UnitTensorLabelsST{0b01011111u},  // (1, 1, 1, 1, 1, 1)
    UnitTensorLabelsST{0b01100111u},  // (1, 2, 0, 1, 1, 1)
    UnitTensorLabelsST{0b01101011u},  // (1, 2, 1, 0, 1, 1)
    UnitTensorLabelsST{0b01101111u},  // (1, 2, 1, 1, 1, 1)
    UnitTensorLabelsST{0b10001100u},  // (2, 0, 1, 1, 0, 0)
    UnitTensorLabelsST{0b10001111u},  // (2, 0, 1, 1, 1, 1)
    UnitTensorLabelsST{0b10011101u},  // (2, 1, 1, 1, 0, 1)
    UnitTensorLabelsST{0b10011110u},  // (2, 1, 1, 1, 1, 0)
    UnitTensorLabelsST{0b10011111u},  // (2, 1, 1, 1, 1, 1)
    UnitTensorLabelsST{0b10101111u}   // (2, 2, 1, 1, 1, 1)
  }};

inline bool UnitTensorAllowed(
    const lgi::LGI::UpstreamLabelsType& ket_labels,
    UnitTensorLabelsST tensor_labels,
    const lgi::LGI::UpstreamLabelsType& bra_labels
  )
{
  if ((tensor_labels.S0() == 0) && (tensor_labels.T0() == 0))
  {
    if ((ket_labels.Sp != bra_labels.Sp) || (ket_labels.Sn != bra_labels.Sn))
      return false;
  }
  return true;
}

}  // namespace spncci::spin

namespace std
{
template<> struct hash<spncci::spin::UnitTensorLabelsST>
{
  inline std::size_t operator()(const spncci::spin::UnitTensorLabelsST& h) const noexcept
  {
    return std::hash<std::bitset<8>>()(std::bitset<8>(h));
  }
};
}  // namespace std

namespace spncci::spin
{
inline std::size_t hash_value(const UnitTensorLabelsST& l)
{
  return std::hash<UnitTensorLabelsST>{}(l);
}
}  // namespace spncci::spin


#endif  // RECURRENCE_INDEXING_SPIN_H_