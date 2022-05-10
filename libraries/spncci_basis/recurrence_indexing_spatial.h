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
  + 05/10/22 (aem): Templatized spatial space.
****************************************************************/
#ifndef RECURRENCE_INDEXING_SPATIAL_H_
#define RECURRENCE_INDEXING_SPATIAL_H_

#include <cppitertools/enumerate.hpp>  // iter
#include <memory>                      // for std::shared_ptr
#include <tuple>
#include <utility>  // for std::forward

#include "am/halfint.h"
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "fmt/format.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"
#include "spncci_basis/recurrence_indexing_spin.h"
#include "u3shell/operator_indexing_spatial.h"


namespace spncci::spatial
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// spncci::spatial::Space() []
// ->spncci::spatial::Sp3RSpace() [sigma]
//     ->spncci::spatial::U3Subspace() [omega]
/////////////////////////////////////////////////////////////////////////////////////////////////////////

using Sp3RSpace = sp3r::Sp3RSpace;
using U3Subspace = sp3r::U3Subspace;

// class Space;

class Space
    : public basis::BaseSpace<Space, Sp3RSpace>
{
 public:
  Space() = default;

  // Construct from list of sigma
  Space(
      const std::vector<u3::U3>& sigma_vector,
      const HalfInt& Nsigma0,
      const unsigned int Nmax
    );

  // Construct from spin::Space
  template<typename tLGIType>
  Space(
      const spin::Space<tLGIType>& spin_space,
      const HalfInt& Nsigma0,
      const unsigned int& Nmax
    )
  {
    for (int i = 0; i < spin_space.size(); ++i)
    {
      const u3::U3& sigma = spin_space.GetSubspace(i).sigma();
      // Only those sigma for which Nmax is >= sigma.N()-Nsigma0
      // are included in spin space.  Thus Nn_max is always >= 0.
      unsigned int Nn_max = Nmax - static_cast<int>(sigma.N() - Nsigma0);
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
template<typename tOperatorStateLabelType> class RecurrenceSpace;
template<typename tOperatorStateLabelType> class RecurrenceSp3RSpace;
template<typename tOperatorStateLabelType> class RecurrenceNnsumSpace;
template<typename tOperatorStateLabelType> class RecurrenceU3Space;

struct OperatorConstraints
{
  OperatorConstraints(unsigned int N1v_, HalfInt Nsigma0_, uint8_t parity_bar_)
      : N1v{N1v_}, Nsigma0{Nsigma0_}, parity_bar{parity_bar_}
  {}

  unsigned int N1v;
  HalfInt Nsigma0;
  uint8_t parity_bar;
};

////////////////////////////////////////////////////////////////////////////////////////
// spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
////////////////////////////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class RecurrenceU3Space
    : public basis::BaseDegenerateSpace<
          RecurrenceU3Space<tOperatorStateLabelType>,
          u3shell::spatial::OperatorU3Subspace<tOperatorStateLabelType>,
          std::tuple<u3::U3, u3::U3>
        >
{
 private:
  using OperatorStateLabelType = tOperatorStateLabelType;
  using BaseDegenerateSpaceType = basis::BaseDegenerateSpace<
      RecurrenceU3Space<OperatorStateLabelType>,
      u3shell::spatial::OperatorU3Subspace<OperatorStateLabelType>,
      std::tuple<u3::U3, u3::U3>
    >;


 public:
  RecurrenceU3Space() = default;

  RecurrenceU3Space(
      const std::tuple<u3::U3, u3::U3>& omega_pair,
      const spncci::spatial::OperatorConstraints& operator_constraints
    )
      : BaseDegenerateSpaceType{omega_pair}
  {
    const auto& [omega, omegap] = omega_pair;

    // Determine limits on Nbar and Nbarp
    unsigned int Nbar_max =
        2 * operator_constraints.N1v
        + std::max(
            static_cast<unsigned int>(omega.N() - operator_constraints.Nsigma0),
            static_cast<unsigned int>(omegap.N() - operator_constraints.Nsigma0)
          );

    int N0 = static_cast<int>(omegap.N() - omega.N());

    //// GenerateSpatialOperators must be defined by type of operator in correct namespace?
    std::map<u3::SU3, std::vector<OperatorStateLabelType>> spatial_operators;
    u3shell::spatial::GenerateSpatialOperators(
        N0, operator_constraints.parity_bar, Nbar_max, spatial_operators
      );

    for (const auto& [x0, state_vector] : spatial_operators)
    {
      unsigned int rho_max = u3::OuterMultiplicity(omega, x0, omegap);
      auto subspace =
          u3shell::spatial::OperatorU3Subspace<OperatorStateLabelType>();
      BaseDegenerateSpaceType::PushSubspace(std::move(subspace), rho_max);
    }
  }
  const u3::U3& omega_ket() const
  {
    return std::get<0>(BaseDegenerateSpaceType::labels());
  }
  const u3::U3& omega_bra() const
  {
    return std::get<1>(BaseDegenerateSpaceType::labels());
  }
  std::string LabelStr() const
  {
    return fmt::format("{} {}", omega_ket().Str(), omega_bra().Str());
  }
};
////////////////////////////////////////////////////////////////////////////////////////
// spatial::RecurrenceNnsumSpace() [Nsum]
////////////////////////////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class RecurrenceNnsumSpace
    : public basis::BaseDegenerateSpace<
          RecurrenceNnsumSpace<tOperatorStateLabelType>,
          RecurrenceU3Space<tOperatorStateLabelType>,
          std::tuple<unsigned int>
        >
{
 private:
  using BaseDegenerateSpaceType = basis::BaseDegenerateSpace<
      RecurrenceNnsumSpace<tOperatorStateLabelType>,
      RecurrenceU3Space<tOperatorStateLabelType>,
      std::tuple<unsigned int>
    >;

 public:
  RecurrenceNnsumSpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceNnsumSpace(
      unsigned int Nnsum,
      const std::vector<std::tuple<std::size_t, std::size_t>> u3subspace_index_pairs,
      const sp3r::Sp3RSpace& sp3r_space_ket,
      const sp3r::Sp3RSpace& sp3r_space_bra,
      const OperatorConstraints& operator_constraints
    )
      : BaseDegenerateSpaceType{Nnsum}
  {
    parity_bar_ = operator_constraints.parity_bar;

    BaseDegenerateSpaceType::reserve(u3subspace_index_pairs.size());
    upsilon_max_pairs_.reserve(u3subspace_index_pairs.size());

    for (const auto& [i_ket, i_bra] : u3subspace_index_pairs)
    {
      const auto& u3subspace_ket = sp3r_space_ket.GetSubspace(i_ket);
      const auto& u3subspace_bra = sp3r_space_bra.GetSubspace(i_bra);
      const u3::U3& omega_ket = u3subspace_ket.omega();
      const u3::U3& omega_bra = u3subspace_bra.omega();

      const unsigned int& upsilon_max_ket = u3subspace_ket.size();
      const unsigned int& upsilon_max_bra = u3subspace_bra.size();

      auto subspace = RecurrenceU3Space<tOperatorStateLabelType>(
          {omega_ket, omega_bra}, operator_constraints
        );

      if (subspace.dimension() == 0)
        continue;

      upsilon_max_pairs_.push_back({upsilon_max_ket, upsilon_max_bra});
      BaseDegenerateSpaceType::PushSubspace(
          std::move(subspace), upsilon_max_bra * upsilon_max_ket
        );
    }
  }

  unsigned int Nnsum() const
  {
    return std::get<0>(BaseDegenerateSpaceType::labels());
  }

  unsigned int upsilon_max_ket(const std::size_t i) const
  {
    return std::get<0>(upsilon_max_pairs_[i]);
  }
  unsigned int upsilon_max_bra(const std::size_t i) const
  {
    return std::get<1>(upsilon_max_pairs_[i]);
  }

  std::size_t GetSubspaceOffset(
      std::size_t i, unsigned int upsilon_ket = 1, unsigned int upsilon_bra = 1
    ) const
  {
    assert(upsilon_ket > 0 && upsilon_bra > 0);
    const std::size_t degeneracy_index =
        (upsilon_ket - 1) * upsilon_max_bra(i) + (upsilon_bra - 1);
    return BaseDegenerateSpaceType::GetSubspaceOffset(i, degeneracy_index);
  }

  uint8_t parity_bar() const { return parity_bar_; }

  inline std::string LabelStr() const { return fmt::format("{}", Nnsum()); }

 private:
  uint8_t parity_bar_;
  std::vector<std::tuple<unsigned int, unsigned int>> upsilon_max_pairs_;
};


////////////////////////////////////////////////////////////////////////////////////////
// spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
////////////////////////////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class RecurrenceSp3RSpace
    : public basis::BaseSpace<
          RecurrenceSp3RSpace<tOperatorStateLabelType>,
          RecurrenceNnsumSpace<tOperatorStateLabelType>,
          std::tuple<u3::U3, u3::U3, uint8_t>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<
      RecurrenceSp3RSpace<tOperatorStateLabelType>,
      RecurrenceNnsumSpace<tOperatorStateLabelType>,
      std::tuple<u3::U3, u3::U3, uint8_t>
    >;

 public:
  RecurrenceSp3RSpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceSp3RSpace(
      std::shared_ptr<const Sp3RSpace> sp3r_space_ket_ptr,
      std::shared_ptr<const Sp3RSpace> sp3r_space_bra_ptr,
      const OperatorConstraints& operator_constraints
    )
    : BaseSpaceType{
        {sp3r_space_ket_ptr->sigma(),
         sp3r_space_bra_ptr->sigma(),
         operator_constraints.parity_bar}
      },
      sp3r_space_ket_ptr_{std::move(sp3r_space_ket_ptr)},
      sp3r_space_bra_ptr_{std::move(sp3r_space_bra_ptr)}
  {
    const u3::U3& sigma_ket = ket_space().sigma();
    const u3::U3& sigma_bra = bra_space().sigma();

    // Partition pairs of omega',omega by Nnsum
    std::map<unsigned int, std::vector<std::tuple<std::size_t, std::size_t>>>
        Nnsum_partition;
    for (std::size_t i_ket = 0; i_ket < ket_space().size(); ++i_ket)
      for (std::size_t i_bra = 0; i_bra < bra_space().size(); ++i_bra)
      {
        const u3::U3& omega_ket = ket_space().GetSubspace(i_ket).omega();
        const u3::U3& omega_bra = bra_space().GetSubspace(i_bra).omega();

        unsigned int Nnsum = static_cast<unsigned int>(
            omega_ket.N() - sigma_ket.N() + omega_bra.N() - sigma_bra.N()
          );
        Nnsum_partition[Nnsum].push_back({i_ket, i_bra});
      }

    // Create RecurrenceNnsumSpaces.  On for each unit tensor state parity
    BaseSpaceType::reserve(Nnsum_partition.size());
    for (const auto& [Nnsum, partition] : Nnsum_partition)
    {
      auto subspace = RecurrenceNnsumSpace<tOperatorStateLabelType>(
          Nnsum, partition, ket_space(), bra_space(), operator_constraints
        );
      if ((Nnsum == 0) && (subspace.dimension() == 0))
        return;  // TODO: is the Nnsum constraint necessary, or is this more general?
      if (subspace.dimension() == 0)
        continue;
      BaseSpaceType::PushSubspace(std::move(subspace));
    }
  }


  const u3::U3& sigma_ket() const
  {
    return std::get<0>(BaseSpaceType::labels());
  }
  const u3::U3& sigma_bra() const
  {
    return std::get<1>(BaseSpaceType::labels());
  }
  uint8_t parity_bar() const { return std::get<2>(BaseSpaceType::labels()); }
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


////////////////////////////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class RecurrenceSpace
    : public basis::BaseSpace<
          RecurrenceSpace<tOperatorStateLabelType>,
          RecurrenceSp3RSpace<tOperatorStateLabelType>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<
      RecurrenceSpace<tOperatorStateLabelType>,
      RecurrenceSp3RSpace<tOperatorStateLabelType>
    >;

 public:
  RecurrenceSpace() = default;

  // spatial_unit_tensors <(x0,Nbar_p,Nbar)>
  RecurrenceSpace(
      const spatial::Space& space_ket,
      const spatial::Space& space_bra,
      const unsigned int& N1v,
      const HalfInt& Nsigma0
    )
      : BaseSpaceType{}
  {
    BaseSpaceType::reserve(2 * space_ket.size() * space_bra.size());
    for (std::size_t i_ket = 0; i_ket < space_ket.size(); ++i_ket)
      for (std::size_t i_bra = 0; i_bra < space_bra.size(); ++i_bra)
      {
        const auto& sp3r_space_ket_ptr = space_ket.GetSubspacePtr(i_ket);
        const auto& sp3r_space_bra_ptr = space_bra.GetSubspacePtr(i_bra);
        for (uint8_t parity_bar : {0, 1})
        {
          auto subspace = RecurrenceSp3RSpace<tOperatorStateLabelType>(
              sp3r_space_ket_ptr, sp3r_space_bra_ptr, {N1v, Nsigma0, parity_bar}
            );
          if (subspace.dimension() == 0)
            continue;

          BaseSpaceType::PushSubspace(std::move(subspace));
        }
      }
  }

  // private:
};


////////////////////////////////////////////////////////////////
// recurrence sectors
////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class RecurrenceU3Sector
    : public basis::BaseSector<RecurrenceU3Space<tOperatorStateLabelType>>
{
 private:
  using OperatorStateLabelType = tOperatorStateLabelType;

 public:
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////
  using BaseSectorType =
      basis::BaseSector<RecurrenceU3Space<OperatorStateLabelType>>;

  using BaseSectorType::BaseSector;  // Querry Patrick: What does this do?

  std::size_t source_subspace_index() const
  {
    return BaseSectorType::ket_subspace_index();
  }
  std::size_t target_subspace_index() const
  {
    return BaseSectorType::bra_subspace_index();
  }
  const auto& source_subspace() const { return BaseSectorType::ket_subspace(); }
  const auto& target_subspace() const { return BaseSectorType::bra_subspace(); }
};


////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class RecurrenceU3Sectors
    : public basis::BaseSectors<
          RecurrenceNnsumSpace<tOperatorStateLabelType>,
          RecurrenceNnsumSpace<tOperatorStateLabelType>,
          RecurrenceU3Sector<tOperatorStateLabelType>
        >
{
 private:
  using OperatorStateLabelType = tOperatorStateLabelType;
  using BaseSectorsType = basis::BaseSectors<
      RecurrenceNnsumSpace<OperatorStateLabelType>,
      RecurrenceNnsumSpace<OperatorStateLabelType>,
      RecurrenceU3Sector<OperatorStateLabelType>
    >;

 public:
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  RecurrenceU3Sectors() = default;

  RecurrenceU3Sectors(
      const RecurrenceSp3RSpace<OperatorStateLabelType>& sp3r_space,
      int target_Nnsum,
      int source_Nnsum
    )
      : BaseSectorsType{
          sp3r_space.LookUpSubspace({target_Nnsum}),
          sp3r_space.LookUpSubspace({source_Nnsum})
        }
  {
    int delta_Nnsum = target_Nnsum - source_Nnsum;
    assert((delta_Nnsum == 2) || (delta_Nnsum == 4));
    const auto& target_Nnsum_space = BaseSectorsType::bra_space();
    const auto& source_Nnsum_space = BaseSectorsType::ket_space();

    if (delta_Nnsum == 2)
    {
      for (auto&& [target_u3_index, target_u3_space] :
           iter::enumerate(target_Nnsum_space))
      {
        for (auto&& [source_u3_index, source_u3_space] :
             iter::enumerate(source_Nnsum_space))
        {
          if (source_u3_space.omega_ket() == sp3r_space.sigma_ket())
          {
            if (target_u3_space.omega_ket() != source_u3_space.omega_ket())
              continue;

            if (
                u3::OuterMultiplicity(
                    {2u, 0u}, source_u3_space.omega_bra(), target_u3_space.omega_bra()
                  )
                == 0
              )
              continue;
          }
          else
          {
            if (target_u3_space.omega_bra() != source_u3_space.omega_bra())
              continue;

            if (
                u3::OuterMultiplicity(
                    {2u, 0u}, source_u3_space.omega_ket(), target_u3_space.omega_ket()
                  )
                == 0
              )
              continue;
          }
          BaseSectorsType::PushSector(target_u3_index, source_u3_index);
        }
      }
    }
    else  // if (delta_Nnsum == 4)
    {
      for (auto&& [target_u3_index, target_u3_space] :
           iter::enumerate(target_Nnsum_space))
      {
        for (auto&& [source_u3_index, source_u3_space] :
             iter::enumerate(source_Nnsum_space))
        {
          if (
              u3::OuterMultiplicity(
                  {2u, 0u}, source_u3_space.omega_bra(), target_u3_space.omega_bra()
                )
              == 0
            )
            continue;

          if (
              u3::OuterMultiplicity(
                  {2u, 0u}, source_u3_space.omega_ket(), target_u3_space.omega_ket()
                )
              == 0
            )
            continue;

          BaseSectorsType::PushSector(target_u3_index, source_u3_index);
        }
      }
    }
  }
};

}  // namespace spncci::spatial

#endif  // RECURRENCE_INDEXING_H_
