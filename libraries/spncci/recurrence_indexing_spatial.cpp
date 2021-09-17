/****************************************************************
  recurrence_indexing.cpp

  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/recurrence_indexing_spatial.h"

#include <algorithm>
#include <fstream>
#include <iostream>

#include "basis/basis.h"
// #include "mcutils/parsing.h"
#include "fmt/format.h"
// #include "am/halfint_fmt.h"

// #include "sp3rlib/vcs.h"

namespace spncci
{
namespace spatial
{
U3Subspace::U3Subspace(const sp3r::U3Subspace& u3subspace)
    : BaseDegenerateSubspace{u3subspace.labels()}
{
  const u3::U3& omega = u3subspace.labels();
  std::map<u3::U3, std::size_t> n_map;
  for (int i = 0; i < u3subspace.size(); ++i)
  {
    const auto& n_rho = u3subspace.GetStateLabels(i);
    if (n_map.count(n_rho.irrep) == 0)
      n_map[n_rho.irrep] = 0;
    n_map[n_rho.irrep] += n_rho.tag;
  }
  for (const auto& [n, rho_max] : n_map)
    PushStateLabels(n, rho_max);
}

LGISpace::LGISpace(const u3::U3& sigma, const int Nn_max)
    : BaseSpace{sigma}
{
  sp3r::Sp3RSpace sp3r_space(sigma, Nn_max);
  for (int i = 0; i < sp3r_space.size(); i++)
  {
    const auto& u3subspace = sp3r_space.GetSubspace(i);
    PushSubspace(U3Subspace(u3subspace));
  }
}

Space::Space(const std::vector<u3::U3>& sigma_vector, const HalfInt& Nsigma0, const int Nmax)
    : BaseSpace{}
{
  for (const auto& sigma : sigma_vector)
  {
    // Nsigma0_=Nsigma0;
    int Nn_max = Nmax - int(sigma.N() - Nsigma0);
    if (Nn_max >= 0)
      PushSubspace(LGISpace(sigma, Nn_max));
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Recurrence indexing (spatial)
///////////////////////////////////////////////////////////////////////////////////////////////////////


RecurrenceOperatorSubspace::RecurrenceOperatorSubspace(
    const u3::SU3& x0, const std::vector<std::tuple<int, int>>& Nbar_pairs
  )
    : BaseSubspace{x0}
{
  reserve(Nbar_pairs.size());
  for (const auto& Nbar_pair : Nbar_pairs) PushStateLabels(Nbar_pair);
}

RecurrenceU3Space::RecurrenceU3Space(
    const std::tuple<u3::U3, u3::U3>& omega_pair,
    const UnitTensorConstraintParameters& unit_tensor_constraints
  )
    : BaseDegenerateSpace{omega_pair}
{
  const auto& [omega, omega_p] = omega_pair;
  std::unordered_map<u3::SU3, std::vector<std::tuple<int, int>>, boost::hash<u3::SU3>> x0_Nbar_pairs;

  ////////////////////////////////////////////////////////////////////////////////
  // Create list of spatial unit tensors to pass through to operator subspace
  int Nbar_max{
      2 * unit_tensor_constraints.N1v + omega.N() - unit_tensor_constraints.Nsigma0
    };

  int Nbar_p_max{
      2 * unit_tensor_constraints.N1v + omega_p.N() - unit_tensor_constraints.Nsigma0
    };

  int N0{omega_p.N() - omega.N()};
  int Nbar_min = unit_tensor_constraints.parity_bar;
  for (int Nbar = Nbar_min; Nbar <= Nbar_max; Nbar += 2)
  {
    int Nbar_p = N0 + Nbar;
    if (Nbar_p <= Nbar_p_max & Nbar_p >= 0)
    {
      MultiplicityTagged<u3::SU3>::vector possible_x0 =
          u3::KroneckerProduct({Nbar_p, 0}, {0, Nbar});

      for (const auto& [x0, rho] : possible_x0)
        x0_Nbar_pairs[x0].push_back({Nbar, Nbar_p});
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Construct operator subspaces
  ////////////////////////////////////////////////////////////////////////////////
  reserve(x0_Nbar_pairs.size());
  for (const auto& [x0, Nbar_pairs] : x0_Nbar_pairs)
  {
    int rho0_max = u3::OuterMultiplicity(omega_ket().SU3(), x0, omega_bra().SU3());
    if (rho0_max > 0)
    {
      auto subspace = RecurrenceOperatorSubspace(x0, Nbar_pairs);
      if (subspace.dimension() == 0)
        continue;
      PushSubspace(std::move(subspace), rho0_max);
    }
  }
}


RecurrenceNnsumSpace::RecurrenceNnsumSpace(
    int Nnsum,
    const std::vector<std::tuple<int, int>> u3subspace_index_pairs,
    const LGISpace& lgi_space_ket,
    const LGISpace& lgi_space_bra,
    const UnitTensorConstraintParameters& unit_tensor_constraints
  )
    : BaseDegenerateSpace{Nnsum}
{
  parity_bar_ = unit_tensor_constraints.parity_bar;

  reserve(u3subspace_index_pairs.size());
  upsilon_pairs_.reserve(u3subspace_index_pairs.size());

  for (const auto& [i_ket, i_bra] : u3subspace_index_pairs)
  {
    const auto& u3subspace_ket = lgi_space_ket.GetSubspace(i_ket);
    const auto& u3subspace_bra = lgi_space_bra.GetSubspace(i_bra);
    const u3::U3& omega_ket = u3subspace_ket.omega();
    const u3::U3& omega_bra = u3subspace_bra.omega();

    const int& upsilon_max_ket = u3subspace_ket.size();
    const int& upsilon_max_bra = u3subspace_bra.size();

    auto subspace = RecurrenceU3Space({omega_ket, omega_bra}, unit_tensor_constraints);
    if (subspace.dimension() == 0)
      continue;

    upsilon_pairs_.push_back({upsilon_max_ket, upsilon_max_bra});
    PushSubspace(std::move(subspace), upsilon_max_bra * upsilon_max_ket);
  }
}


RecurrenceLGISpace::RecurrenceLGISpace(
    const LGISpace& lgi_space_ket,
    const LGISpace& lgi_space_bra,
    const UnitTensorConstraintParameters& unit_tensor_constraints
  )
    : BaseSpace{
        {lgi_space_ket.sigma(), lgi_space_bra.sigma(), unit_tensor_constraints.parity_bar}
      }
{
  const u3::U3& sigma_ket = lgi_space_ket.sigma();
  const u3::U3& sigma_bra = lgi_space_bra.sigma();
  // std::cout<<sigma_bra.Str()<<"  "<<sigma_ket.Str()<<std::endl;

  // Partition pairs of omega',omega by Nnsum
  std::map<int, std::vector<std::tuple<int, int>>> Nnsum_partition;
  for (int i_ket = 0; i_ket < lgi_space_ket.size(); ++i_ket)
    for (int i_bra = 0; i_bra < lgi_space_bra.size(); ++i_bra)
    {
      const u3::U3& omega_ket = lgi_space_ket.GetSubspace(i_ket).omega();
      const u3::U3& omega_bra = lgi_space_bra.GetSubspace(i_bra).omega();

      int Nnsum =
          int(omega_ket.N() - sigma_ket.N() + omega_bra.N() - sigma_bra.N());
      Nnsum_partition[Nnsum].push_back({i_ket, i_bra});
    }

  // Create RecurrenceNnsumSpaces.  On for each unit tensor state parity
  reserve(Nnsum_partition.size());
  for (const auto& [Nnsum, partition] : Nnsum_partition)
  {
    auto subspace = RecurrenceNnsumSpace(
        Nnsum, partition, lgi_space_ket, lgi_space_bra, unit_tensor_constraints
      );
    if ((Nnsum == 0) && (subspace.dimension() == 0))
      return;  // TODO: is the Nnsum constraint necessary, or is this more general?
    if (subspace.dimension() == 0)
      continue;
    PushSubspace(std::move(subspace));
  }
}

RecurrenceSpace::RecurrenceSpace(
    const spatial::Space& space_ket,
    const spatial::Space& space_bra,
    const int& N1v,
    const HalfInt& Nsigma0
  )
    : BaseSpace{}
{
  reserve(2*space_ket.size()*space_bra.size());
  for (int i_ket = 0; i_ket < space_ket.size(); ++i_ket)
    for (int i_bra = 0; i_bra < space_bra.size(); ++i_bra)
    {
      const auto& lgi_space_ket = space_ket.GetSubspace(i_ket);
      const auto& lgi_space_bra = space_bra.GetSubspace(i_bra);
      for (uint8_t parity_bar : {0, 1})
      {
        auto subspace = RecurrenceLGISpace(
            lgi_space_ket, lgi_space_bra, {N1v, Nsigma0, parity_bar}
          );
        if (subspace.dimension() == 0) continue;
        PushSubspace(std::move(subspace));
     }
    }
}


}  // namespace spatial


}  // namespace spncci
