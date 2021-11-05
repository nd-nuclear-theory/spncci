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
#include <cppitertools/enumerate.hpp>
#include <fstream>
#include <iostream>
#include <numeric>

#include "basis/basis.h"
// #include "mcutils/parsing.h"
#include "fmt/format.h"
// #include "am/halfint_fmt.h"

// #include "sp3rlib/vcs.h"

namespace spncci::spatial
{
Sp3RSpace::Sp3RSpace(const u3::U3& sigma, const int Nn_max)
    : BaseSpace{sigma}
{
  // find all raising polynomials
  std::vector<u3::U3> n_vector = sp3r::RaisingPolynomialLabels(Nn_max);

  // temporary container
  std::map<u3::U3, MultiplicityTagged<u3::U3>::vector> states;
  std::unordered_map<u3::U3, std::array<basis::OperatorBlock<double>, 2>> K_matrices;

  // For each raising polynomial n obtain all allowed couplings
  // omega (sigma x n -> omega) with outer multiplicities rho_max.
  for (const auto& n : n_vector)
  {
    MultiplicityTagged<u3::U3>::vector omega_rhomax_vector =
        u3::KroneckerProduct(sigma, n);
    for (const auto& [omega, rho_max] : omega_rhomax_vector)
      states[omega].emplace_back(n, rho_max);
  }

  // TODO(aem): populate K and Kinv
  for (const auto& [omega, n_rho_vector] : states)
  {
    std::size_t dim = std::accumulate(
        n_rho_vector.begin(), n_rho_vector.end(), 0,
        [](const int& s, const MultiplicityTagged<u3::U3>& n_rho)
        {
          const auto& [n,rho] = n_rho; return s+rho;
        }
        );
    K_matrices[omega] = {
        basis::OperatorBlock<double>::Identity(dim, dim),
        basis::OperatorBlock<double>::Identity(dim, dim)
      };
  }

  // Push subpaces
  for (const auto& [omega, n_rho_vector] : states)
    PushSubspace(U3Subspace(
        omega,
        n_rho_vector,
        std::move(K_matrices[omega][0]),
        std::move(K_matrices[omega][1])
      ));
}

Space::Space(
    const std::vector<u3::U3>& sigma_vector, const HalfInt& Nsigma0, const int Nmax
  )
    : BaseSpace{}
{
  for (const auto& sigma : sigma_vector)
  {
    // Nsigma0_=Nsigma0;
    int Nn_max = Nmax - int(sigma.N() - Nsigma0);
    if (Nn_max >= 0)
      PushSubspace(Sp3RSpace(sigma, Nn_max));
  }
}

////////////////////////////////////////////////////////////////////////////////
// Recurrence indexing (spatial)
////////////////////////////////////////////////////////////////////////////////


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
  std::unordered_map<u3::SU3, std::vector<std::tuple<int, int>>, boost::hash<u3::SU3>>
      x0_Nbar_pairs;

  //////////////////////////////////////////////////////////////////////////////
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

  //////////////////////////////////////////////////////////////////////////////
  // Construct operator subspaces
  //////////////////////////////////////////////////////////////////////////////
  reserve(x0_Nbar_pairs.size());
  for (const auto& [x0, Nbar_pairs] : x0_Nbar_pairs)
  {
    int rho0_max =
        u3::OuterMultiplicity(omega_ket().SU3(), x0, omega_bra().SU3());
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
    const Sp3RSpace& sp3r_space_ket,
    const Sp3RSpace& sp3r_space_bra,
    const UnitTensorConstraintParameters& unit_tensor_constraints
  )
    : BaseDegenerateSpace{Nnsum}
{
  parity_bar_ = unit_tensor_constraints.parity_bar;

  reserve(u3subspace_index_pairs.size());
  upsilon_max_pairs_.reserve(u3subspace_index_pairs.size());

  for (const auto& [i_ket, i_bra] : u3subspace_index_pairs)
  {
    const auto& u3subspace_ket = sp3r_space_ket.GetSubspace(i_ket);
    const auto& u3subspace_bra = sp3r_space_bra.GetSubspace(i_bra);
    const u3::U3& omega_ket = u3subspace_ket.omega();
    const u3::U3& omega_bra = u3subspace_bra.omega();

    const int& upsilon_max_ket = u3subspace_ket.size();
    const int& upsilon_max_bra = u3subspace_bra.size();

    auto subspace =
        RecurrenceU3Space({omega_ket, omega_bra}, unit_tensor_constraints);
    if (subspace.dimension() == 0)
      continue;

    upsilon_max_pairs_.push_back({upsilon_max_ket, upsilon_max_bra});
    PushSubspace(std::move(subspace), upsilon_max_bra * upsilon_max_ket);
  }
}


RecurrenceSp3RSpace::RecurrenceSp3RSpace(
    const Sp3RSpace& sp3r_space_ket,
    const Sp3RSpace& sp3r_space_bra,
    const UnitTensorConstraintParameters& unit_tensor_constraints
  )
    : BaseSpace{
        {sp3r_space_ket.sigma(),
         sp3r_space_bra.sigma(),
         unit_tensor_constraints.parity_bar}
      }
{
  const u3::U3& sigma_ket = sp3r_space_ket.sigma();
  const u3::U3& sigma_bra = sp3r_space_bra.sigma();
  // std::cout<<sigma_bra.Str()<<"  "<<sigma_ket.Str()<<std::endl;

  // Partition pairs of omega',omega by Nnsum
  std::map<int, std::vector<std::tuple<int, int>>> Nnsum_partition;
  for (int i_ket = 0; i_ket < sp3r_space_ket.size(); ++i_ket)
    for (int i_bra = 0; i_bra < sp3r_space_bra.size(); ++i_bra)
    {
      const u3::U3& omega_ket = sp3r_space_ket.GetSubspace(i_ket).omega();
      const u3::U3& omega_bra = sp3r_space_bra.GetSubspace(i_bra).omega();

      int Nnsum =
          int(omega_ket.N() - sigma_ket.N() + omega_bra.N() - sigma_bra.N());
      Nnsum_partition[Nnsum].push_back({i_ket, i_bra});
    }

  // Create RecurrenceNnsumSpaces.  On for each unit tensor state parity
  reserve(Nnsum_partition.size());
  for (const auto& [Nnsum, partition] : Nnsum_partition)
  {
    auto subspace = RecurrenceNnsumSpace(
        Nnsum, partition, sp3r_space_ket, sp3r_space_bra, unit_tensor_constraints
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
  reserve(2 * space_ket.size() * space_bra.size());
  for (int i_ket = 0; i_ket < space_ket.size(); ++i_ket)
    for (int i_bra = 0; i_bra < space_bra.size(); ++i_bra)
    {
      const auto& sp3r_space_ket = space_ket.GetSubspace(i_ket);
      const auto& sp3r_space_bra = space_bra.GetSubspace(i_bra);
      for (uint8_t parity_bar : {0, 1})
      {
        auto subspace = RecurrenceSp3RSpace(
            sp3r_space_ket, sp3r_space_bra, {N1v, Nsigma0, parity_bar}
          );
        if (subspace.dimension() == 0)
          continue;
        PushSubspace(std::move(subspace));
      }
    }
}

RecurrenceU3Sectors::RecurrenceU3Sectors(
    const RecurrenceSp3RSpace& sp3r_space, int target_Nnsum, int source_Nnsum
  )
    : BaseSectors{
        sp3r_space.LookUpSubspace({target_Nnsum}),
        sp3r_space.LookUpSubspace({source_Nnsum})
      }
{
  int delta_Nnsum = target_Nnsum - source_Nnsum;
  assert((delta_Nnsum == 2) || (delta_Nnsum == 4));
  const auto& target_Nnsum_space = bra_space();
  const auto& source_Nnsum_space = ket_space();

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
                  u3::U3{2, {2, 0}},
                  source_u3_space.omega_bra(),
                  target_u3_space.omega_bra()
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
                  u3::U3{2, {2, 0}},
                  source_u3_space.omega_ket(),
                  target_u3_space.omega_ket()
                )
              == 0
            )
            continue;
        }
        EmplaceSector(
            target_u3_index, source_u3_index, target_u3_space, source_u3_space
          );
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
                u3::U3{2, {2, 0}},
                source_u3_space.omega_bra(),
                target_u3_space.omega_bra()
              )
            == 0
          )
          continue;
        if (
            u3::OuterMultiplicity(
                u3::U3{2, {2, 0}},
                source_u3_space.omega_ket(),
                target_u3_space.omega_ket()
              )
            == 0
          )
          continue;
        EmplaceSector(
            target_u3_index, source_u3_index, target_u3_space, source_u3_space
          );
      }
    }
  }
}

}  // namespace spncci::spatial
