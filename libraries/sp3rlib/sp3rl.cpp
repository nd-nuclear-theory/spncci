/****************************************************************
  sp3rl.cpp

  SPDX-License-Identifier: MIT

  10/14/23 (cvc): Created.

****************************************************************/

#include "sp3rlib/sp3rl.h"
//#include "sp3rlib/sp3r.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>

#include "sp3rlib/u3boson.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"
namespace sp3r
{

Sp3RLSpace::Sp3RLSpace(
    const u3::U3& sigma,
    unsigned int Nn_max,
    unsigned int l_const, // currently passing definite l, can try more complex selection rule later
    const bool cache_Kmatrices,
    const bool branch_to_so3
  )
    : BaseSpace{sigma}, Nn_max_(Nn_max), l_const_(l_const)
{
  // Check that sigma is an LGI of a unitary Sp(3,R) irrep
  assert(IsUnitary(sigma));

  // Construct the assocaited u3boson space
  u3boson::U3BosonSpace u3boson_space(sigma, Nn_max);

  // If constructing the full space, then compute K matrices and
  // use K matrices to get upsilon_max.  If subspace labels only,
  // K matrices only need to be computed if branching must be restricted.
  vcs::KmatrixMap K_matrices;
  bool get_upsilon_from_K = false;
  if (cache_Kmatrices == true || sp3r::ModifySp3RBranching(sigma))
  {
    double zero_threshold = 1e-12;
    K_matrices = vcs::GenerateKmatrices(sigma, u3boson_space, zero_threshold);
    get_upsilon_from_K = true;
  }

  for (std::size_t i = 0; i < u3boson_space.size(); ++i)
  {
    const auto& u3boson_subspace = u3boson_space.GetSubspace(i);
    const u3::U3& omega = u3boson_subspace.omega();
    unsigned int upsilon_max =
        get_upsilon_from_K
            ? static_cast<unsigned int>(K_matrices[omega][0].rows())
            : static_cast<unsigned int>(u3boson_subspace.dimension());

    if (upsilon_max == 0)
      continue;

    const std::pair<unsigned int,unsigned int>& l_pair={l_const,l_const};
    // Check branching
    bool non_empty_subspace = false;
    const auto& subspace_so3_labels = u3::BranchingSO3(omega.SU3());
    // Check if subspace is empty
    for (const auto& [L, kappa_max] : subspace_so3_labels)
    {
      if (L == l_const)
        non_empty_subspace = true;
    };
    // Push subspace if non-empty
    if (non_empty_subspace)
    {
    // Labels only
    if (cache_Kmatrices)
      PushSubspace(U3Subspace(
        omega,
        upsilon_max,
        u3boson_space.GetSubspacePtr(i),
        std::move(K_matrices[omega][0]),
        std::move(K_matrices[omega][1]),
        branch_to_so3,
        l_pair
      ));


    // Full space
    else
      PushSubspace(U3Subspace(
        omega,
        upsilon_max,
        branch_to_so3,
        l_pair
      ));
    }
  }
}

std::string Sp3RLSpace::DebugStr() const
{
  std::ostringstream ss;

  // print space labels
  ss << "space sigma " << sigma().Str() << " Nn_max " << Nn_max_ << " L " << l_const_ << std::endl;

  // iterate over subspaces
  for (std::size_t i_subspace = 0; i_subspace < size(); ++i_subspace)
  {
    const sp3r::U3Subspace& subspace = GetSubspace(i_subspace);
    ss << subspace.DebugStr();
  }
  return ss.str();
}

////////////////////////////////////////////////////////////////
// space and subspace indexing
////////////////////////////////////////////////////////////////

Sp3RLSectors::Sp3RLSectors(
    const Sp3RLSpace& space, const u3::U3& omega0, const bool& su3_generator
  )
    : BaseSectors{space, space}, omega0_(omega0)
{
  // omega0_ = omega0;
  for (std::size_t bra_index = 0; bra_index < space.size(); ++bra_index)
    for (std::size_t ket_index = 0; ket_index < space.size(); ++ket_index)
    {
      const u3::U3& omega_bra = space.GetSubspace(bra_index).U3();
      const u3::U3& omega_ket = space.GetSubspace(ket_index).U3();

      if (su3_generator && (omega_bra != omega_ket))
        continue;

      unsigned int multiplicity =
          u3::OuterMultiplicity(omega_ket, omega0, omega_bra);
      if (multiplicity > 0)
        PushSector(bra_index, ket_index, multiplicity);
    }
}

}