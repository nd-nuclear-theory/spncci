/****************************************************************
  sp3r.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include "sp3rlib/sp3r.h"

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

//////////////////////////////////////////////////////////////////
//  Sp(3,R) raising polynomial
//////////////////////////////////////////////////////////////////
// TauConjugate computes the first two rows of the partition conjugate
// of tau = [f1-f3,f2-f3] = [lambda+mu,mu], defined in
// jpa-18-1985-939-Rowe. Note, in ref., tau is referred to as lambda.
std::array<int, 2> TauConjugate(const u3::U3& sigma)
{
  std::array<int, 2> tau_tilde = {0};
  const auto& lambda = sigma.SU3().lambda();
  const auto& mu = sigma.SU3().mu();

  if ((lambda + mu) >= 1)
    tau_tilde[0]++;

  if ((lambda + mu) >= 2)
    tau_tilde[1]++;

  if (mu >= 1)
    tau_tilde[0]++;

  if (mu >= 2)
    tau_tilde[1]++;

  return tau_tilde;
}

bool IsUnitary(const u3::U3& sigma)
// Based on criterion for unitary irrep in jpa-18-1985-939-Rowe
// equations (2.7) and (2.8).  Note, first criterion of
// conjugate{tau}_1<=2 always met for U(3) irrep, so we only check the
// second criterion conjugate{tau}_1+conjugate{tau}_2 <=2*f3.
{
  HalfInt f3 = sigma.f3();
  auto tau_conjugate = TauConjugate(sigma);
  return (tau_conjugate[0] + tau_conjugate[1]) <= 2 * f3;
}


bool ModifySp3RBranching(const u3::U3& sigma)
// From section 5 of jpa-18-1985-939-Rowe,
// there are U(3) irreps in the U(3) boson basis
// obtained via laddering which have no counter-part
// in the Sp(3,R) basis and thus must be removed.
// Such cases occur when f3<3.
{
  return sigma.f3() < 3;
}


Sp3RSpace::Sp3RSpace(
    const u3::U3& sigma,
    unsigned int Nn_max,
    const bool cache_Kmatrices,
    const bool branch_to_so3
  )
    : BaseSpace{sigma}, Nn_max_(Nn_max)
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

    // Labels only
    if (cache_Kmatrices)
      PushSubspace(U3Subspace(
        omega,
        upsilon_max,
        u3boson_space.GetSubspacePtr(i),
        std::move(K_matrices[omega][0]),
        std::move(K_matrices[omega][1]),
        branch_to_so3
      ));


    // Full space
    else
      PushSubspace(U3Subspace(omega, upsilon_max));
  }
}

////////////////////////////////////////////////////////////////////////
std::string U3Subspace::DebugStr() const
{
  std::ostringstream ss;
  ss << fmt::format(
      "subspace: {}\nupsilon_max: {}  size: {}  dimension: {}",
      omega(),
      upsilon_max(),
      size(),
      dimension()
    ) << std::endl;


  for (std::size_t i_state = 0; i_state < size(); ++i_state)
  {
    ss << fmt::format(
        "{} : {}", GetState(i_state).L(), GetState(i_state).kappa_max()
      ) << std::endl;
  }

  return ss.str();
}


std::string Sp3RSpace::DebugStr() const
{
  std::ostringstream ss;

  // print space labels
  ss << "space sigma " << sigma().Str() << " Nn_max " << Nn_max_ << std::endl;

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

Sp3RSectors::Sp3RSectors(
    const Sp3RSpace& space, const u3::U3& omega0, const bool& su3_generator
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

}  // namespace sp3r
