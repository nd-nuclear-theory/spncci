/****************************************************************
  sp3rsj.cpp

  Colin V. Coane

  SPDX-License-Identifier: MIT

  10/17/23 (cvc): Created.
  01/25/24 (cvc): Added class for operator space and subspace structure.
  06/03/24 (cvc): Added class for Sp(3,R)xSU(2)>SU(2) observable.
  06/10/24 (cvc): Fixed multiplicity indexing for hypersectors.
  07/05/24 (cvc): Added bookkeeping for observables as linear combinations of operators.
  07/08/24 (cvc): Added ability to compute Sp(3,R)xSU(2)>SU(2) RMEs to OperatorInfoSp3RSJ.

****************************************************************/

#include "sp3rlib/sp3rsj.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>

#include "sp3rlib/u3boson.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"
#include "am/am.h"
#include "am/halfint.h"
#include "basis/hypersector.h"

namespace sp3r
{

Sp3RSJSpace::Sp3RSJSpace(
    const u3::U3& sigma,
    unsigned int Nn_max,
    const HalfInt& s_const,
    const HalfInt& j_const,
    const bool cache_Kmatrices,
    const bool branch_to_so3
  )
    : BaseSpace{sigma}, Nn_max_(Nn_max), s_const_(s_const), j_const_(j_const)//, sp3r_space_(sp3r::Sp3RSpace(sigma,Nn_max))
{
  // Check that sigma is an LGI of a unitary Sp(3,R) irrep
  assert(IsUnitary(sigma));
  const HalfInt::pair l_range = am::ProductAngularMomentumRange(s_const,j_const);
  // Check allowed coupling of S and J results in integer values of L
  assert(IsInteger(l_range.second));

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

    // Allowed angular momentum range
    const std::pair<unsigned int,unsigned int>& l_pair={l_range.first.TwiceValue()/2,l_range.second.TwiceValue()/2};
    const auto& [Lmin,Lmax] = l_pair;
    // Check branching
    bool non_empty_subspace = false;
    const auto& subspace_so3_labels = u3::BranchingSO3(omega.SU3());
    // Check if subspace is empty
    for (const auto& [L, kappa_max] : subspace_so3_labels)
    {
      if (L >= Lmin && L<= Lmax)
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

  // Construct associated untruncated Sp3RSpace
  //sp3r::Sp3RSpace sp3r_space_(sigma,Nn_max);
};

std::string Sp3RSJSpace::DebugStr() const
{
  std::ostringstream ss;

  // print space labels
  ss << "space sigma " << sigma().Str() << " Nn_max " << Nn_max_ << " S " << s_const_ << " J " << j_const_ << std::endl;

  // iterate over subspaces
  for (std::size_t i_subspace = 0; i_subspace < size(); ++i_subspace)
  {
    const sp3r::U3Subspace& subspace = GetSubspace(i_subspace);
    ss << subspace.DebugStr();
  }
  return ss.str();
}



std::string OperatorU3Subspace::DebugStr() const
{
  std::ostringstream ss;
  ss << fmt::format(
      "subspace: {}\nsize: {}  dimension: {}",
      omega(),
      size(),
      dimension()
    ) << std::endl;


  for (std::size_t i_state = 0; i_state < size(); ++i_state)
  {
    ss << fmt::format(
        "L {} : kappa {}", GetState(i_state).L(), GetState(i_state).kappa()
      ) << std::endl;
  }

  return ss.str();
}


std::string OperatorSpaceSp3RSJ::DebugStr() const
{
  std::ostringstream ss;

  // print space labels
  ss << fmt::format(
      "space: J0 {}  S0 {}",
      J0(),
      S0()
    ) << std::endl;

  // iterate over subspaces
  for (std::size_t i_subspace = 0; i_subspace < size(); ++i_subspace)
  {
    const sp3r::OperatorU3Subspace& subspace = GetSubspace(i_subspace);
    ss << subspace.DebugStr();
  }
  return ss.str();
}

////////////////////////////////////////////////////////////////
// space and subspace indexing
////////////////////////////////////////////////////////////////

// Sectors
Sp3RSJSectors::Sp3RSJSectors(
    const Sp3RSJSpace& space, const u3::U3& omega0, const bool& su3_generator
  )
    : BaseSectors{space, space}, omega0_(omega0)
{
  // omega0_ = omega0;
  for (std::size_t bra_index = 0; bra_index < space.size(); bra_index++)
    for (std::size_t ket_index = 0; ket_index < space.size(); ket_index++)
    {
      const SubspaceType& bra_subspace = space.GetSubspace(bra_index);
      const SubspaceType& ket_subspace = space.GetSubspace(ket_index);
      const u3::U3& omega_bra = bra_subspace.U3();
      const u3::U3& omega_ket = ket_subspace.U3();

      if (su3_generator && (omega_bra != omega_ket))
        continue;

      // verify selection rules
      bool allowed = true;

      unsigned int multiplicity = u3::OuterMultiplicity(omega_ket, omega0, omega_bra);
      allowed &= (multiplicity > 0);

      if (allowed)
        PushSector(bra_index, ket_index, multiplicity);
    }
};

// Hypersectors
Sp3RSJHypersectors::Sp3RSJHypersectors(
    const Sp3RSJSpace& space,
    const OperatorSpaceSp3RSJ& operator_space,
    const bool& su3_generator
  )
{
  //
  for (std::size_t bra_index = 0; bra_index < space.size(); bra_index++)
    for (std::size_t ket_index = 0; ket_index < space.size(); ket_index++)
      for (int op_index=0; op_index<operator_space.size(); op_index++)
      {

        const SubspaceType& bra_subspace = space.GetSubspace(bra_index);
        const SubspaceType& ket_subspace = space.GetSubspace(ket_index);
        const OperatorSubspaceType& op_subspace = operator_space.GetSubspace(op_index);

        const u3::U3& omega_bra = bra_subspace.U3();
        const u3::U3& omega_ket = ket_subspace.U3();
        const u3::U3& omega_op = op_subspace.U3();

        // verify selection rules
        bool allowed = true;

        unsigned int multiplicity = u3::OuterMultiplicity(omega_ket, omega_op, omega_bra);
        allowed &= (multiplicity > 0);

        if (allowed)
          for (int mult_idx = 1; mult_idx <= multiplicity; mult_idx++)
            PushHypersector(HypersectorType(
                              bra_index,ket_index,op_index,
                              bra_subspace,ket_subspace,op_subspace,mult_idx));
        
        // SHOULD THERE BE MULTIPLICITY HERE OR NOT????????

      }

}

}