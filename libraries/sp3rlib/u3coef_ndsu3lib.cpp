/****************************************************************
  u3coef.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include <algorithm>
#include <cassert>

#include "c++wrappers.h"  // from ndsu3lib
#include "fmt/format.h"
#include "sp3rlib/u3coef.h"


namespace u3
{

////////////////////////////////////////////////////////////////
// direct access to su3lib FORTRAN subroutines
////////////////////////////////////////////////////////////////

namespace ndsu3lib
{


ndsu3lib::su3irrep ndsu3(const u3::SU3& x)
{
  return ndsu3lib::su3irrep(x.lambda(), x.mu());
}


void U3CoefInit(int max_lambda_plus_mu)
{
  ndsu3lib_init(true, max_lambda_plus_mu);
}


std::size_t WCoefBlock::CoefIndex(
    const unsigned int kappa1,
    const unsigned int kappa2,
    const unsigned int kappa3,
    const unsigned int rho
  ) const
{
  // SU3lib convention is mulitplicity index is 0<=index<multiplicity
  std::size_t index = rho - 1;
  index = index * kappa3_max_ + (kappa3 - 1);
  index = index * kappa2_max_ + (kappa2 - 1);
  index = index * kappa1_max_ + (kappa1 - 1);
  return index;
}

WCoefBlock::WCoefBlock(const u3::WCoefLabels& labels)
{
  // calculate multiplicities
  const auto& [x1, L1, x2, L2, x3, L3] = labels.Key();
  const auto& [lambda1, mu1] = x1.Key();
  const auto& [lambda2, mu2] = x2.Key();
  const auto& [lambda3, mu3] = x3.Key();
  std::tie(kappa1_max_, kappa2_max_, kappa3_max_, rho_max_) =
      WMultiplicity(x1, L1, x2, L2, x3, L3);

  int dimension =
      static_cast<int>(kappa1_max_ * kappa2_max_ * kappa3_max_ * rho_max_);
  coefs_.resize(dimension);
  auto& wigner = coefs_.data();
  ndsu3lib::calculate_wigner_su3so3(
      ndsu3(x1),
      static_cast<int>(L1),
      ndsu3(x2),
      static_cast<int>(L2),
      ndsu3(x3),
      static_cast<int>(L3),
      static_cast<int>(kappa1_max_),
      static_cast<int>(kappa2_max_),
      static_cast<int>(kappa3_max_),
      static_cast<int>(rho_max_),
      dimension,
      wigner
    );
}

//////////////////////////////////////////////////////////////////
// Same for U and Z
//////////////////////////////////////////////////////////////////

std::size_t RecouplingCoefBlock::CoefIndex(
    unsigned int r12, unsigned int r12_3, unsigned int r23, unsigned int r1_23
  ) const
{
  std::size_t index = (r1_23 - 1);
  index = index * r23_max_ + (r23 - 1);
  index = index * r12_3_max_ + (r12_3 - 1);
  index = index * r12_max_ + (r12 - 1);
  return index;
}


RecouplingCoefBlock::RecouplingCoefBlock(
    const u3::UCoefLabels& labels, const u3::RecouplingMode mode
  )
{
  // calculate multiplicities
  const auto& [x1, x2, x, x3, x12, x23] = labels.Key();
  std::tie(r12_max_, r12_3_max_, r23_max_, r1_23_max_) =
      UMultiplicity(x1, x2, x, x3, x12, x23);

  mode_ = mode;

  // pre-size vector
  int dimension = static_cast<int>(r12_max_ * r12_3_max_ * r23_max_ * r1_23_max_);
  assert(dimension > 0);

  coefs_.resize(dimension);
  auto& racah = coefs_.data();
  int info;  // for MKL subroutine.  Not used here

  // populate vector
  if (mode == RecouplingMode::kZ)
  {
    ndsu3lib::calculate_z_coeff(
        ndsu3(x1),
        ndsu3(x2),
        ndsu3(x),
        ndsu3(x3),
        ndsu3(x12),
        ndsu3(x23),
        static_cast<int>(r12_max_),
        static_cast<int>(r12_3_max_),
        static_cast<int>(r23_max_),
        static_cast<int>(r1_23_max_),
        dimension,
        racah,
        info
      );
  }
  // mode is kU
  else
  {
    ndsu3lib::calculate_u_coeff(
        ndsu3(x1),
        ndsu3(x2),
        ndsu3(x),
        ndsu3(x3),
        ndsu3(x12),
        ndsu3(x23),
        static_cast<int>(r12_max_),
        static_cast<int>(r12_3_max_),
        static_cast<int>(r23_max_),
        static_cast<int>(r1_23_max_),
        dimension,
        racah,
        info
      );
  }
}


PhiCoefBlock::PhiCoefBlock(const u3::PhiCoefLabels& labels)
{
  // calculate multiplicities
  u3::SU3 x1, x2, x3;
  std::tie(x1, x2, x3) = labels.Key();
  rho_max_ = u3::OuterMultiplicity(x1, x2, x3);

  auto& racah = coefs_.data();
  int info;  // for MKL subroutine.  Not used here.

  ndsu3lib::calculate_z_coeff(
      ndsu3(x1),
      ndsu3({0, 0}),
      ndsu3(x3),
      ndsu3(x2),
      ndsu3(x1),
      ndsu3(x2),
      1,
      rho_max_,
      1,
      rho_max_,
      rho_max_ * rho_max_,
      racah,
      info
    );
}

}  // namespace ndsu3lib

}  // namespace u3
