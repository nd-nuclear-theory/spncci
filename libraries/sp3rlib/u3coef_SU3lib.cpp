/****************************************************************
  u3coef.cpp

  SU(3) coupling coefficient wrappers for Dytrych and Langr

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include <algorithm>
#include <cassert>

#include "fmt/format.h"
#include "sp3rlib/u3coef.h"
#include "su3.h"

namespace u3
{

void U3CoefInit(int max_lambda_plus_mu)
{
  // max_lambda_plus_mu*=2;
  int jjmax = 2 * max_lambda_plus_mu;
  su3::init_thread(2 * max_lambda_plus_mu, jjmax);
}

std::size_t WCoefBlock::CoefIndex(
    const unsigned int kappa1,
    const unsigned int kappa2,
    const unsigned int kappa3,
    const unsigned int rho
  ) const
{
  // SU3lib convention is mulitplicity index is 0<=index<multiplicity
  std::size_t index = (kappa3 - 1);
  index = index * kappa2_max_ + (kappa2 - 1);
  index = index * kappa1_max_ + (kappa1 - 1);
  index = index * rho_max_ + (rho - 1);
  return index;
}


WCoefBlock::WCoefBlock(const u3::WCoefLabels& labels)
{
  // calculate multiplicities
  const auto& [x1, L1, x2, L2, x3, L3] = labels.Key();
  std::tie(kappa1_max_, kappa2_max_, kappa3_max_, rho_max_) =
      WMultiplicity(x1, L1, x2, L2, x3, L3);

  if (am::AllowedTriangle(L1, L2, L3))
    su3::wu3r3w(
        x1.lambda(),
        x1.mu(),
        x2.lambda(),
        x2.mu(),
        x3.lambda(),
        x3.mu(),
        static_cast<int>(L1),
        static_cast<int>(L2),
        static_cast<int>(L3),
        static_cast<int>(kappa1_max_),
        static_cast<int>(kappa2_max_),
        static_cast<int>(kappa3_max_),
        static_cast<int>(rho_max_),
        coefs_
      );
  else
    coefs_ = std::vector<double>(
        kappa1_max_ * kappa2_max_ * kappa3_max_ * rho_max_, 0.0
      );
}

// ////////////////////////////////////////////////////////////////

// Same for U and Z
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
  unsigned int r_max = r12_max_ * r12_3_max_ * r23_max_ * r1_23_max_;
  assert(r_max > 0);

  // populate vector
  if (mode == RecouplingMode::kZ)
  {
    su3::wzu3(
        x1.lambda(),
        x1.mu(),
        x2.lambda(),
        x2.mu(),
        x.lambda(),
        x.mu(),
        x3.lambda(),
        x3.mu(),
        x12.lambda(),
        x12.mu(),
        x23.lambda(),
        x23.mu(),
        static_cast<int>(r12_max_),
        static_cast<int>(r12_3_max_),
        static_cast<int>(r23_max_),
        static_cast<int>(r1_23_max_),
        coefs_
      );
  }
  // mode is kU
  else
  {
    su3::wru3(
        x1.lambda(),
        x1.mu(),
        x2.lambda(),
        x2.mu(),
        x.lambda(),
        x.mu(),
        x3.lambda(),
        x3.mu(),
        x12.lambda(),
        x12.mu(),
        x23.lambda(),
        x23.mu(),
        static_cast<int>(r12_max_),
        static_cast<int>(r12_3_max_),
        static_cast<int>(r23_max_),
        static_cast<int>(r1_23_max_),
        coefs_
      );
  }
}


PhiCoefBlock::PhiCoefBlock(const u3::PhiCoefLabels& labels)
{
  // calculate multiplicities
  u3::SU3 x1, x2, x3;
  std::tie(x1, x2, x3) = labels.Key();
  rho_max_ = u3::OuterMultiplicity(x1, x2, x3);
  // zero initialize array
  su3::wzu3(
  // 05/11/23 (jh): Change wru3 to wzu3
      x1.lambda(),
      x1.mu(),
      0,
      0,
      x3.lambda(),
      x3.mu(),
      x2.lambda(),
      x2.mu(),
      x1.lambda(),
      x1.mu(),
      x2.lambda(),
      x2.mu(),
      1,
      static_cast<int>(rho_max_),
      1,
      static_cast<int>(rho_max_),
      coefs_
    );
}


}  // namespace u3
