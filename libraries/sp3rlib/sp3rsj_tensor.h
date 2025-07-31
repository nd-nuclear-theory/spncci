/****************************************************************
  sp3rsj_tensor.h

  RMEs for Sp(3,R)xSU(2)>SU(2) operators.

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/
#ifndef SP3RSJ_TENSOR_H_
#define SP3RSJ_TENSOR_H_

#include "sp3rlib/sp3r_operator.h"

#include <functional>
#include <fstream>
#include <tuple>
#include <vector>

#include "am/halfint.h"

#include "sp3rlib/u3.h"
#include "mcutils/eigen.h"

#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3_racah_product.h"
#include "sp3rlib/u3_tensor.h"

namespace sp3r_operator
{

// Wrapper for J^2
basis::OperatorBlock<double> JdotJSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache,
  const HalfInt& S_bra,
  const HalfInt& S_ket,
  const HalfInt& J_bra,
  const HalfInt& J_ket
);

// Wrapper for S^2
basis::OperatorBlock<double> SdotSSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache,
  const HalfInt& S_bra,
  const HalfInt& S_ket,
  const HalfInt& J_bra,
  const HalfInt& J_ket
);

// Wrapper for LdotS = (1/2)*(J^2 - S^2 - L^2)
basis::OperatorBlock<double> LdotSSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache,
  const HalfInt& S_bra,
  const HalfInt& S_ket,
  const HalfInt& J_bra,
  const HalfInt& J_ket
);


}

#endif