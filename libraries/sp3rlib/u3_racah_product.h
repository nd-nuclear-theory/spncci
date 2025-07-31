/****************************************************************
  u3_racah_product.h

  RMEs for coupled product of U(3) tensor operators.

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  06/03/24 (cvc): Created. Split from u3_racah_product.cpp. 
  06/03/24 (cvc): Refactored to .h/.cpp file format.
  07/02/24 (cvc): Functions for specific tensors moved to u3_tensor.h.

****************************************************************/
#ifndef U3_RACAH_PRODUCT_H_
#define U3_RACAH_PRODUCT_H_

#include "sp3rlib/sp3r_operator.h"

#include <functional>
#include <fstream>
#include <tuple>
#include <vector>

#include "am/halfint.h"

#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "mcutils/eigen.h"

#include "sp3rlib/sp3r.h"

namespace sp3r_operator
{

// Couple Sp(3,R)>U(3) operators [T^wT x S^wS]^w0
//
// Function inputs currently must have the form below:
// Tfunc(wT, rhoT, base_space, bra_space, ket_space, u_coef_cache) = <bra ||T^wT|| ket>_(rhoT)
// Sfunc(wS, rhoS, base_space, bra_space, ket_space, u_coef_cache) = <bra ||S^wS|| ket>_(rhoS)
//
basis::OperatorBlock<double> CoupledU3Operators(
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const sp3r::Sp3RSpace&,
    const sp3r::U3Subspace&,
    const sp3r::U3Subspace&,
    u3::UCoefCache&
    )> Tfunc,
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const sp3r::Sp3RSpace&,
    const sp3r::U3Subspace&,
    const sp3r::U3Subspace&,
    u3::UCoefCache&
    )> Sfunc,
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

}

#endif
