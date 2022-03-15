  /****************************************************************
  sp3r_operators.h                       

  Sp(3,R) generator and operator reduced matrix elements.

  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT
 
  3/9/16 (aem): Created based on prototype vcs.py and sp3r.py.

****************************************************************/
#ifndef SP3R_OPERATORS_H_
#define SP3R_OPERATORS_H_

#include "sp3rlib/u3.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/u3coef.h"

namespace sp3r
{

 
  basis::OperatorBlock<double> Sp3rRaisingOperator(
    const u3::U3& sigma,
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    );
  // Reduced matrix elements of symplectic raising operator between states in
  // omegap and omega subspace of Sp(3,R) irrep defiend by sp3r_space 


    basis::OperatorBlock<double> Sp3rLoweringOperator(
    const u3::U3& sigma,
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    );

    basis::OperatorBlock<double> SU3Generator(
      const u3::U3& sigma,
      const sp3r::U3Subspace& bra_subspace,
      const sp3r::U3Subspace& ket_subspace
    );
}

#endif
