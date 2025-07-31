/****************************************************************
  u3_tensor.h

  RMEs for various U(3) tensor operators, including Sp(3,R) generators,
  coupled products of generators, Q, L^2, and Q dot Q. 

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/02/24 (cvc): Created. Functions moved from u3_racah_product.h.
  07/02/24 (cvc): Fixed issue with multiplicity of SU3GeneratorSp3R().
  07/03/24 (cvc): Q dot Q implemented for mass quadrupole.

****************************************************************/
#ifndef U3_TENSOR_H_
#define U3_TENSOR_H_

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

namespace sp3r_operator
{

// Wrapper for Identity
basis::OperatorBlock<double> IdentitySp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Wrapper for symplectic raising and lowering operators
// laddering between states in w', w subspaces of Sp(3,R) irrep defined by base_space
// Raising Operator, A^(2,0)
basis::OperatorBlock<double> RaisingOperatorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Lowering Operator, B^(0,2)
basis::OperatorBlock<double> LoweringOperatorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// SU(3) Generators, C^(1,1)
basis::OperatorBlock<double> SU3GeneratorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Harmonic Oscillator Hamiltonian, H_0^(0,0)
basis::OperatorBlock<double> HarmonicOscillatorSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Mass Quadrupole Operator Q_2 = sqrt(3)*[C^(1,1) + A^(2,0) + B^(0,2)]_2
basis::OperatorBlock<double> MassQuadrupoleSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// SU(3) Casimir
basis::OperatorBlock<double> SU3CasimirSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

///////////////////////////////////////////
//////////// COUPLED OPERATORS ////////////
///////////////////////////////////////////

// Coupled SU(3) Generators [C^(1,1) x C^(1,1)]
basis::OperatorBlock<double> CtimesCSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Coupled Raising Operators [A^(2,0) x A^(2,0)]
basis::OperatorBlock<double> AtimesASp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Coupled Lowering Operators [B^(0,2) x B^(0,2)]
basis::OperatorBlock<double> BtimesBSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Coupled Symplectic Raising/Lowering Operators [A^(2,0) x B^(0,2)]
basis::OperatorBlock<double> AtimesBSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Coupled Symplectic Loweding/Raising Operators [B^(0,2) x A^(2,0)]
basis::OperatorBlock<double> BtimesASp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

// Coupled Generator/Raising Operators [C^(1,1) x A^(2,0)]
basis::OperatorBlock<double> CtimesASp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

// Coupled Raising/Generator Operators [A^(2,0) x C^(1,1)]
basis::OperatorBlock<double> AtimesCSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

// Coupled Generator/Lowering Operators [C^(1,1) x B^(0,2)]
basis::OperatorBlock<double> CtimesBSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

// Coupled Lowering/Generator Operators [B^(0,2) x C^(1,1)]
basis::OperatorBlock<double> BtimesCSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Angular momentum operator L^2 = -sqrt(3)*[C_1 x C_1]_00
// L^2 = alpha*[C^(1,1) x C^(1,1)]^(0,0) + beta*[C^(1,1) x C^(1,1)]^(2,2)
basis::OperatorBlock<double> LdotLSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Elliott Quadrupole operator QdotQ = 3*sqrt(5)*[C_2 x C_2]_00
basis::OperatorBlock<double> QdotQElliottSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Sp(3,R) Mass Quadrupole operator QdotQ = 3*sqrt(5)*[(C_2 + A_2 + B_2) x (C_2 + A_2 + B_2)]_00
basis::OperatorBlock<double> QdotQSp3R(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

}

#endif