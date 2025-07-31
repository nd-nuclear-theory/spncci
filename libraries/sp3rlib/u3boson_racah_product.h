/****************************************************************
  u3boson_racah_product.h

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  06/03/24 (cvc): Created. Refactored to .h/.cpp file format
****************************************************************/
#ifndef U3BOSON_RACAH_PRODUCT_H_
#define U3BOSON_RACAH_PRODUCT_H_

#include "sp3rlib/sp3r_operator.h"

#include <iostream>
#include <functional>

#include <fstream>
#include <iterator>
#include <map>
#include <tuple>
#include <vector>
#include <utility>

#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"

#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3.h"

// Ask about what to do with wrappers

namespace u3boson_operator
{

// Wrapper for Identity
basis::OperatorBlock<double> IdentityU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Wrapper for Creation and Annihiliation Operators
// FOR a and a^\dagger (ladder 2 U(3) quanta)!!!! 
// NOT b and b^\dagger (Weyl boson operators, 1 quanta)!!!!
// Need better notation
// Creation Operator <bra ||a^\dagger^(2,0)|| ket>_(rho)
basis::OperatorBlock<double> CreationU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Annhiliation Operator <bra ||a^(0,2)|| ket>_(rho)
basis::OperatorBlock<double> AnnihiliationU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Wrapper for SU(3) Casimir Operator
basis::OperatorBlock<double> SU3CasimirU3Boson(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Couple U(3) operators [T^wT x S^wS]^w0
// Matrix elements of the form <bra || [T^wT x S^wS]^rho0,w0 || ket>_(outer_mult)
// 
// Function inputs currently must have the form below:
// Tfunc(wT, rhoT, base_space, bra_space, ket_space, u_coef_cache) = <bra ||T^wT|| ket>_(rhoT)
// Sfunc(wS, rhoS, base_space, bra_space, ket_space, u_coef_cache) = <bra ||S^wS|| ket>_(rhoS)
//
// Note: currently U(3)-boson<Sp(3,R) operators within a single Sp(3,R) irrep
basis::OperatorBlock<double> CoupledU3BosonOperators(
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const u3boson::U3BosonSpace&,
    const u3boson::U3Subspace&,
    const u3boson::U3Subspace&,
    u3::UCoefCache&
    )> Tfunc,
  std::function<basis::OperatorBlock<double>(
    const u3::U3&,
    int,
    const u3boson::U3BosonSpace&,
    const u3boson::U3Subspace&,
    const u3boson::U3Subspace&,
    u3::UCoefCache&
    )> Sfunc,
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache  
);


// commutator [a,a^\dagger]^omega_0
// no restrictions on labels here, can use to test failure if incorrectly assigned
basis::OperatorBlock<double> CommutatorBosonQuanta(
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// commutator [a^\dagger,a^\dagger]^omega
// can be passed into coupling function CoupledU3BosonOperators()
basis::OperatorBlock<double> CommutatorCreationCreation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// commutator [a,a]^omega
// can be passed into coupling function CoupledU3BosonOperators()
basis::OperatorBlock<double> CommutatorAnnihilationAnnihilation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// commutator [a,a^\dagger]^omega
// can be passed into coupling function CoupledU3BosonOperators()
basis::OperatorBlock<double> CommutatorAnnihilationCreation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


/////////////////////////////////////////////
/////////////////////////////////////////////
/////////////////  TESTING  /////////////////
/////////////////////////////////////////////
/////////////////////////////////////////////

// number operator [a^\dagger x a]^omega
basis::OperatorBlock<double> CreationAnnihilation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// commutator [[a^\dagger,a]^0(0,0), a^\dagger]^omega
basis::OperatorBlock<double> CommutatorNumberCreation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// commutator [[a^\dagger,a]^0(0,0), a]^omega
basis::OperatorBlock<double> CommutatorNumberAnnihilation(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// Generators
basis::OperatorBlock<double> SU3GeneratorC(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);


// SU3 Casimir
basis::OperatorBlock<double> SU3CasimirCoupling(
  const u3::U3& omega,
  int rho,
  const u3boson::U3BosonSpace& base_space,
  const u3boson::U3Subspace& bra_subspace,
  const u3boson::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
);

}

#endif