/****************************************************************
  u3_racah_product_test.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  08/21/23 (cvc): Created. Split from u3_boson_racah_product.cpp.
  06/03/24 (cvc): Split from main function for obtaining RMEs.
  06/06/24 (cvc): Commutator tests ported over.

****************************************************************/

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

#include "sp3rlib/sp3r.h"
#include "sp3rlib/vcs.h"

#include "sp3rlib/u3_racah_product.h"
#include "sp3rlib/u3_tensor.h"



namespace sp3r_operator
{

// commutator [B,A]^omega
// no restrictions on labels here, can use to test failure if incorrectly assigned
basis::OperatorBlock<double> CommutatorLoweringRaising(
  const u3::U3& omega_T,
  const u3::U3& omega_S,
  const u3::U3& omega0,
  int rho0,
  int outer_multiplicity,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&LoweringOperatorSp3R,&RaisingOperatorSp3R,
                          omega_T,omega_S,omega0,rho0,outer_multiplicity,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&RaisingOperatorSp3R,&LoweringOperatorSp3R,
                          omega_S,omega_T,omega0,rho0,outer_multiplicity,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}

// commutator [B,A]^omega
// can be passed into coupling function CoupledU3Operators()
basis::OperatorBlock<double> CommutatorBA(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  int rho0 = 1;
  //
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&LoweringOperatorSp3R,&RaisingOperatorSp3R,
                          omega_B,omega_A,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&RaisingOperatorSp3R,&LoweringOperatorSp3R,
                          omega_A,omega_B,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}

// nested commutator [[B^(0,2),A^(2,0)]^(0,0),A^(2,0)]^omega = [sqrt(8/3)*H_0^(0,0),A^(2,0)]^omega
// can be passed into coupling function CoupledU3Operators()
basis::OperatorBlock<double> CommutatorOscillatorA(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_H(0,{0u,0u});
  int rho0 = 1;
  //
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&CommutatorBA,&RaisingOperatorSp3R,
                          omega_H,omega_A,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&RaisingOperatorSp3R,&CommutatorBA,
                          omega_A,omega_H,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}

// nested commutator [[B^(0,2),A^(2,0)]^(0,0),B^(0,2)]^omega = [sqrt(8/3)*H_0^(0,0),B^(0,2)]^omega
// can be passed into coupling function CoupledU3Operators()
basis::OperatorBlock<double> CommutatorOscillatorB(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_H(0,{0u,0u});
  int rho0 = 1;
  //
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&CommutatorBA,&LoweringOperatorSp3R,
                          omega_H,omega_B,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&LoweringOperatorSp3R,&CommutatorBA,
                          omega_B,omega_H,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}

// commutator [H_0,A]^omega
// can be passed into coupling function CoupledU3Operators()
basis::OperatorBlock<double> CommutatorHA(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_H(0,{0u,0u});
  int rho0 = 1;
  //
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&HarmonicOscillatorSp3R,&RaisingOperatorSp3R,
                          omega_H,omega_A,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&RaisingOperatorSp3R,&HarmonicOscillatorSp3R,
                          omega_A,omega_H,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}

// commutator [H_0,B]^omega
// can be passed into coupling function CoupledU3Operators()
basis::OperatorBlock<double> CommutatorHB(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_H(0,{0u,0u});
  int rho0 = 1;
  //
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&HarmonicOscillatorSp3R,&LoweringOperatorSp3R,
                          omega_H,omega_B,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&LoweringOperatorSp3R,&HarmonicOscillatorSp3R,
                          omega_B,omega_H,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


// commutator [B,C]^omega
// can be passed into coupling function CoupledU3Operators()
basis::OperatorBlock<double> CommutatorBC(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_C(0,{1u,1u});
  int rho0 = 1;
  //
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&LoweringOperatorSp3R,&SU3GeneratorSp3R,
                          omega_B,omega_C,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&SU3GeneratorSp3R,&LoweringOperatorSp3R,
                          omega_C,omega_B,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}

// commutator [A,C]^omega
// can be passed into coupling function CoupledU3Operators()
basis::OperatorBlock<double> CommutatorAC(
  const u3::U3& omega,
  int rho,
  const sp3r::Sp3RSpace& base_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  u3::UCoefCache& u_coef_cache
)
{
  const u3::U3 sigma = base_space.sigma();
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_C(0,{1u,1u});
  int rho0 = 1;
  //
  // commutator matrix
  basis::OperatorBlock<double> commutator_matrix = basis::OperatorBlock<double>::Zero(
      bra_subspace.dimension(),
      ket_subspace.dimension()
      );

  commutator_matrix += CoupledU3Operators(&RaisingOperatorSp3R,&SU3GeneratorSp3R,
                          omega_A,omega_C,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  commutator_matrix -= CoupledU3Operators(&SU3GeneratorSp3R,&RaisingOperatorSp3R,
                          omega_C,omega_A,omega,rho0,rho,
                          base_space,bra_subspace,ket_subspace,u_coef_cache);
  
  return commutator_matrix;
}


}




// Test couplings
int main(int argc, char **argv)
{
  u3::U3CoefInit(39); // what size should max_lambda_plus_mu be???
  
  //u3::U3 sigma(0,{0u,0u});
  //int Nn_max=20;
  u3::U3 sigma(16,{2u,1u});
  //u3::U3 sigma(5,{2u,0u});
  int Nn_max=8;
  //int Nn_max=4;
  sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
  u3::UCoefCache u_coef_cache;
  int rho0 = 1;
  int outer_multiplicity = 1;
  u3::U3 omega_Id(0,{0u,0u});
  u3::U3 omega0(0,{0u,0u});
  u3::U3 omega1(0,{1u,1u});
  u3::U3 omega2(0,{2u,2u});
  // check these labels
  u3::U3 omega_A(2,{2u,0u});
  u3::U3 omega_B(-2,{0u,2u});
  u3::U3 omega_C(0,{1u,1u});
  std::cout<<sp3r_space.DebugStr()<<std::endl;

  for (const auto& bra_subspace : sp3r_space)
  {
    for (const auto& ket_subspace : sp3r_space)
    {
      const auto& bra_labels = bra_subspace.omega();
      const auto& ket_labels = ket_subspace.omega();
      // Skip highest N values in tests
      if (((bra_labels.N() - sigma.N()) == Nn_max) || ((ket_labels.N() - sigma.N()) == Nn_max)) continue;
      // compute rme of Identity coupled to itself
      std::cout<<fmt::format("Testing I x I")<<std::endl;
      basis::OperatorBlock<double> IxI = sp3r_operator::CoupledU3Operators(
        &sp3r_operator::IdentitySp3R,
        &sp3r_operator::IdentitySp3R,
        omega_Id,
        omega_Id,
        omega0,
        rho0,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> validation_id = sp3r_operator::IdentitySp3R(omega0,rho0,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(IxI-validation_id,1e-8))
      {
        std::cout<<fmt::format("I x I incorrect rme")<<std::endl;
        std::cout<<IxI<<std::endl;
        std::cout<<validation_id<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("I x I correct rme")<<std::endl;
        //std::cout<<IxI<<std::endl;
        //std::cout<<validation_id<<std::endl;
      }
      // compute rme of [B,A]^omega
      // TODO FIX DEFINITIONS AND CALCULATE EXPLICIT COUPLINGS
      std::cout<<fmt::format("Testing [B,A]^(0,0)")<<std::endl;
      basis::OperatorBlock<double> BA00 = sp3r_operator::CommutatorLoweringRaising(
        omega_B,
        omega_A,
        omega0,
        rho0,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> val_BA00 = std::sqrt(2.0/3.0)*ket_labels.N().TwiceValue()*sp3r_operator::IdentitySp3R(omega0,rho0,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(BA00-val_BA00,1e-8))
      {
        std::cout<<fmt::format("[B,A]^(0,0) incorrect rme")<<std::endl;
        std::cout<<BA00<<std::endl;
        std::cout<<val_BA00<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[B,A]^(0,0) correct rme")<<std::endl;
        //std::cout<<IxI<<std::endl;
        //std::cout<<validation_id<<std::endl;
      }
      // compute rme of [B,A]^omega
      // TODO FIX DEFINITIONS AND CALCULATE EXPLICIT COUPLINGS
      std::cout<<fmt::format("Testing [B,A]^(0,0)")<<std::endl;
      basis::OperatorBlock<double> comm_BA00 = sp3r_operator::CommutatorBA(
        omega0,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> validation_BA00 = std::sqrt(8.0/3.0)*sp3r_operator::HarmonicOscillatorSp3R(omega0,rho0,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(comm_BA00-validation_BA00,1e-8))
      {
        std::cout<<fmt::format("[B,A]^(0,0) incorrect rme")<<std::endl;
        std::cout<<comm_BA00<<std::endl;
        std::cout<<validation_BA00<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[B,A]^(0,0) correct rme")<<std::endl;
        //std::cout<<IxI<<std::endl;
        //std::cout<<validation_id<<std::endl;
      }
      std::cout<<fmt::format("Testing [B,A]^(1,1)")<<std::endl;
      basis::OperatorBlock<double> comm_BA11 = sp3r_operator::CommutatorBA(
        omega1,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> validation_BA11 = std::sqrt(5.0/2.0)*sp3r_operator::SU3GeneratorSp3R(omega1,rho0,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(comm_BA11-validation_BA11,1e-8))
      {
        std::cout<<fmt::format("[B,A]^(1,1) incorrect rme")<<std::endl;
        std::cout<<comm_BA11<<std::endl;
        std::cout<<validation_BA11<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[B,A]^(1,1) correct rme")<<std::endl;
        //std::cout<<comm_BA11<<std::endl;
        //std::cout<<validation_BA11<<std::endl;
      }
      std::cout<<fmt::format("Testing [B,A]^(2,2)")<<std::endl;
      basis::OperatorBlock<double> comm_BA22 = sp3r_operator::CommutatorBA(
        omega2,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> validation_BA22 = (0.0)*sp3r_operator::IdentitySp3R(omega2,rho0,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(comm_BA22-validation_BA22,1e-8))
      {
        std::cout<<fmt::format("[B,A]^(2,2) incorrect rme")<<std::endl;
        std::cout<<comm_BA22<<std::endl;
        std::cout<<validation_BA22<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[B,A]^(2,2) correct rme")<<std::endl;
        //std::cout<<comm_BA22<<std::endl;
        //std::cout<<validation_BA22<<std::endl;
      }
      // compute rme of [H,A]^omega
      // TODO FIX DEFINITIONS AND CALCULATE EXPLICIT COUPLINGS
      basis::OperatorBlock<double> validation_A = sp3r_operator::RaisingOperatorSp3R(omega_A,rho0,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      basis::OperatorBlock<double> validation_B = sp3r_operator::LoweringOperatorSp3R(omega_B,rho0,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      std::cout<<fmt::format("Testing [H_0,A]^(2,0)")<<std::endl;
      basis::OperatorBlock<double> comm_HA = sp3r_operator::CommutatorHA(
        omega_A,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(comm_HA-(2.0)*validation_A,1e-8))
      {
        std::cout<<fmt::format("[H_0,A]^(2,0) incorrect rme")<<std::endl;
        std::cout<<comm_HA<<std::endl;
        std::cout<<(2.0)*validation_A<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[H_0,A]^(0,0) correct rme")<<std::endl;
        //std::cout<<comm_HA<<std::endl;
        //std::cout<<(2.0)*validation_A<<std::endl;
      }
      std::cout<<fmt::format("Testing [[B,A]^(0,0),A]^(2,0)")<<std::endl;
      basis::OperatorBlock<double> comm_OscA = sp3r_operator::CommutatorOscillatorA(
        omega_A,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(comm_OscA-(2.0)*std::sqrt(8.0/3.0)*validation_A,1e-8))
      {
        std::cout<<fmt::format("[[B,A]^(0,0),A]^(2,0) incorrect rme")<<std::endl;
        std::cout<<comm_OscA<<std::endl;
        std::cout<<(2.0)*std::sqrt(8.0/3.0)*validation_A<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[[B,A]^(0,0),A]^(0,0) correct rme")<<std::endl;
        //std::cout<<comm_OscA<<std::endl;
        //std::cout<<(2.0)*std::sqrt(8.0/3.0)*validation_A<<std::endl;
      }
      std::cout<<fmt::format("Testing [H_0,B]^(0,2)")<<std::endl;
      basis::OperatorBlock<double> comm_HB = sp3r_operator::CommutatorHB(
        omega_B,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(comm_HB-(-2.0)*validation_B,1e-8))
      {
        std::cout<<fmt::format("[H_0,B]^(0,2) incorrect rme")<<std::endl;
        std::cout<<comm_HA<<std::endl;
        std::cout<<(-2.0)*validation_B<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[H_0,B]^(0,0) correct rme")<<std::endl;
        //std::cout<<comm_HA<<std::endl;
        //std::cout<<(-2.0)*validation_B<<std::endl;
      }
      std::cout<<fmt::format("Testing [[B,A]^(0,0),B]^(0,2)")<<std::endl;
      basis::OperatorBlock<double> comm_OscB = sp3r_operator::CommutatorOscillatorB(
        omega_B,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(comm_OscB-(-2.0)*std::sqrt(8.0/3.0)*validation_B,1e-8))
      {
        std::cout<<fmt::format("[[B,A]^(0,0),B]^(0,2) incorrect rme")<<std::endl;
        std::cout<<comm_OscB<<std::endl;
        std::cout<<(-2.0)*std::sqrt(8.0/3.0)*validation_B<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[[B,A]^(0,0),B]^(0,0) correct rme")<<std::endl;
        //std::cout<<comm_OscB<<std::endl;
        //std::cout<<(-2.0)*std::sqrt(8.0/3.0)*validation_B<<std::endl;
      }
      std::cout<<fmt::format("Testing [B,C]^(0,2)")<<std::endl;
      basis::OperatorBlock<double> comm_BC = sp3r_operator::CommutatorBC(
        omega_B,
        outer_multiplicity,
        sp3r_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(comm_BC-std::sqrt(40.0/3.0)*validation_B,1e-8))
      {
        std::cout<<fmt::format("[B,C]^(0,2) incorrect rme")<<std::endl;
        std::cout<<comm_BC<<std::endl;
        std::cout<<std::sqrt(1.0)*validation_B<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[H_0,B]^(0,0) correct rme")<<std::endl;
        //std::cout<<comm_HA<<std::endl;
        //std::cout<<(-2.0)*validation_B<<std::endl;
      }
    }
  }


  

  /*
  u3::U3 sigma2(16,{2u,1u});

  u3::U3 omegap(0,{3u,0u});
  u3::U3 omegam(0,{0u,3u});
  //u3::U3 sigma(5,{2u,0u});
  int Nn_max2=0;
  //int Nn_max=4;
  sp3r::Sp3RSpace sp3r_space2(sigma2,Nn_max2);
  
  outer_multiplicity = 1;

  for (const auto& bra_subspace : sp3r_space2)
    //for (const auto& ket_subspace : sp3r_space2)
    {
      const auto& bra_labels = bra_subspace.omega();
      //const auto& ket_labels = ket_subspace.omega();
      // Test CxC^(0,0)
      basis::OperatorBlock<double> CxC0 = sp3r_operator::CtimesCSp3R(
            omega0,
            outer_multiplicity,
            sp3r_space2,
            bra_subspace,
            bra_subspace,
            u_coef_cache
          );

      std::cout<<CxC0<<std::endl;

      // Test CxC^(1,1)
      basis::OperatorBlock<double> CxC1 = sp3r_operator::CtimesCSp3R(
            omega1,
            outer_multiplicity,
            sp3r_space2,
            bra_subspace,
            bra_subspace,
            u_coef_cache
          );

      std::cout<<CxC1<<std::endl;

      // Test CxC^(2,2)
      basis::OperatorBlock<double> CxC2 = sp3r_operator::CtimesCSp3R(
            omega2,
            outer_multiplicity,
            sp3r_space2,
            bra_subspace,
            bra_subspace,
            u_coef_cache
          );
      
      std::cout<<CxC2<<std::endl;

      // Test CxC^(2,2)
      basis::OperatorBlock<double> CxCp = sp3r_operator::CtimesCSp3R(
            omegap,
            outer_multiplicity,
            sp3r_space2,
            bra_subspace,
            bra_subspace,
            u_coef_cache
          );

      std::cout<<CxCp<<std::endl;

      // Test CxC^(2,2)
      basis::OperatorBlock<double> CxCm = sp3r_operator::CtimesCSp3R(
            omegam,
            outer_multiplicity,
            sp3r_space2,
            bra_subspace,
            bra_subspace,
            u_coef_cache
          );
      
      std::cout<<CxCm<<std::endl;
    }
  */
    

}
