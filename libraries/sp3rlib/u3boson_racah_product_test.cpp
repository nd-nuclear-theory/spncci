/****************************************************************
  u3boson_racah_product_test.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  08/18/23 (cvc): Created. Split from racah_product.cpp
  06/03/24 (cvc): Renamed u3boson_racah_product_test.cpp
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
#include "sp3rlib/u3.h"

#include "sp3rlib/u3boson_racah_product.h"


// Test couplings
int main(int argc, char **argv)
{
  u3::U3CoefInit(100); // what size should max_lambda_plus_mu be???
  
  //u3::U3 sigma(0,{0u,0u});
  //int Nn_max=20;
  u3::U3 sigma(16,{2u,1u});
  int Nn_max=8;
  u3boson::U3BosonSpace u3boson_space(sigma,Nn_max);
  u3::UCoefCache u_coef_cache;
  int rho0 = 1;
  int outer_multiplicity = 1;
  u3::U3 omega_Id(0,{0u,0u});
  u3::U3 omega0(0,{0u,0u});
  u3::U3 omega1(0,{1u,1u});
  u3::U3 omega2(0,{2u,2u});
  // check these labels
  u3::U3 omega_adag(2,{2u,0u});
  u3::U3 omega_a(-2,{0u,2u});
  u3::U3 omega_cas(0,{1u,1u});
  u3::U3 omega_adag_adag(4,{2u,1u});
  u3::U3 omega_a_a(-4,{1u,2u});
  std::cout<<u3boson_space.DebugStr()<<std::endl;

  // diagnostics
  int failed_IxI = 0;
  int failed_comm_aa = 0;
  int failed_Ixad = 0;
  int failed_adxI = 0;
  int failed_Ixa = 0;
  int failed_axI = 0;
  int tot_tests = 0;
  //std::cout<<u3boson_space.size()<<std::endl;
  for (const auto& bra_subspace : u3boson_space)
  {
    for (const auto& ket_subspace : u3boson_space)
    {
      const auto& bra_labels = bra_subspace.omega();
      const auto& ket_labels = ket_subspace.omega();
      tot_tests++;
      // Skip highest N values in tests
      if (((bra_labels.N() - sigma.N()) == Nn_max) || ((ket_labels.N() - sigma.N()) == Nn_max)) continue;
      //std::cout<<bra_subspace.DebugStr()<<std::endl;
      //std::cout<<ket_subspace.DebugStr()<<std::endl;
      // compute rme of Identity coupled to itself
      //std::cout<<fmt::format("Testing I x I")<<std::endl;
      basis::OperatorBlock<double> IxI = u3boson_operator::CoupledU3BosonOperators(
        &u3boson_operator::IdentityU3Boson,
        &u3boson_operator::IdentityU3Boson,
        omega_Id,
        omega_Id,
        omega0,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> validation_id = u3boson_operator::IdentityU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(IxI-validation_id,1e-8))
      {
        std::cout<<fmt::format("I x I incorrect rme")<<std::endl;
        std::cout<<IxI<<std::endl;
        std::cout<<validation_id<<std::endl;
        failed_IxI++;
      }
      else
      {
        //std::cout<<fmt::format("I x I correct rme")<<std::endl;
        //std::cout<<IxI<<std::endl;
        //std::cout<<validation_id<<std::endl;
      }
      // compute rme of [a,a^dagger], check if proportional to identity
      //std::cout<<fmt::format("Testing [a,a^dagger]")<<std::endl;
      basis::OperatorBlock<double> adaggera = u3boson_operator::CommutatorBosonQuanta(
        omega_a,
        omega_adag,
        omega0,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> validation_comm = std::sqrt(6)*u3boson_operator::IdentityU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(adaggera-validation_comm,1e-8))
      {
        std::cout<<fmt::format("[a,a^dagger] incorrect rme")<<std::endl;
        std::cout<<adaggera<<std::endl;
        std::cout<<validation_comm<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
        failed_comm_aa++;
      }
      else
      {
        //std::cout<<fmt::format("[a,a^dagger] correct rme")<<std::endl;
        //std::cout<<adaggerb<<std::endl;
        //std::cout<<validation_comm<<std::endl;
      }
      // compute rme of [a,a^dagger], check if proportional to identity
      //std::cout<<fmt::format("Testing [a,a^dagger]^(0,0)")<<std::endl;
      basis::OperatorBlock<double> aad00 = u3boson_operator::CommutatorAnnihilationCreation(
        omega0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      //basis::OperatorBlock<double> validation_comm = std::sqrt(6)*IdentityU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(aad00-validation_comm,1e-8))
      {
        std::cout<<fmt::format("[a,a^dagger]^(0,0) incorrect rme")<<std::endl;
        std::cout<<aad00<<std::endl;
        std::cout<<validation_comm<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
        failed_comm_aa++;
      }
      else
      {
        //std::cout<<fmt::format("[a,a^dagger] correct rme")<<std::endl;
        //std::cout<<adaggerb<<std::endl;
        //std::cout<<validation_comm<<std::endl;
      }
      //std::cout<<fmt::format("Testing [a,a^dagger]^(1,1)")<<std::endl;
      basis::OperatorBlock<double> aad11 = u3boson_operator::CommutatorAnnihilationCreation(
        omega1,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      //basis::OperatorBlock<double> validation_comm = std::sqrt(6)*IdentityU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(aad11-0.0*validation_comm,1e-8))
      {
        std::cout<<fmt::format("[a,a^dagger]^(1,1) incorrect rme")<<std::endl;
        std::cout<<aad11<<std::endl;
        std::cout<<0.0*validation_comm<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
        failed_comm_aa++;
      }
      else
      {
        //std::cout<<fmt::format("[a,a^dagger] correct rme")<<std::endl;
        //std::cout<<adaggerb<<std::endl;
        //std::cout<<validation_comm<<std::endl;
      }
      //std::cout<<fmt::format("Testing [a,a^dagger]^(2,2)")<<std::endl;
      basis::OperatorBlock<double> aad22 = u3boson_operator::CommutatorAnnihilationCreation(
        omega2,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      //basis::OperatorBlock<double> validation_comm = std::sqrt(6)*IdentityU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(aad22-0.0*validation_comm,1e-8))
      {
        std::cout<<fmt::format("[a,a^dagger]^(2,2) incorrect rme")<<std::endl;
        std::cout<<aad22<<std::endl;
        std::cout<<0.0*validation_comm<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
        failed_comm_aa++;
      }
      else
      {
        //std::cout<<fmt::format("[a,a^dagger] correct rme")<<std::endl;
        //std::cout<<adaggerb<<std::endl;
        //std::cout<<validation_comm<<std::endl;
      }
      // Check commutators with self [a,a]^(1,2)
      //std::cout<<fmt::format("Testing [a,a]^(1,2)")<<std::endl;
      basis::OperatorBlock<double> comm_a_a = u3boson_operator::CommutatorAnnihilationAnnihilation(
        omega_a_a,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      //basis::OperatorBlock<double> validation_comm = std::sqrt(6)*IdentityU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(comm_a_a-0.0*validation_comm,1e-8))
      {
        std::cout<<fmt::format("[a,a]^(1,2) incorrect rme")<<std::endl;
        std::cout<<comm_a_a<<std::endl;
        std::cout<<0.0*validation_comm<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
        failed_comm_aa++;
      }
      else
      {
        //std::cout<<fmt::format("[a,a]^(1,2) correct rme")<<std::endl;
        //std::cout<<comm_a_a<<std::endl;
        //std::cout<<0.0*validation_comm<<std::endl;
      }
      // Check commutators with self [a^\dagger,a^\dagger]^(2,1)
      //std::cout<<fmt::format("Testing [a^dagger,a^dagger]^(2,1)")<<std::endl;
      basis::OperatorBlock<double> comm_adag_adag = u3boson_operator::CommutatorCreationCreation(
        omega_adag_adag,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      //basis::OperatorBlock<double> validation_comm = std::sqrt(6)*IdentityU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(comm_adag_adag-0.0*validation_comm,1e-8))
      {
        std::cout<<fmt::format("[a^dagger,a^dagger]^(2,1) incorrect rme")<<std::endl;
        std::cout<<comm_adag_adag<<std::endl;
        std::cout<<0.0*validation_comm<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
        failed_comm_aa++;
      }
      else
      {
        //std::cout<<fmt::format("[a^dagger,a^dagger]^(2,1) correct rme")<<std::endl;
        //std::cout<<comm_adag_adag<<std::endl;
        //std::cout<<0.0*validation_comm<<std::endl;
      }
      // validation for a^dagger and a
      basis::OperatorBlock<double> validation_adag = u3boson_operator::CreationU3Boson(omega_adag,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      basis::OperatorBlock<double> validation_a = u3boson_operator::AnnihiliationU3Boson(omega_a,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      // compute rme of I x a^dagger
      //std::cout<<fmt::format("Testing I x a^dagger")<<std::endl;
      basis::OperatorBlock<double> Ixadag = u3boson_operator::CoupledU3BosonOperators(
        &u3boson_operator::IdentityU3Boson,
        &u3boson_operator::CreationU3Boson,
        omega_Id,
        omega_adag,
        omega_adag,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(Ixadag-validation_adag,1e-8))
      {
        std::cout<<fmt::format("I x a^dagger incorrect rme")<<std::endl;
        std::cout<<Ixadag<<std::endl;
        std::cout<<validation_adag<<std::endl;
        failed_Ixad++;
      }
      else
      {
        //std::cout<<fmt::format("I x I correct rme")<<std::endl;
        //std::cout<<Ixadag<<std::endl;
        //std::cout<<validation_adag<<std::endl;
      }
      // compute rme of a^dagger x I
      //std::cout<<fmt::format("Testing a^dagger x I")<<std::endl;
      basis::OperatorBlock<double> adagxI = u3boson_operator::CoupledU3BosonOperators(
        &u3boson_operator::CreationU3Boson,
        &u3boson_operator::IdentityU3Boson,
        omega_adag,
        omega_Id,
        omega_adag,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(adagxI-validation_adag,1e-8))
      {
        std::cout<<fmt::format("a^dagger x I incorrect rme")<<std::endl;
        std::cout<<adagxI<<std::endl;
        std::cout<<validation_adag<<std::endl;
        failed_adxI++;
      }
      else
      {
        //std::cout<<fmt::format("I x I correct rme")<<std::endl;
        //std::cout<<adagxI<<std::endl;
        //std::cout<<validation_adag<<std::endl;
      }
      // compute rme of I x a
      //std::cout<<fmt::format("Testing I x a")<<std::endl;
      basis::OperatorBlock<double> Ixa = u3boson_operator::CoupledU3BosonOperators(
        &u3boson_operator::IdentityU3Boson,
        &u3boson_operator::AnnihiliationU3Boson,
        omega_Id,
        omega_a,
        omega_a,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(Ixa-validation_a,1e-8))
      {
        std::cout<<fmt::format("I x a incorrect rme")<<std::endl;
        std::cout<<Ixa<<std::endl;
        std::cout<<validation_a<<std::endl;
        failed_Ixa++;
      }
      else
      {
        //std::cout<<fmt::format("I x I correct rme")<<std::endl;
        //std::cout<<Ixa<<std::endl;
        //std::cout<<validation_a<<std::endl;
      }
      // compute rme of a x I
      //std::cout<<fmt::format("Testing a x I")<<std::endl;
      basis::OperatorBlock<double> axI = u3boson_operator::CoupledU3BosonOperators(
        &u3boson_operator::AnnihiliationU3Boson,
        &u3boson_operator::IdentityU3Boson,
        omega_a,
        omega_Id,
        omega_a,
        rho0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(axI-validation_a,1e-8))
      {
        std::cout<<fmt::format("a x I incorrect rme")<<std::endl;
        std::cout<<axI<<std::endl;
        std::cout<<validation_a<<std::endl;
        failed_adxI++;
      }
      else
      {
        //std::cout<<fmt::format("I x I correct rme")<<std::endl;
        //std::cout<<axI<<std::endl;
        //std::cout<<validation_a<<std::endl;
      }
      // compute rme of [[a^dagger x a], a^dagger], check if proportional to a^dagger
      //std::cout<<fmt::format("Testing [[a^dagger x a], a^dagger]^(2,0)")<<std::endl;
      basis::OperatorBlock<double> nest_adag = u3boson_operator::CommutatorNumberCreation(
        omega_adag,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(nest_adag-std::sqrt(1.0/6.0)*validation_adag,1e-8))
      {
        std::cout<<fmt::format("[[a^dagger x a], a^dagger]^(2,0) incorrect rme")<<std::endl;
        std::cout<<nest_adag<<std::endl;
        std::cout<<std::sqrt(1.0/6.0)*validation_adag<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[a,a^dagger] correct rme")<<std::endl;
        //std::cout<<adaggerb<<std::endl;
        //std::cout<<validation_comm<<std::endl;
      }
      // compute rme of [[a^dagger x a], a], check if proportional to a
      //std::cout<<fmt::format("Testing [[a^dagger x a], a]^(0,2)")<<std::endl;
      basis::OperatorBlock<double> nest_a = u3boson_operator::CommutatorNumberAnnihilation(
        omega_a,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      if (!mcutils::IsZero(nest_a-(-1.0)*std::sqrt(1.0/6.0)*validation_a,1e-8))
      {
        std::cout<<fmt::format("[[a^dagger x a], a]^(0,2) incorrect rme")<<std::endl;
        std::cout<<nest_a<<std::endl;
        std::cout<<(-1.0)*std::sqrt(1.0/6.0)*validation_a<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[a,a^dagger] correct rme")<<std::endl;
        //std::cout<<adaggerb<<std::endl;
        //std::cout<<validation_comm<<std::endl;
      }
      // compute rme of [a^dagger x a], check if proportional to N*I
      //std::cout<<fmt::format("Testing [a^dagger x a]^(0,0)")<<std::endl;
      basis::OperatorBlock<double> number_op = u3boson_operator::CreationAnnihilation(
        omega0,
        outer_multiplicity,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      int boson_excitations = (ket_labels.N().TwiceValue()-sigma.N().TwiceValue())/2;
      if (!mcutils::IsZero(number_op-(boson_excitations/2.0)*std::sqrt(1.0/6.0)*validation_id,1e-8))
      {
        std::cout<<fmt::format("[a^dagger x a]^(0,0) incorrect rme")<<std::endl;
        std::cout<<number_op<<std::endl;
        std::cout<<(boson_excitations/2.0)*std::sqrt(1.0/6.0)*validation_id<<std::endl;
        //std::cout<<bra_subspace.DebugStr()<<std::endl;
        //std::cout<<ket_subspace.DebugStr()<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("[a,a^dagger] correct rme")<<std::endl;
        //std::cout<<adaggerb<<std::endl;
        //std::cout<<validation_comm<<std::endl;
      }
      /*
      // compute rme of SU(3) Casimir operator
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      ////// CHECK THIS AND POTENTIALLY TRY SOME SORT OF CASIMIR COUPLING //////
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      std::cout<<fmt::format("Testing SU(3) Casimir")<<std::endl;
      basis::OperatorBlock<double> su3casimir = u3boson_operator::SU3CasimirCoupling(
        omega0,
        rho0,
        u3boson_space,
        bra_subspace,
        ket_subspace,
        u_coef_cache
      );
      basis::OperatorBlock<double> validation_cas = u3boson_operator::CasimirU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache);
      if (!mcutils::IsZero(su3casimir-validation_cas,1e-8))
      {
        std::cout<<fmt::format("SU(3) Casimir incorrect rme")<<std::endl;
        std::cout<<su3casimir<<std::endl;
        std::cout<<CasimirU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache)<<std::endl;
        std::cout<<bra_subspace.DebugStr()<<std::endl;
        std::cout<<ket_subspace.DebugStr()<<std::endl;
      }
      else
      {
        //std::cout<<fmt::format("SU(3) Casimir correct rme")<<std::endl;
        std::cout<<su3casimir<<std::endl;
        std::cout<<CasimirU3Boson(omega0,rho0,u3boson_space,bra_subspace,ket_subspace,u_coef_cache)<<std::endl;
      }
      */

    }
  }
}
