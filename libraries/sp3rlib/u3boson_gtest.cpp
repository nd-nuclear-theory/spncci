/****************************************************************
  u3boson_gtest.cpp

  Anna E. McCoy
  Institute for Nuclear Theory
  
  SPDX-License-Identifier: MIT 

  3/31/22 (aem): Created.
****************************************************************/
#include <iostream>

#include "gtest/gtest.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3boson.h"

namespace u3boson
{
  using Matrix = basis::OperatorBlock<double>;


void CheckingU3Subspace(
  const u3::U3& sigma,
  const unsigned int Nn_max,
  const u3boson::U3Subspace& subspace
  )
  {
    const auto& omega = subspace.omega();

    // All raising polynomials which can couple
    // with sigma to omega should be included in the basis.
    const auto raising_polynomials = u3boson::RaisingPolynomialLabels(Nn_max);

    int num_state_labels=0;
    int num_states=0;
    for(const auto& n : raising_polynomials)
    {
      unsigned int rho_max = u3::OuterMultiplicity(sigma,n,omega);
      if(rho_max>0)
        {
          EXPECT_TRUE(subspace.LookUpStateIndex(n) != basis::kNone)
          <<fmt::format("Raising polynomial {} not in basis for sigma = {}, omega = {}",
            n,sigma,omega
          );
          num_state_labels++;
          num_states+=rho_max;
        }
    }

    EXPECT_TRUE(num_state_labels == subspace.size());
    EXPECT_TRUE(num_states == subspace.dimension());
  }


void CheckingU3BosonSpace(
  const u3::U3& sigma,
  const unsigned int Nn_max,
  const u3boson::U3BosonSpace& irrep
  )
  {
    // Checking accessors
    EXPECT_TRUE(irrep.sigma()==sigma)<<fmt::format("sigma does not match: {} vs {}",irrep.sigma(),sigma);
    EXPECT_TRUE(irrep.Nn_max()==Nn_max)<<fmt::format("Nn_max does not match: {} vs {}",irrep.Nn_max(),Nn_max);

    // Checking subspaces
    for(int i=0; i<irrep.size(); ++i)
      {
        const auto& subspace = irrep.GetSubspace(i);
        CheckingU3Subspace(sigma, Nn_max,subspace);
      }
  }

void CheckU3BosonCommutator(
    const u3boson::U3BosonSpace& boson_space,
    u3::UCoefCache& u_coef_cache
  )
  // Check commutation relations given in jmp-24-1983-2461-Rosenteel
  // Eq.3
  {
    unsigned int Nn_max = boson_space.Nn_max();
    const u3::U3& sigma = boson_space.sigma();

    for(const auto& bra_subspace : boson_space)
      for(const auto& ket_subspace : boson_space)
        {
          const auto& omega_bra = bra_subspace.omega();
          const auto& omega_ket = ket_subspace.omega();
          if(omega_bra!=omega_ket) continue;
          if((omega_bra.N()-sigma.N())==Nn_max) continue;

          u3boson::Matrix validation_matrix
            = std::sqrt(6)*u3boson::Matrix::Identity(bra_subspace.dimension(), ket_subspace.dimension());

          u3boson::Matrix operator_block
            = u3boson::Matrix::Zero(bra_subspace.dimension(),ket_subspace.dimension());

          for(const auto& bar_subspace : boson_space)
            {
              const auto& omega_bar = bar_subspace.omega();
              // a x a^\dagger term
              if(u3::OuterMultiplicity(omega_bar,{-2,{0u,2u}},omega_bra)>0)
                if(u3::OuterMultiplicity(omega_ket,{2,{2u,0u}},omega_bar)>0)
                  {
                    operator_block
                      += u3::UCached(u_coef_cache,omega_ket.SU3(),{2u,0u},omega_bra.SU3(),{0u,2u},omega_bar.SU3(),1,1,{0u,0u},1,1)
                          *u3boson::U3BosonAnnihilationOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                          *u3boson::U3BosonCreationOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                  }

              // a^\dagger x a term
              if(u3::OuterMultiplicity(omega_ket,{-2,{0u,2u}},omega_bar)>0)
                if(u3::OuterMultiplicity(omega_bar,{2,{2u,0u}},omega_bra)>0)
                  {
                    operator_block
                      -= u3::UCached(u_coef_cache,omega_ket.SU3(),{0u,2u},omega_bra.SU3(),{2u,0u},omega_bar.SU3(),1,1,{0u,0u},1,1)
                          *u3boson::U3BosonCreationOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                          *u3boson::U3BosonAnnihilationOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                  }

            }
          EXPECT_TRUE(mcutils::IsZero(operator_block-validation_matrix,1e-8));
        }
  }

}

class U3BosonTest
    : public testing::TestWithParam<std::pair<u3::U3, unsigned int>>
{
 protected:
  void SetUp() override
  {
    // initialization
    auto& [sigma, Nn_max] = GetParam();
    const auto lambda_plus_mu = (sigma.SU3().lambda() + sigma.SU3().mu() + 2*Nn_max);
    fmt::print("lambda+mu max: {:d}\n", lambda_plus_mu);
    u3::U3CoefInit(lambda_plus_mu);
  }

  void TearDown() override
  {
    // finalization -- TODO
    // u3::U3CoefFinalize();
  }
};


INSTANTIATE_TEST_SUITE_P(
    U3BosonTestValues,
    U3BosonTest,
    testing::Values(
      std::pair{u3::U3(16,{2u,1u}), 6},
      std::pair{u3::U3({29,2},{2u,1u}), 6},
      std::pair{u3::U3({3,2},{0u,0u}),6},
      std::pair{u3::U3({5,{2u,0u}}),6}
    )
  );

TEST_P(U3BosonTest, BasisContent)
{
  const auto& [sigma, Nn_max] = GetParam();
  u3boson::U3BosonSpace boson_space(sigma,Nn_max);
  u3boson::CheckingU3BosonSpace(sigma,Nn_max,boson_space);
}


//Validating coupled commutator of boson creation and annihilation operators
//Commutator given by eq. 3 of G. Rosensteel and D. J. Rowe. J. Math Phys. 24 (1983) 2461.
TEST_P(U3BosonTest, OperatorCommutator)
{
  const auto& [sigma, Nn_max] = GetParam();
  u3boson::U3BosonSpace boson_space(sigma,Nn_max);
  u3::UCoefCache u_coef_cache;
  u3boson::CheckU3BosonCommutator(boson_space,u_coef_cache);
}


// Comparison of RME value.  Test values from D. J. Rowe. J. Math. Phys. 27 (1984) 2662.
std::map<std::tuple<u3::U3,u3::U3>,u3boson::Matrix> u3boson_rmes_test_map={
    {{{22,13,10},{20,13,10}},u3boson::Matrix{{1}}},
    {{{21,14,10},{20,13,10}},u3boson::Matrix{{-1}}},
    {{{22,15,10},{22,13,10}},u3boson::Matrix{{std::sqrt(5./6)},{std::sqrt(7./36)}}},
    {{{22,15,10},{21,14,10}},u3boson::Matrix{{std::sqrt(2./3)},{-std::sqrt(80./63)}}}
  };

TEST(U3BosonRME,OperatorValues)
{
  u3::U3 sigma(20,13,10);
  unsigned int Nn_max =4;
  const auto lambda_plus_mu = (sigma.SU3().lambda() + sigma.SU3().mu() + 2*Nn_max);
  fmt::print("lambda+mu max: {:d}\n", lambda_plus_mu);
  u3::U3CoefInit(lambda_plus_mu);

  u3boson::U3BosonSpace boson_space(sigma,Nn_max);
  u3::UCoefCache u_coef_cache;
  for(const auto& [irreps,test_block] : u3boson_rmes_test_map)
    {
      const auto& [omega_bra,omega_ket] = irreps;
      const auto& bra_subspace = boson_space.LookUpSubspace(omega_bra);
      const auto& ket_subspace = boson_space.LookUpSubspace(omega_ket);
      u3boson::Matrix block
        = u3boson::U3BosonCreationOperator(sigma,bra_subspace,ket_subspace,u_coef_cache);

      EXPECT_TRUE(mcutils::IsZero(block-test_block,1e-8));
    }
}

