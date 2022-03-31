/****************************************************************
  sp3r_gtest.cpp

  Anna E. McCoy
  Institute for Nuclear Theory
  
  SPDX-License-Identifier: MIT 

  3/30/22 (aem): Created.
****************************************************************/
#include "sp3rlib/sp3r.h"
#include "gtest/gtest.h"

#include <iostream>
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "fmt/format.h"

class Sp3RSpaceTest
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
    Sp3RTest,
    Sp3RSpaceTest,
    testing::Values(
      std::pair{u3::U3(16,{2u,1u}), 6},
      std::pair{u3::U3({29,2},{2u,1u}), 6},
      std::pair{u3::U3({3,2},{0u,0u}),6},
      std::pair{u3::U3({5,{2u,0u}}),6}
    )
  );


std::map<u3::U3,bool> modify_basis_test_map =
  {
    {{16,{2u,1u}},false},
    {{{29,2},{2u,1u}},false},
    {{{3,2},{0u,0u}},true},
    {{5,{2u,0u}},true}
  };

// std::map<u3::U3,bool> unitary_test_map=
//   {
//     {{16,{2u,1u}},true},
//     {{3,{0u,0u}},true},
//     {{6,5,1},false}
//   };


void CheckingU3Subspace(
  const u3::U3& sigma,
  const unsigned int Nn_max,
  const sp3r::U3Subspace& subspace
  )
  {
    const auto& omega = subspace.omega();
    // Check SO3States
    const auto& L_kappa_vector = u3::BranchingSO3(omega.SU3());

    // All of the possible L values hould be included in the set of SO3States
    for(const auto&[L,kappa_max] : L_kappa_vector)
      EXPECT_TRUE(subspace.LookUpStateIndex(L) != basis::kNone)
        <<fmt::format("L = {} not in omega = {} subspace",L,omega);

    // Check non-orthogonal basis.  All raising polynomials which can couple
    // with sigma to omega should be included in the basis.
    const auto raising_polynomials = u3boson::RaisingPolynomialLabels(Nn_max);
    const auto& nonorthogonal_basis = subspace.nonorthogonal_basis();
    int num_state_labels=0;
    int num_states=0;
    for(const auto& n : raising_polynomials)
    {
      unsigned int rho_max = u3::OuterMultiplicity(sigma,n,omega);
      if(rho_max>0)
        {
          EXPECT_TRUE(nonorthogonal_basis.LookUpStateIndex(n) != basis::kNone)
          <<fmt::format("Raising polynomial {} not in nonorthogonal_basis for sigma = {}, omega = {}",
            n,sigma,omega
          );
          num_state_labels++;
          num_states+=rho_max;
        }
    }

    EXPECT_TRUE(num_state_labels == nonorthogonal_basis.size());
    EXPECT_TRUE(num_states == nonorthogonal_basis.dimension());
    EXPECT_TRUE(subspace.nonorthogonal_basis().dimension()==subspace.K_matrix().cols())
      <<fmt::format("Kmatrix columns do not match nonorthogonal_basis dimension: {} vs {}",
          subspace.nonorthogonal_basis().dimension(),subspace.K_matrix().cols()
        );

    EXPECT_TRUE(subspace.nonorthogonal_basis_dimension()==subspace.nonorthogonal_basis().dimension());

    EXPECT_TRUE(subspace.upsilon_max()==subspace.K_matrix().rows())
    <<fmt::format("Kmatrix rows do not match upsilon_max: {} vs {}",
        subspace.upsilon_max(),subspace.K_matrix().rows()
      );

    EXPECT_TRUE(subspace.omega()==subspace.U3());

    EXPECT_TRUE(subspace.upsilon_max()==subspace.Kinv_matrix().cols());
    EXPECT_TRUE(subspace.nonorthogonal_basis_dimension()==subspace.K_matrix().cols());
    EXPECT_TRUE(subspace.nonorthogonal_basis_dimension()==subspace.Kinv_matrix().rows());
  }

void CheckingSp3RSpace(
  const u3::U3& sigma,
  const unsigned int Nn_max,
  const sp3r::Sp3RSpace& irrep
  )
  {
    // Checking accessors
    EXPECT_TRUE(irrep.sigma()==sigma)<<fmt::format("sigma does not match: {} vs {}",irrep.sigma(),sigma);
    EXPECT_TRUE(irrep.Nn_max()==Nn_max)<<fmt::format("Nn_max does not match: {} vs {}",irrep.Nn_max(),Nn_max);

    EXPECT_TRUE(modify_basis_test_map[sigma]==sp3r::ModifySp3RBranching(sigma))
    <<fmt::format("Irrep {} needs to be modfied? {} vs ",
      sigma,modify_basis_test_map[sigma],sp3r::ModifySp3RBranching(sigma)
    );

    // Checking subspaces
    for(int i=0; i<irrep.size(); ++i)
      {
        const auto& subspace = irrep.GetSubspace(i);
        CheckingU3Subspace(sigma, Nn_max,subspace);
      }
  }

TEST_P(Sp3RSpaceTest, BasisContent)
{
  const auto& [sigma, Nn_max] = GetParam();
  sp3r::Sp3RSpace irrep(sigma,Nn_max);
  CheckingSp3RSpace(sigma,Nn_max,irrep);
}


