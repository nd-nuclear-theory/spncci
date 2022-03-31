/****************************************************************
  vcs_gtest.cpp

  Unit test for vcs module

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT 
 
   3/30/22 (aem) : Created.
****************************************************************/
#include "cppitertools/itertools.hpp"
#include "gtest/gtest.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "sp3rlib/vcs.h"




namespace vcs
{
	using Matrix =  basis::OperatorBlock<double>;

// Test Omega function for irreps
// 12(1,1) and 17/2(2,1) up to Nn_max=4
std::map<std::tuple<u3::U3,u3::U3>,double>
Omega_validation_table ={
    {{u3::U3(0,0,0),u3::U3(5,4,3)},27},
    {{u3::U3(2,0,0),u3::U3(5,5,4)},30},
    {{u3::U3(2,0,0),u3::U3(6,5,3)},34},
    {{u3::U3(2,0,0),u3::U3(6,4,4)},32},
    {{u3::U3(2,0,0),u3::U3(7,4,3)},37},
    {{u3::U3(2,2,0),u3::U3(6,6,4)},39},
    {{u3::U3(2,2,0),u3::U3(6,5,5)},37},
    {{u3::U3(2,2,0),u3::U3(7,6,3)},44},
    {{u3::U3(2,2,0),u3::U3(7,5,4)},41},
    {{u3::U3(4,0,0),u3::U3(7,5,4)},38},
    {{u3::U3(4,0,0),u3::U3(8,5,3)},44},
    {{u3::U3(4,0,0),u3::U3(8,4,4)},42},
    {{u3::U3(4,0,0),u3::U3(9,4,3)},49},
    {{u3::U3(0,0,0),u3::U3(4.5_hi,2.5_hi,1.5_hi)},17.375},
    {{u3::U3(2,0,0),u3::U3(4.5_hi,4.5_hi,1.5_hi)},20.375},
    {{u3::U3(2,0,0),u3::U3(4.5_hi,3.5_hi,2.5_hi)},17.375},
    {{u3::U3(2,0,0),u3::U3(5.5_hi,3.5_hi,1.5_hi)},22.375},
    {{u3::U3(2,0,0),u3::U3(5.5_hi,2.5_hi,2.5_hi)},20.375},
    {{u3::U3(2,0,0),u3::U3(6.5_hi,2.5_hi,1.5_hi)},26.375},
    {{u3::U3(2,2,0),u3::U3(4.5_hi,4.5_hi,3.5_hi)},20.375},
    {{u3::U3(2,2,0),u3::U3(5.5_hi,4.5_hi,2.5_hi)},24.375},
    {{u3::U3(4,0,0),u3::U3(5.5_hi,4.5_hi,2.5_hi)},21.375},
    {{u3::U3(2,2,0),u3::U3(5.5_hi,3.5_hi,3.5_hi)},22.375},
    {{u3::U3(2,2,0),u3::U3(6.5_hi,4.5_hi,1.5_hi)},30.375},
    {{u3::U3(4,0,0),u3::U3(6.5_hi,4.5_hi,1.5_hi)},27.375},
    {{u3::U3(2,2,0),u3::U3(6.5_hi,3.5_hi,2.5_hi)},27.375},
    {{u3::U3(4,0,0),u3::U3(6.5_hi,3.5_hi,2.5_hi)},24.375},
    {{u3::U3(4,0,0),u3::U3(7.5_hi,3.5_hi,1.5_hi)},31.375},
    {{u3::U3(4,0,0),u3::U3(7.5_hi,2.5_hi,2.5_hi)},29.375},
    {{u3::U3(4,0,0),u3::U3(8.5_hi,2.5_hi,1.5_hi)},37.375}
  };

std::map<u3::U3,std::map<u3::U3,vcs::Matrix>> Kmatrix_validation_table = {
{
  u3::U3(14,{3u,1u}),
  {
    {u3::U3(14,{3u,1u}),vcs::Matrix{{ 1.00000000}}},
    {u3::U3(16,{3u,2u}),vcs::Matrix{{ 3.00000000}}},
    {u3::U3(16,{4u,0u}),vcs::Matrix{{ 2.64575131}}},
    {u3::U3(18,{2u,2u}),vcs::Matrix{{ 6.80213819,-0.25347458},{-0.25347458,  4.96008239}}},
    {u3::U3(20,{1u,2u}),vcs::Matrix{{12.72792206,-0.00000000},{-0.00000000, 12.72792206}}},
    {u3::U3(20,{3u,4u}),vcs::Matrix{{33.82805696,-0.68011935},{-0.68011935, 25.58783769}}},
    {u3::U3(20,{2u,3u}),vcs::Matrix{{19.62822980, 2.16141450,-0.24674332},{ 2.16141450, 19.15191779, -0.72961120},{-0.24674332,-0.72961120, 13.68965979}}},
    {u3::U3(20,{3u,1u}),vcs::Matrix{{20.33536264,-1.99295318, 0.03411572},{-1.99295318, 16.44081762, -0.63570756},{ 0.03411572,-0.63570756, 15.37589807}}}
  }
},
{
  u3::U3(HalfInt(29,2),{2u,1u}),
  {
    {u3::U3(HalfInt(29,2),{2u,1u}),vcs::Matrix{{ 1.00000000}}},
    {u3::U3(HalfInt(33,2),{0u,3u}),vcs::Matrix{{ 2.64575131}}},
    {u3::U3(HalfInt(37,2),{1u,2u}),vcs::Matrix{{ 7.39228476,-0.14419682},{-0.14419682,5.68441207}}},
    {u3::U3(HalfInt(37,2),{3u,1u}),vcs::Matrix{{ 8.20961975,-0.51846909},{-0.51846909,7.00460728}}},
    {u3::U3(HalfInt(41,2),{2u,1u}),vcs::Matrix{{23.07368495,-1.24931673, 0.21040579},{-1.24931673,18.44600494,-0.42907995},{0.21040579,-0.42907995,17.93799375}}}
  }
}
};

void CheckOmega(const u3::U3& sigma, const unsigned int Nn_max)
{

  if(vcs::Omega_validation_table.count({{0,0,0},sigma}))
  {
    u3boson::U3BosonSpace u3boson_space(sigma,Nn_max);
    for(const auto& subspace : u3boson_space)
      {
        const auto& omega=subspace.omega();
        for(const auto& state : subspace)
          {
            const auto& n = state.n();
            double Omega_value = vcs::Omega(n,omega);
            const double test_value = vcs::Omega_validation_table.at({n,omega});
            ASSERT_EQ(test_value, Omega_value)
            <<fmt::format("Omega value does not match for n = {} and omega = {}\n Expected: {}  Calculated: {}",
                n,omega,Omega_value,vcs::Omega_validation_table.at({n,omega})
              );
          }
      }
  }
}

void KMatrixTest(const u3::U3& sigma, const unsigned int Nn_max)
{
  // Compute K matrices
  u3boson::U3BosonSpace u3boson_space(sigma,Nn_max);
  auto K_matrix_map = vcs::GenerateKmatrices(sigma,u3boson_space,1e-12);

  // Get test values
  const auto& test_values = vcs::Kmatrix_validation_table[sigma];

  for(const auto& subspace : u3boson_space)
    {
      const auto& omega = subspace.omega();

      if(K_matrix_map.count(omega)==0)
      {
        EXPECT_TRUE(sp3r::ModifySp3RBranching(sigma))
        <<fmt::format("Omega = {} not in irrep {}",omega,sigma);
        continue;
      }

      const vcs::Matrix K = K_matrix_map.at(omega)[0];
      const vcs::Matrix Kinv = K_matrix_map.at(omega)[1];

      vcs::Matrix Id = vcs::Matrix::Identity(K.rows(),Kinv.cols());

      // Check Kmatrix dimensions
      EXPECT_TRUE(K.cols()==subspace.dimension());
      EXPECT_TRUE(Kinv.rows()==subspace.dimension());
      EXPECT_TRUE(K.cols()>=K.rows());
      EXPECT_TRUE(Kinv.rows()>=Kinv.cols());

      // Check Kinv is right inverse of K.
      EXPECT_TRUE(mcutils::IsZero(K*Kinv-Id,1e-10))
        <<fmt::format("K*Kinv!=I\n")<< K*Kinv-Id;

      // If omega in test_values, check against stored values
      if(test_values.count(omega))
        {
          EXPECT_TRUE(mcutils::IsZero(K-test_values.at(omega),1e-8))
            <<fmt::format("Calculated K matrix for {} does not match test value\n",omega)
            <<"Test:"<<std::endl
            <<test_values.at(omega)<<std::endl
            <<"Calculated:"<<std::endl
            <<K;
        }
    }
}



}//namespace

class VCSTestFixture
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
    VCSTestValues,
    VCSTestFixture,
    testing::Values(
      std::pair{u3::U3(5,4,3), 4},
      std::pair{u3::U3(4.5_hi,2.5_hi,1.5_hi), 4},
      std::pair{u3::U3(14,{3u,1u}), 10},
      std::pair{u3::U3(HalfInt(29,2),{2u,1u}), 10}
    )
  );

TEST_P(VCSTestFixture, VCSTest)
{
  const auto& [sigma, Nn_max] = GetParam();
  // vcs::CheckOmega(sigma,Nn_max);
  vcs::KMatrixTest(sigma,Nn_max);
}
