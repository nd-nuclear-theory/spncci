/****************************************************************
  su3coupling_gtest.cpp

  SU(3) coupling coefficient unit testing

  10/31/21 (cvc) Created.

****************************************************************/

#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <tuple>
#include <vector>

#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"
#include "gtest/gtest.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
// #include "utilities/utilities.h"

namespace u3
{
void PrintTo(const SU3& x, std::ostream* os)
{
  *os << x.Str();
}
};  // namespace u3


// Test fixture for W coefficients
class WCoeffTest
    : public testing::TestWithParam<std::pair<u3::SU3, u3::SU3>>
{
 protected:
  void SetUp() override
  {
    // initialization
    auto& [x1, x2] = GetParam();
    const auto lambda_plus_mu = (x1.lambda() + x1.mu() + x2.lambda() + x2.mu());
    fmt::print("lambda+mu max: {:d}\n", lambda_plus_mu);
    u3::U3CoefInit(lambda_plus_mu);
  }

  void TearDown() override
  {
    // finalization -- TODO
    // u3::U3CoefFinalize();
  }
};

// Set up test suite
INSTANTIATE_TEST_SUITE_P(
    SU3Test,
    WCoeffTest,
    testing::Values(
        std::pair{u3::SU3{1, 2}, u3::SU3{3, 4}},
        std::pair{u3::SU3{2, 0}, u3::SU3{2, 0}},
        std::pair{u3::SU3{4, 3}, u3::SU3{2, 4}}
      )
  );


/************************************************/
/*************** Helper functions ***************/
/************************************************/

// struct to hold x, L, and kappa
struct L_kappa
{
  // constructor
  L_kappa(const u3::SU3& x, int L_, int kappa_)
      : irrep(x), L(L_), kappa(kappa_)
  {}
  // default constructor
  L_kappa()
      : L(0), kappa(0)
  {}

  u3::SU3 irrep;
  int L;
  int kappa;
};


std::vector<L_kappa> get_coupling_labels(u3::SU3& x)
// get vector of allowed L, kappa values for coupling
{
  // branching
  MultiplicityTagged<unsigned int>::vector branch = u3::BranchingSO3(x);
  // vector containing labels
  std::vector<L_kappa> lk;
  // add L, kappa to vector of labels
  for (int i = 0; i < branch.size(); i++)
    for (int j = 1; j <= branch[i].tag; j++)
    {
      lk.push_back(L_kappa(x, branch[i].irrep, j));
    }
  return lk;
}


double compute_W_sum_SU3(
    u3::SU3& x1,
    u3::SU3& x2,
    int kappa1,
    int kappa1prime,
    int L1,
    int L1prime,
    int kappa2,
    int kappa2prime,
    int L2,
    int L2prime,
    int L
  )
// compute orthogonality sum over (lambda,mu), rho, kappa, for a given L
{
  double sum = 0;
  // couple x1 to x2
  MultiplicityTagged<u3::SU3>::vector coupling = u3::KroneckerProduct(x1, x2);
  for (int u = 0; u < coupling.size(); u++)
  {
    //(lambda,mu)
    u3::SU3 x = coupling[u].irrep;
    // outer multiplicity
    int rho_max = u3::OuterMultiplicity(x1, x2, x);
    std::vector<L_kappa> lk = get_coupling_labels(x);
    // sum over rho, L, kappa

    for (int rho = 1; rho <= rho_max; rho++)
    {
      for (L_kappa& tLK : lk)
      {
        if (L == tLK.L)  // add to sum otherwise skip
        {
          sum += u3::W(x1, kappa1, L1, x2, kappa2, L2, x, tLK.kappa, tLK.L, rho)
                 * u3::W(
                     x1,
                     kappa1prime,
                     L1prime,
                     x2,
                     kappa2prime,
                     L2prime,
                     x,
                     tLK.kappa,
                     tLK.L,
                     rho
                   );
        }
      }
    }
  }
  return sum;
}


double compute_W_sum_alpha(
    u3::SU3& x1,
    u3::SU3& x2,
    u3::SU3& x,
    u3::SU3& xprime,
    int kappa,
    int L,
    int rho,
    int rhoprime
  )
// compute orthogonality sum over L1, kappa1, L2, kappa2 for given kappa, L
{
  double sum = 0;
  std::vector<L_kappa> lk1 = get_coupling_labels(x1);
  std::vector<L_kappa> lk2 = get_coupling_labels(x2);
  for (L_kappa& tLK1 : lk1)
  {
    for (L_kappa& tLK2 : lk2)
    {
      double coef =
          u3::W(x1, tLK1.kappa, tLK1.L, x2, tLK2.kappa, tLK2.L, x, kappa, L, rho)
          * u3::W(
              x1, tLK1.kappa, tLK1.L, x2, tLK2.kappa, tLK2.L, xprime, kappa, L, rhoprime
            );
      EXPECT_TRUE(std::isfinite(coef)) << fmt::format(
          "coefficient not finite: W({} {} {}; {} {} {} || {} {} {}){} = {}",
          x1.Str(),
          tLK1.kappa,
          tLK1.L,
          x2.Str(),
          tLK2.kappa,
          tLK2.L,
          xprime.Str(),
          kappa,
          L,
          rhoprime,
          coef
        );
      sum += coef;
    }
  }
  return sum;
}

/*
// SU3 test class
class SU3TestClass : public ::testing::Test
{
protected:

    SU3TestClass() : su3_1(1,2), su3_2(3,4) {}
    ~SU3TestClass() override {}

    void SetUp() override
    {
        // initialization
        u3::U3CoefInit();
    }

    void TearDown() override {}

    u3::SU3 su3_1;
    u3::SU3 su3_2;
};
*/


/************************************************/
/****************** Unit Tests ******************/
/************************************************/


TEST_P(WCoeffTest, WvsWCached)
// check coupling (lambda_1, mu_1) x (lambda_2, mu_2)
// check that W = WCached
{
  auto [x1, x2] = GetParam();
  // vector of tuples of L, kappa
  std::vector<L_kappa> lk1 = get_coupling_labels(x1);
  std::vector<L_kappa> lk2 = get_coupling_labels(x2);
  // couple x1 to x2
  MultiplicityTagged<u3::SU3>::vector coupling = u3::KroneckerProduct(x1, x2);
  // check coupling for x3 in vector coupling
  for (int u = 0; u < coupling.size(); u++)
  {
    u3::SU3 x3 = coupling[u].irrep;
    // outer multiplicity
    int rho = u3::OuterMultiplicity(x1, x2, x3);
    // branching and vector of tuples of L, kappa
    std::vector<L_kappa> lk3 = get_coupling_labels(x3);
    // compute coefficients
    for (L_kappa& tupleLK1 : lk1)
      for (L_kappa& tupleLK2 : lk2)
        for (L_kappa& tupleLK3 : lk3)
        {
          int L1 = tupleLK1.L;
          int k1 = tupleLK1.kappa;
          int L2 = tupleLK2.L;
          int k2 = tupleLK2.kappa;
          int L3 = tupleLK3.L;
          int k3 = tupleLK3.kappa;

          // check triangle inequality
          if (L3 >= std::abs(L1 - L2) && L3 <= (L1 + L2))
          {
            // calculate coefficients
            double w_coeff = u3::W(x1, k1, L1, x2, k2, L2, x3, k3, L3, rho);
            ASSERT_TRUE(std::isfinite(w_coeff)) << fmt::format(
                "coefficient not finite: W({} {} {}; {} {} {} || {} {} {}){} = {}",
                x1.Str(),
                k1,
                L1,
                x2.Str(),
                k2,
                L2,
                x3.Str(),
                k3,
                L3,
                rho,
                w_coeff
              );
            // define empty cache
            u3::WCoefCache cache;
            double w_coeff_cached =
                u3::WCached(cache, x1, k1, L1, x2, k2, L2, x3, k3, L3, rho);
            EXPECT_EQ(w_coeff, w_coeff_cached) << fmt::format(
                "W({} {} {}; {} {} {} || {} {} {}){} = {} != WCached = ",
                x1.Str(),
                k1,
                L1,
                x2.Str(),
                k2,
                L2,
                x3.Str(),
                k3,
                L3,
                rho,
                w_coeff,
                w_coeff_cached
              );
          }
        }
  }
}

TEST_P(WCoeffTest, Orthonormality)
// tests orthonormality of (l1, m1) x (l2, m2)
{
  auto [x1, x2] = GetParam();
  // get L, kappa for x1, x2
  std::vector<L_kappa> lk1 = get_coupling_labels(x1);
  std::vector<L_kappa> lk2 = get_coupling_labels(x2);
  // couple x1 to x2
  MultiplicityTagged<u3::SU3>::vector coupling = u3::KroneckerProduct(x1, x2);
  // loop over kappa1, kappa2, L1, L2 labels
  for (L_kappa& tLK1 : lk1)
    for (L_kappa& tLK2 : lk2)
    {
      int kappa1 = tLK1.kappa;
      int kappa2 = tLK2.kappa;
      int L1 = tLK1.L;
      int L2 = tLK2.L;
      // L obeys triangle inequality
      for (int L = std::abs(L1 - L2); L <= (L1 + L2); L++)
      {
        // orthonormality check for summing over SU(3) > SO(3) labels
        double w_sum = compute_W_sum_SU3(
            x1, x2, kappa1, kappa1, L1, L1, kappa2, kappa2, L2, L2, L
          );
        ASSERT_TRUE(std::isfinite(w_sum)) << fmt::format(
            "non-finite sum W({} {} {}; {} {} {}; {}) -> {}",
            x1.Str(),
            kappa1,
            L1,
            x2.Str(),
            kappa2,
            L2,
            L,
            w_sum
          ) << std::endl;
        EXPECT_FLOAT_EQ(w_sum, 1) << fmt::format(
            "Orthogonality failed W({} {} {}; {} {} {}; {}) -> {}",
            x1.Str(),
            kappa1,
            L1,
            x2.Str(),
            kappa2,
            L2,
            L,
            w_sum
          ) << std::endl;
      }
    }
  // loop over x = (l1, m1) x (l2, m2)
  for (int u = 0; u < coupling.size(); u++)
  {
    u3::SU3 x = coupling[u].irrep;
    // outer multiplicity
    int rho_max = u3::OuterMultiplicity(x1, x2, x);
    std::vector<L_kappa> lk = get_coupling_labels(x);
    for (int rho = 1; rho <= rho_max; rho++)
    {
      for (L_kappa& tLK : lk)
      {
        int kappa = tLK.kappa;
        int L = tLK.L;
        // orthonormality check for summing over alpha = Li, kappa_i
        double w_sum_2 = compute_W_sum_alpha(x1, x2, x, x, kappa, L, rho, rho);
        ASSERT_TRUE(std::isfinite(w_sum_2)) << fmt::format(
            "non-finite sum W({}; {}| {} {} {}){} -> {}",
            x1.Str(),
            x2.Str(),
            x.Str(),
            kappa,
            L,
            rho,
            w_sum_2
          ) << std::endl;
        EXPECT_FLOAT_EQ(w_sum_2, 1) << fmt::format(
            "Orthogonality failed W({}; {}| {} {} {}){} -> {}",
            x1.Str(),
            x2.Str(),
            x.Str(),
            kappa,
            L,
            rho,
            w_sum_2
          ) << std::endl;
        ;
      }
    }
  }
}

/*
int main(int argc, char **argv)
{
    // initialization
    // u3::U3CoefInit();

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/
