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
#include "utilities/utilities.h"


// Test fixture for W coefficients
class WCoeffTest
    : public testing::TestWithParam<std::pair<u3::SU3, u3::SU3>>
{
 protected:
  void SetUp() override
  {
    // initialization
    u3::U3CoefInit();
  }

  void TearDown() override {}
};

class UCoeffTest
    : public testing::TestWithParam<std::tuple<u3::SU3, u3::SU3, u3::SU3>>
{
 protected:
  void SetUp() override
  {
    // initialization
    u3::U3CoefInit();
  }

  void TearDown() override {}
};


class ZCoeffTest
    : public testing::TestWithParam<std::tuple<u3::SU3, u3::SU3, u3::SU3>>
{
 protected:
  void SetUp() override
  {
    // initialization
    u3::U3CoefInit();
  }

  void TearDown() override {}
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

INSTANTIATE_TEST_SUITE_P(
    SU3Test,
    UCoeffTest,
    testing::Values(
      std::tuple{u3::SU3{2, 0}, u3::SU3{2, 0}, u3::SU3{2, 0}},
      std::tuple{u3::SU3{1, 2}, u3::SU3{3, 4}, u3::SU3{5, 6}},
      std::tuple{u3::SU3{4, 3}, u3::SU3{2, 4}, u3::SU3{0, 2}}
      )
);

INSTANTIATE_TEST_SUITE_P(
    SU3Test,
    ZCoeffTest,
    testing::Values(
      std::tuple{u3::SU3{2, 0}, u3::SU3{2, 0}, u3::SU3{2, 0}},
      std::tuple{u3::SU3{1, 2}, u3::SU3{3, 4}, u3::SU3{5, 6}},
      std::tuple{u3::SU3{4, 3}, u3::SU3{2, 4}, u3::SU3{0, 2}}
      )
);



/************************************************/
/*************** Helper functions ***************/
/************************************************/

// struct to hold x, L, and kappa
struct L_kappa
{
  // constructor
  L_kappa(const u3::SU3 &x, int L_, int kappa_)
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


std::vector<L_kappa> get_coupling_labels(u3::SU3 &x)
// get vector of allowed L, kappa values for coupling
{
  // branching
  MultiplicityTagged<int>::vector branch = u3::BranchingSO3(x);
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
    u3::SU3 &x1,
    u3::SU3 &x2,
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
      for (L_kappa &tLK : lk)
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
    u3::SU3 &x1,
    u3::SU3 &x2,
    u3::SU3 &x,
    u3::SU3 &xprime,
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
  for (L_kappa &tLK1 : lk1)
  {
    for (L_kappa &tLK2 : lk2)
    {
      sum +=
          u3::W(x1, tLK1.kappa, tLK1.L, x2, tLK2.kappa, tLK2.L, x, kappa, L, rho)
          * u3::W(
              x1, tLK1.kappa, tLK1.L, x2, tLK2.kappa, tLK2.L, xprime, kappa, L, rhoprime
            );
    }
  }
  return sum;
}

double compute_U_sum_x12(u3::SU3 &x1, u3::SU3 &x2, u3::SU3 &x, u3::SU3 &x3, u3::SU3 &x23, 
                    u3::SU3 &x23prime, int r23, int r23prime, int r1_23, int r1_23prime)
// compute orthogonality sum over x12, rho12, rho12_3
{
    double sum = 0;
    // couple x1 to x2 -> x12
    MultiplicityTagged<u3::SU3>::vector x12couple = u3::KroneckerProduct(x1,x2);
    // sum over x12
    for (int u = 0; u < x12couple.size(); u++)
    {
        // x12
        u3::SU3 x12 = x12couple[u].irrep;
        // get rho12, rho12_3 values
        u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12,x23);
        int r12_max = std::get<0>(UMult);
        int r12_3_max = std::get<1>(UMult);
        // loop over r12, r12_3
        for (int r12 = 1; r12 <= r12_max; r12++)
        {
            for (int r12_3 = 1; r12_3 <= r12_3_max; r12_3++)
            {
                // compute sum
                sum += u3::U(x1, x2, x, x3, x12, r12, r12_3, x23, r23, r1_23)
                        *u3::U(x1, x2, x, x3, x12, r12, r12_3, x23prime, r23prime, r1_23prime);
            }
        }

    }
    return sum;
}


double compute_U_sum_x23(u3::SU3 &x1, u3::SU3 &x2, u3::SU3 &x, u3::SU3 &x3, u3::SU3 &x12, 
                    u3::SU3 &x12prime, int r12, int r12prime, int r12_3, int r12_3prime)
// compute orthogonality sum over x23, rho23, rho1_23
{
    double sum = 0;
    // couple x2 to x3 -> x23
    MultiplicityTagged<u3::SU3>::vector x23couple = u3::KroneckerProduct(x2,x3);
    // sum over x12
    for (int u = 0; u < x23couple.size(); u++)
    {
        // x12
        u3::SU3 x23 = x23couple[u].irrep;
        // get rho23, rho1_23 values
        u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12,x23);
        int r23_max = std::get<2>(UMult);
        int r1_23_max = std::get<3>(UMult);
        // loop over r23, r1_23
        for (int r23 = 1; r23 <= r23_max; r23++)
        {
            for (int r1_23 = 1; r1_23 <= r1_23_max; r1_23++)
            {
                // compute sum
                sum += u3::U(x1, x2, x, x3, x12, r12, r12_3, x23, r23, r1_23)
                        *u3::U(x1, x2, x, x3, x12prime, r12prime, r12_3prime, x23, r23, r1_23);
            }
        }

    }
    return sum;
}


double compute_Z_sum_x12(u3::SU3 &x1, u3::SU3 &x2, u3::SU3 &x, u3::SU3 &x3, u3::SU3 &x13, 
                    u3::SU3 &x13prime, int r13, int r13prime, int r13_2, int r13_2prime)
// compute orthogonality sum over x12, rho12, rho12_3
{
    double sum = 0;
    // couple x1 to x2 -> x12
    MultiplicityTagged<u3::SU3>::vector x12couple = u3::KroneckerProduct(x1,x2);
    // sum over x12
    for (int u = 0; u < x12couple.size(); u++)
    {
        // x12
        u3::SU3 x12 = x12couple[u].irrep;
        // get rho12, rho12_3 values
        u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12,x13);
        int r12_max = std::get<0>(UMult);
        int r12_3_max = std::get<1>(UMult);
        // loop over r12, r12_3
        for (int r12 = 1; r12 <= r12_max; r12++)
        {
            for (int r12_3 = 1; r12_3 <= r12_3_max; r12_3++)
            {
                // compute sum
                sum += u3::Z(x1, x2, x, x3, x12, r12, r12_3, x13, r13, r13_2)
                        *u3::Z(x1, x2, x, x3, x12, r12, r12_3, x13prime, r13prime, r13_2prime);
            }
        }

    }
    return sum;
}


double compute_Z_sum_x13(u3::SU3 &x1, u3::SU3 &x2, u3::SU3 &x, u3::SU3 &x3, u3::SU3 &x12, 
                    u3::SU3 &x12prime, int r12, int r12prime, int r12_3, int r12_3prime)
// compute orthogonality sum over x13, rho13, rho13_2
{
    double sum = 0;
    // couple x1 to x3 -> x13
    MultiplicityTagged<u3::SU3>::vector x13couple = u3::KroneckerProduct(x1,x3);
    // sum over x12
    for (int u = 0; u < x13couple.size(); u++)
    {
        // x12
        u3::SU3 x13 = x13couple[u].irrep;
        // get rho13, rho1_13 values
        u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12,x13);
        int r13_max = std::get<2>(UMult);
        int r13_2_max = std::get<3>(UMult);
        // loop over r13, r13_2
        for (int r13 = 1; r13 <= r13_max; r13++)
        {
            for (int r13_2 = 1; r13_2 <= r13_2_max; r13_2++)
            {
                // compute sum
                sum += u3::Z(x1, x2, x, x3, x12, r12, r12_3, x13, r13, r13_2)
                        *u3::Z(x1, x2, x, x3, x12prime, r12prime, r12_3prime, x13, r13, r13_2);
            }
        }

    }
    return sum;
}

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
    for (L_kappa &tupleLK1 : lk1)
      for (L_kappa &tupleLK2 : lk2)
        for (L_kappa &tupleLK3 : lk3)
        {
          int L1 = tupleLK1.L;
          int k1 = tupleLK1.kappa;
          int L2 = tupleLK2.L;
          int k2 = tupleLK2.kappa;
          int L3 = tupleLK3.L;
          int k3 = tupleLK3.kappa;

          // check triangle inequality
          if (L3 >= std::abs(L1 - L2) || L3 <= (L1 + L2))
          {
            // calculate coefficients
            double w_coeff = u3::W(x1, k1, L1, x2, k2, L2, x3, k3, L3, rho);
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
  for (L_kappa &tLK1 : lk1)
    for (L_kappa &tLK2 : lk2)
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
      for (L_kappa &tLK : lk)
      {
        int kappa = tLK.kappa;
        int L = tLK.L;
        // orthonormality check for summing over alpha = Li, kappa_i
        double w_sum_2 = compute_W_sum_alpha(x1, x2, x, x, kappa, L, rho, rho);
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

TEST_P(UCoeffTest, Orthonormality)
// tests orthonormality of (l1, m1) x (l2, m2) x (l3, m3)
{

    auto [x1, x2, x3] = GetParam();

    // couple x1 to x2 -> vector of x12
    MultiplicityTagged<u3::SU3>::vector x12couple = u3::KroneckerProduct(x1,x2);
    // couple x2 to x3 -> vector of x23
    MultiplicityTagged<u3::SU3>::vector x23couple = u3::KroneckerProduct(x2,x3);

    // Loop over x23, r23, r1_23, (lambda, mu) in coupling of (l1, m1) x (l23, m23)
    for (int u = 0; u <x23couple.size(); u++)
    {
        // extract x23
        u3::SU3 x23 = x23couple[u].irrep;
        // couple x1 with x23 to get x
        MultiplicityTagged<u3::SU3>::vector xvec = u3::KroneckerProduct(x1,x23);
        // loop over xvec
        for (int i = 0; i < xvec.size(); i++)
        {
            // extract x
            u3::SU3 x = xvec[i].irrep;
            // calculate max r23, r1_23, and loop over them
            u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12couple[0].irrep,x23);
            int r23_max = std::get<2>(UMult);
            int r1_23_max = std::get<3>(UMult);
            for (int r23 = 1; r23 <= r23_max; r23++)
            {
                for (int r1_23 = 1; r1_23 <= r1_23_max; r1_23++)
                {
                    // calculate orthogonality sum over x12, r12, r12_3
                    double sum12 = compute_U_sum_x12(x1,x2,x,x3,x23,x23,r23,r23,r1_23,r1_23);

                    EXPECT_FLOAT_EQ(sum12, 1)<<fmt::format("orthogonality failed U[{}{}{}{}; (l12,m12)r12 r12_3 {} {} {}] -> {}",
                                    x1.Str(),x2.Str(),x.Str(),x3.Str(),x23.Str(),r23,r1_23,sum12)<<std::endl;
                }
            }
        }
    }


    // Loop over x12, r12, r12_3, (lambda, mu) in coupling of (l12, m12) x (l3, m3)
    for (int u = 0; u <x12couple.size(); u++)
    {
        // extract x12
        u3::SU3 x12 = x12couple[u].irrep;
        // couple x12 with x3 to get x
        MultiplicityTagged<u3::SU3>::vector xvec = u3::KroneckerProduct(x12,x3);
        // loop over xvec
        for (int i = 0; i < xvec.size(); i++)
        {
            // extract x
            u3::SU3 x = xvec[i].irrep;
            // calculate max r12, r12_3, and loop over them
            u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12,x23couple[0].irrep);
            int r12_max = std::get<0>(UMult);
            int r12_3_max = std::get<1>(UMult);
            for (int r12 = 1; r12 <= r12_max; r12++)
            {
                for (int r12_3 = 1; r12_3 <= r12_3_max; r12_3++)
                {
                    // calculate orthogonality sum over x23, r23, r1_23
                    double sum23 = compute_U_sum_x23(x1,x2,x,x3,x12,x12,r12,r12,r12_3,r12_3);

                    EXPECT_FLOAT_EQ(sum23, 1)<<fmt::format("orthogonality failed U[{}{}{}{}; {} {} {} (l23,m23)r23 r1_23] -> {}",
                                    x1.Str(),x2.Str(),x.Str(),x3.Str(),x12.Str(),r12,r12_3,sum23)<<std::endl;
                }
            }
        }
    }


}




TEST_P(ZCoeffTest, Orthonormality)
// tests orthonormality of (l1, m1) x (l2, m2) x (l3, m3)
{
    auto [x1, x2, x3] = GetParam();

    // couple x1 to x2 -> vector of x12
    MultiplicityTagged<u3::SU3>::vector x12couple = u3::KroneckerProduct(x1,x2);
    // couple x1 to x3 -> vector of x13
    MultiplicityTagged<u3::SU3>::vector x13couple = u3::KroneckerProduct(x1,x3);

    // Loop over x13, r13, r13_2, (lambda, mu) in coupling of (l13, m13) x (l2, m2)
    for (int u = 0; u <x13couple.size(); u++)
    {
        // extract x13
        u3::SU3 x13 = x13couple[u].irrep;
        // couple x13 with x2 to get x
        MultiplicityTagged<u3::SU3>::vector xvec = u3::KroneckerProduct(x13,x2);
        // loop over xvec
        for (int i = 0; i < xvec.size(); i++)
        {
            // extract x
            u3::SU3 x = xvec[i].irrep;
            // calculate max r13, r13_2, and loop over them
            u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12couple[0].irrep,x13);
            int r13_max = std::get<2>(UMult);
            int r13_2_max = std::get<3>(UMult);
            for (int r13 = 1; r13 <= r13_max; r13++)
            {
                for (int r13_2 = 1; r13_2 <= r13_2_max; r13_2++)
                {
                    // calculate orthogonality sum over x12, r12, r12_3
                    double sum12 = compute_Z_sum_x12(x1,x2,x,x3,x13,x13,r13,r13,r13_2,r13_2);

                    EXPECT_FLOAT_EQ(sum12, 1)<<fmt::format("orthogonality failed Z[{}{}{}{}; (l12,m12)r12 r12_3 {} {} {}] -> {}",
                                    x1.Str(),x2.Str(),x.Str(),x3.Str(),x13.Str(),r13,r13_2,sum12)<<std::endl;
                }
            }
        }
    }


    // Loop over x12, r12, r12_3, (lambda, mu) in coupling of (l13, m13) x (l2, m2)
    for (int u = 0; u <x12couple.size(); u++)
    {
        // extract x12
        u3::SU3 x12 = x12couple[u].irrep;
        // couple x12 with x3 to get x
        MultiplicityTagged<u3::SU3>::vector xvec = u3::KroneckerProduct(x12,x3);
        // loop over xvec
        for (int i = 0; i < xvec.size(); i++)
        {
            // extract x
            u3::SU3 x = xvec[i].irrep;
            // calculate max r12, r12_3, and loop over them
            u3::UMultiplicityTuple UMult = u3::UMultiplicity(x1,x2,x,x3,x12,x13couple[0].irrep);
            int r12_max = std::get<0>(UMult);
            int r12_3_max = std::get<1>(UMult);
            for (int r12 = 1; r12 <= r12_max; r12++)
            {
                for (int r12_3 = 1; r12_3 <= r12_3_max; r12_3++)
                {
                    // calculate orthogonality sum over x13, r13, r13_2
                    double sum13 = compute_Z_sum_x13(x1,x2,x,x3,x12,x12,r12,r12,r12_3,r12_3);

                    EXPECT_FLOAT_EQ(sum13, 1)<<fmt::format("orthogonality failed Z[{}{}{}{}; {} {} {} (l13,m13)r13 r13_2] -> {}",
                                    x1.Str(),x2.Str(),x.Str(),x3.Str(),x12.Str(),r12,r12_3,sum13)<<std::endl;
                }
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