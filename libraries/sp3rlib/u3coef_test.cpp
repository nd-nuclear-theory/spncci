/****************************************************************
  u3coef_test.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/10/16 (aem,mac): Created based on prototype u3.py and
  T. Dytrych CSU3Master.
****************************************************************/
#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"
// #include "utilities/utilities.h"
#include <map>

#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"

void basic_test()
{
  u3::SU3 x1(2u, 0u);
  u3::SU3 x2(2u, 0u);
  u3::SU3 x12(4u, 0u);
  u3::SU3 x3(2u, 0u);
  u3::SU3 x4(1u, 0u);
  u3::SU3 x23(4u, 0u);
  u3::SU3 x34(3u, 0u);
  u3::SU3 x13(4u, 0u);
  u3::SU3 x24(3u, 0u);
  u3::SU3 x(6u, 0u);
  u3::SU3 xx(7u, 0u);


  std::cout << U(x1, x2, x, x3, x12, 1u, 1u, x23, 1u, 1u) << std::endl;
  std::cout << U(x1, x2, x, x3, x12, 1u, 1u, x23, 1u, 1u) << std::endl;
  std::cout << Z(x1, x2, x, x3, x12, 1u, 1u, x23, 1u, 1u) << std::endl;

  std::cout << W(x1, 1u, 2u, x2, 1u, 2u, x12, 1u, 4u, 1u) << std::endl;
  std::cout << W(x1, 1u, 2u, x2, 1u, 2u, x12, 1u, 2u, 1u) << std::endl;
  std::cout << W(x1, 1u, 2u, x2, 1u, 2u, x12, 1u, 0u, 1u) << std::endl;

}

void block_test()
{
  u3::SU3 x1(2u, 0u);
  u3::SU3 x2(2u, 0u);
  u3::SU3 x12(4u, 0u);
  u3::SU3 x3(2u, 0u);
  u3::SU3 x23(4u, 0u);
  u3::SU3 x(6u, 0u);

  // block access
  std::cout << "U block test" << std::endl;
  u3::UCoefLabels labels(x1, x2, x, x3, x12, x23);
  u3::RecouplingCoefBlock block(labels, u3::RecouplingMode::kU);
  int r12_max, r12_3_max, r23_max, r1_23_max;
  std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();
  std::cout << "multiplicities " << r12_max << " " << r12_3_max << " "
            << r23_max << " " << r1_23_max << std::endl;
  std::cout << block.GetCoef(1, 1, 1, 1) << std::endl;
  std::cout << std::endl;
}

void iteration_test()
{
  // label set
  std::vector<u3::UCoefLabels> label_set;
  u3::SU3 x1(10u, 0u);
  u3::SU3 x2(8u, 5u);
  u3::SU3 x3(5u, 8u);
  MultiplicityTagged<u3::SU3>::vector x12_vector = KroneckerProduct(x1, x2);
  MultiplicityTagged<u3::SU3>::vector x23_vector = KroneckerProduct(x2, x3);
  for (std::size_t i = 0; i < x12_vector.size(); i++)
  {
    u3::SU3 x12 = x12_vector[i].irrep;
    MultiplicityTagged<u3::SU3>::vector x_vector = KroneckerProduct(x12, x3);
    for (std::size_t j = 0; j < x23_vector.size(); j++)
    {
      u3::SU3 x23 = x23_vector[j].irrep;
      for (std::size_t k = 0; k < x_vector.size(); k++)
      {
        u3::SU3 x = x_vector[k].irrep;
        if (u3::OuterMultiplicity(x1, x23, x) > 0)
          label_set.push_back(u3::UCoefLabels(x1, x2, x, x3, x12, x23));
      }
    }
  }

  int countHash = 0;
  std::map<std::size_t, int> uniqueHash;
  for (std::size_t a = 0; a < label_set.size(); a++)
  {
    std::size_t newHash = hash_value(label_set[a]);
    std::cout << label_set[a].Str() << " " << newHash << std::endl;
    countHash++;
    uniqueHash[newHash]++;
  }

  int collisionHash = 0;
  for (std::map<std::size_t, int>::const_iterator it = uniqueHash.begin();
       it != uniqueHash.end();
       ++it)
  {
    //      if (it->second > 1)
    std::cout << it->first << "\t" << it->second << std::endl;
    //        collisionHash++;
  }

  std::cout << "Total number of hashes made: " << std::to_string(countHash)
            << std::endl;
  // std::cout << "Number of hash collision: " << std::to_string(collisionHash) << std::endl;
}


void caching_test()
// Test use of caching wrapper for U coefficients
{
  // generate label set for testing
  u3::SU3 x1(4u, 0u);
  u3::SU3 x2(3u, 2u);
  u3::SU3 x3(2u, 3u);
  // u3::SU3 x1(0u,0u);
  // u3::SU3 x2(0u,0u);
  // u3::SU3 x3(0u,0u);
  std::vector<u3::UCoefLabels> label_set;
  MultiplicityTagged<u3::SU3>::vector x12_values = KroneckerProduct(x1, x2);
  MultiplicityTagged<u3::SU3>::vector x23_values = KroneckerProduct(x2, x3);
  for (auto it12 = x12_values.begin(); it12 != x12_values.end(); ++it12)
    for (auto it23 = x23_values.begin(); it23 != x23_values.end(); ++it23)
    {
      u3::SU3 x12 = it12->irrep;
      u3::SU3 x23 = it23->irrep;
      MultiplicityTagged<u3::SU3>::vector x_values = KroneckerProduct(x12, x3);
      for (auto it = x_values.begin(); it != x_values.end(); ++it)
      {
        u3::SU3 x = it->irrep;
        if (u3::OuterMultiplicity(x1, x23, x) > 0)
        {
          u3::UCoefLabels labels(x1, x2, x, x3, x12, x23);

          label_set.push_back(labels);
        }
      }
    }

  // cache coefficients
  std::cout << "Caching coefficients" << std::endl;
  u3::UCoefCache u_coef_cache;
  for (const u3::UCoefLabels& labels : label_set)
  {
    // if (u_coef_cache.count(labels)>0)
    //   std::cout << "  duplicate " << labels.Str() << std::endl;
    if ((u_coef_cache.size() % 100) == 0)
      std::cout << "  cache size " << u_coef_cache.size() << "..." << std::endl;
    u_coef_cache[labels] =
        u3::RecouplingCoefBlock(labels, u3::RecouplingMode::kU);
  }
  std::cout << "  cached " << u_coef_cache.size() << std::endl;

  // retrieve labels and compare with on-the-fly values
  std::cout << "Checking cached values" << std::endl;
  for (const u3::UCoefLabels& labels : label_set)
  {
    // retrieve labels
    u3::SU3 x1, x2, x, x3, x12, x23;
    std::tie(x1, x2, x, x3, x12, x23) = labels.Key();

    // retrieve coefficient block
    const auto& block = u_coef_cache[labels];

    // retrieve multiplicities
    unsigned int r12_max, r12_3_max, r23_max, r1_23_max;
    std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();

    // loop over multiplicity indices
    // std::cout << "   " << labels.Str() << std::endl;
    for (unsigned int r12 = 1; r12 <= r12_max; ++r12)
      for (unsigned int r12_3 = 1; r12_3 <= r12_3_max; ++r12_3)
        for (unsigned int r23 = 1; r23 <= r23_max; ++r23)
          for (unsigned int r1_23 = 1; r1_23 <= r1_23_max; ++r1_23)
          {
            double coef_direct =
                u3::U(x1, x2, x, x3, x12, r12, r12_3, x23, r23, r1_23);
            double coef_cached = u3::UCached(
                u_coef_cache, x1, x2, x, x3, x12, r12, r12_3, x23, r23, r1_23
              );
            bool compare_ok = (coef_direct == coef_cached);
            if (!compare_ok)
              std::cout << " " << coef_direct << " " << coef_cached << " "
                        << compare_ok << std::endl;
          }
  }
  std::cout << "Done." << std::endl;
}

void OrthogonalitySum1(const u3::SU3& x1, const u3::SU3& x2)
{
  MultiplicityTagged<u3::SU3>::vector product = KroneckerProduct(x1, x2);
  MultiplicityTagged<unsigned int>::vector branch1 = BranchingSO3(x1);
  MultiplicityTagged<unsigned int>::vector branch2 = BranchingSO3(x2);
  // std::cout<<fmt::format("{} {}",x1.Str(),x2.Str())<<std::endl;
  for (std::size_t l1 = 0; l1 < branch1.size(); ++l1)
  {
    unsigned int L1 = branch1[l1].irrep;
    unsigned int kappa1_max = branch1[l1].tag;
    for (std::size_t l2 = 0; l2 < branch2.size(); ++l2)
    {
      unsigned int L2 = branch2[l2].irrep;
      unsigned int kappa2_max = branch2[l2].tag;
      for (unsigned int kappa1 = 1; kappa1 <= kappa1_max; ++kappa1)
        for (unsigned int kappa2 = 1; kappa2 <= kappa2_max; ++kappa2)
        {
          for (unsigned int L =
                   static_cast<unsigned int>(abs(static_cast<int>(L1 - L2)));
               L <= (L1 + L2);
               ++L)
          {
            // sum over x, rho, kappa
            double coef = 0;
            for (std::size_t i = 0; i < product.size(); ++i)
            {
              u3::SU3 x(product[i].irrep);
              // std::cout<<x.Str()<<std::endl;
              unsigned int rho_max = product[i].tag;
              unsigned int kappa_max = u3::BranchingMultiplicitySO3(x, L);
              for (unsigned int rho = 1; rho <= rho_max; ++rho)
                for (unsigned int kappa = 1; kappa <= kappa_max; ++kappa)
                  coef +=
                      u3::W(x1, kappa1, L1, x2, kappa2, L2, x, kappa, L, rho)
                      * u3::W(x1, kappa1, L1, x2, kappa2, L2, x, kappa, L, rho);
            }
            if (fabs(coef - 1) > 10e-10)
              std::cout << fmt::format(
                  "test1: W({} {} {}; {} {} {}; {})  {}",
                  x1.Str(),
                  kappa1,
                  L1,
                  x2.Str(),
                  kappa2,
                  L2,
                  L,
                  coef
                ) << std::endl;
          }
        }
    }
  }
}


void OrthogonalitySum2(
    const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, unsigned int rho
  )
{
  MultiplicityTagged<unsigned int>::vector branch1 = BranchingSO3(x1);
  MultiplicityTagged<unsigned int>::vector branch2 = BranchingSO3(x2);
  MultiplicityTagged<unsigned int>::vector branch = BranchingSO3(x);
  for (std::size_t l = 0; l < branch.size(); ++l)
  {
    unsigned int L = branch[l].irrep;
    unsigned int kappa_max = branch[l].tag;
    for (unsigned int kappa = 1; kappa <= kappa_max; ++kappa)
    {
      double coef = 0;
      // sum over kappa1, L1, kappa2, L2
      for (std::size_t l1 = 0; l1 < branch1.size(); ++l1)
      {
        unsigned int L1 = branch1[l1].irrep;
        unsigned int kappa1_max = branch1[l1].tag;
        for (std::size_t l2 = 0; l2 < branch2.size(); ++l2)
        {
          unsigned int L2 = branch2[l2].irrep;
          unsigned int kappa2_max = branch2[l2].tag;
          for (unsigned int kappa1 = 1; kappa1 <= kappa1_max; ++kappa1)
            for (unsigned int kappa2 = 1; kappa2 <= kappa2_max; ++kappa2)
            {
              coef += u3::W(x1, kappa1, L1, x2, kappa2, L2, x, kappa, L, rho)
                      * u3::W(x1, kappa1, L1, x2, kappa2, L2, x, kappa, L, rho);
            }
        }
      }
      if (fabs(coef - 1) > 10e-10)
        std::cout << fmt::format(
            "test2: W({}; {}| {} {} {}){}  {}",
            x1.Str(),
            x2.Str(),
            x.Str(),
            kappa,
            L,
            rho,
            coef
          ) << std::endl;
    }
  }
}


void TestOrthogonalityW(
    unsigned int lm_min, unsigned int lm_max, unsigned int mu_min, unsigned int mu_max
  )
{
  for (unsigned int l1 = lm_min; l1 <= lm_max; l1++)
    for (unsigned int m1 = mu_min; m1 <= mu_max; m1++)
    {
      u3::SU3 x1(l1, m1);
      for (unsigned int l2 = lm_min; l2 <= lm_max; l2++)
        for (unsigned int m2 = mu_min; m2 <= mu_max; m2++)
        {
          u3::SU3 x2(l2, m2);
          MultiplicityTagged<u3::SU3>::vector product = KroneckerProduct(x1, x2);
          // std::cout<<"Testing Orthogonality sum over kappa1,L1, kappa2,L2"<<std::endl;
          OrthogonalitySum1(x1, x2);
          for (std::size_t i = 0; i < product.size(); ++i)
          {
            u3::SU3 x(product[i].irrep);
            unsigned int rho_max = product[i].tag;
            // std::cout<<"Testing Orthogonality sum over (lambda,mu)rho kappa"<<std::endl;
            for (unsigned int rho = 1; rho <= rho_max; ++rho)
              OrthogonalitySum2(x1, x2, x, rho);
          }
        }
    }
}

void TestWSymmetries13(
    const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x3, unsigned int rho
  )
{
  MultiplicityTagged<unsigned int>::vector branch1 = BranchingSO3(x1);
  MultiplicityTagged<unsigned int>::vector branch2 = BranchingSO3(x2);
  MultiplicityTagged<unsigned int>::vector branch3 = BranchingSO3(x3);
  for (std::size_t l1 = 0; l1 < branch1.size(); ++l1)
  {
    unsigned int L1 = branch1[l1].irrep;
    unsigned int kappa1_max = branch1[l1].tag;
    for (unsigned int l2 = 0; l2 < branch2.size(); ++l2)
    {
      unsigned int L2 = branch2[l2].irrep;
      unsigned int kappa2_max = branch2[l2].tag;
      for (unsigned int l3 = 0; l3 < branch3.size(); ++l3)
      {
        unsigned int L3 = branch3[l3].irrep;
        unsigned int kappa3_max = branch3[l3].tag;
        for (unsigned int kappa1 = 1; kappa1 <= kappa1_max; ++kappa1)
          for (unsigned int kappa2 = 1; kappa2 <= kappa2_max; ++kappa2)
            for (unsigned int kappa3 = 1; kappa3 <= kappa3_max; ++kappa3)
            {
              double coef1 =
                  u3::W(x1, kappa1, L1, x2, kappa2, L2, x3, kappa3, L3, rho);
              double coef12 = u3::W(
                  x3, kappa3, L3, Conjugate(x2), kappa2, L2, x1, kappa1, L1, rho
                );
              double coef2 =
                  ParitySign(
                      x1.lambda() + x1.mu() - x3.lambda() - x3.mu()
                      + static_cast<int>(L1 + L2 - L3)
                    )
                  * std::sqrt(
                      u3::dim(x3) * (2 * L1 + 1.) / (u3::dim(x1) * (2 * L3 + 1.))
                    )
                  * u3::W(
                      x3, kappa3, L3, Conjugate(x2), kappa2, L2, x1, kappa1, L1, rho
                    );

              if (fabs(coef1 - coef2) > 10e-13)
              {
                std::cout << fmt::format(
                    "w({} {} {}; {} {} {}| {} {} {}){}  {}  {}  {}",
                    x1.Str(),
                    kappa1,
                    L1,
                    x2.Str(),
                    kappa2,
                    L2,
                    x3.Str(),
                    kappa3,
                    L3,
                    rho,
                    coef1,
                    coef12,
                    coef2
                  ) << std::endl;
              }
            }
      }
    }
  }
}

void TestWSymmetries(unsigned int lm_max)
{
  for (unsigned int l1 = 0; l1 <= lm_max; l1++)
    for (unsigned int m1 = 0; m1 <= lm_max; m1++)
    {
      u3::SU3 x1(l1, m1);
      for (unsigned int l2 = 0; l2 <= lm_max; l2++)
        for (unsigned int m2 = 0; m2 <= lm_max; m2++)
        {
          u3::SU3 x2(l2, m2);
          MultiplicityTagged<u3::SU3>::vector product = KroneckerProduct(x1, x2);
          for (unsigned int i = 0; i < product.size(); ++i)
          {
            u3::SU3 x3(product[i].irrep);
            unsigned int rho_max = product[i].tag;
            for (unsigned int rho = 1; rho <= rho_max; ++rho)
              TestWSymmetries13(x1, x2, x3, rho);
          }
        }
    }
}

void caching_W_test()
// Test use of caching wrapper for U coefficients
{
  // generate label set for testing
  // u3::SU3 x1(7,8);
  // u3::SU3 x2(6,4);
  u3::SU3 x1(3, 2);
  u3::SU3 x2(2, 2);
  MultiplicityTagged<u3::SU3>::vector x3_values = u3::KroneckerProduct(x1, x2);
  std::vector<u3::WCoefLabels> label_set;
  for (auto it = x3_values.begin(); it != x3_values.end(); ++it)
  {
    u3::SU3 x3(it->irrep);
    MultiplicityTagged<unsigned int>::vector L1_kappa1_values =
        u3::BranchingSO3(x1);
    MultiplicityTagged<unsigned int>::vector L2_kappa2_values =
        u3::BranchingSO3(x2);
    MultiplicityTagged<unsigned int>::vector L3_kappa3_values =
        u3::BranchingSO3(x3);

    for (auto it1 = L1_kappa1_values.begin(); it1 != L1_kappa1_values.end(); ++it1)
      for (auto it2 = L2_kappa2_values.begin(); it2 != L2_kappa2_values.end();
           ++it2)
        for (auto it3 = L3_kappa3_values.begin(); it3 != L3_kappa3_values.end();
             ++it3)
        {
          unsigned int L1 = it1->irrep;
          unsigned int kappa1_max = it1->tag;
          unsigned int L2 = it2->irrep;
          unsigned int kappa2_max = it2->tag;
          unsigned int L3 = it3->irrep;
          unsigned int kappa3_max = it3->tag;
          u3::WCoefLabels labels(x1, L1, x2, L2, x3, L3);
          label_set.push_back(labels);
        }
  }
  // cache coefficients
  std::cout << "Caching coefficients" << std::endl;
  u3::WCoefCache w_coef_cache;
  for (const u3::WCoefLabels& labels : label_set)
  {
    // if (u_coef_cache.count(labels)>0)
    //   std::cout << "  duplicate " << labels.Str() << std::endl;
    //  if ((w_coef_cache.size()%1000)==0)
    //    std::cout << "  cache size " << w_coef_cache.size() << "..." << std::endl;
    w_coef_cache[labels] = u3::WCoefBlock(labels);
  }
  std::cout << "  cached " << w_coef_cache.size() << std::endl;

  // retrieve labels and compare with on-the-fly values
  std::cout << "Checking cached values" << std::endl;
  for (const u3::WCoefLabels& labels : label_set)
  {
    // retrieve labels
    u3::SU3 x1, x2, x3;
    unsigned int L1, L2, L3;
    std::tie(x1, L1, x2, L2, x3, L3) = labels.Key();

    // retrieve coefficient block
    const u3::WCoefBlock& block = w_coef_cache[labels];

    // retrieve multiplicities
    unsigned int rho_max, kappa1_max, kappa2_max, kappa3_max;
    std::tie(kappa1_max, kappa2_max, kappa3_max, rho_max) = block.Key();
    // loop over multiplicity indices
    for (unsigned int kappa1 = 1; kappa1 <= kappa1_max; ++kappa1)
      for (unsigned int kappa2 = 1; kappa2 <= kappa2_max; ++kappa2)
        for (unsigned int kappa3 = 1; kappa3 <= kappa3_max; ++kappa3)
          for (unsigned int rho = 1; rho <= rho_max; ++rho)
          {
            double coef_direct =
                u3::W(x1, kappa1, L1, x2, kappa2, L2, x3, kappa3, L3, rho);
            double coef_cached = u3::WCached(
                w_coef_cache, x1, kappa1, L1, x2, kappa2, L2, x3, kappa3, L3, rho
              );
            bool compare_ok = (coef_direct == coef_cached);
            if (!compare_ok)
              std::cout << " " << coef_direct << " " << coef_cached << " "
                        << compare_ok << std::endl;
          }
  }
  std::cout << "Done." << std::endl;
}


double UTest(
    const u3::SU3& x1,
    const u3::SU3& x2,
    const u3::SU3& x,
    const u3::SU3& x3,
    const u3::SU3& x12,
    unsigned int r12,
    unsigned int r12_3,
    const u3::SU3& x23,
    unsigned int r23,
    unsigned int r1_23
  )
{
  double sum = 0;
  MultiplicityTagged<unsigned int>::vector L1_set = u3::BranchingSO3(x1);
  MultiplicityTagged<unsigned int>::vector L2_set = u3::BranchingSO3(x2);
  MultiplicityTagged<unsigned int>::vector L3_set = u3::BranchingSO3(x3);
  MultiplicityTagged<unsigned int>::vector L12_set = u3::BranchingSO3(x12);
  MultiplicityTagged<unsigned int>::vector L23_set = u3::BranchingSO3(x23);
  MultiplicityTagged<unsigned int>::vector L_set = u3::BranchingSO3(x);
  for (unsigned int i1 = 0; i1 < L1_set.size(); ++i1)
    for (unsigned int i2 = 0; i2 < L2_set.size(); ++i2)
      for (unsigned int i3 = 0; i3 < L3_set.size(); ++i3)
        for (unsigned int i12 = 0; i12 < L12_set.size(); ++i12)
          for (unsigned int i23 = 0; i23 < L23_set.size(); ++i23)
          // for(int i=0; i<L_set.size(); ++i)
          {
            unsigned int L1 = L1_set[i1].irrep;
            unsigned int kappa1_max = L1_set[i1].tag;
            unsigned int L2 = L2_set[i2].irrep;
            unsigned int kappa2_max = L2_set[i2].tag;
            unsigned int L3 = L3_set[i3].irrep;
            unsigned int kappa3_max = L3_set[i3].tag;
            unsigned int L12 = L12_set[i12].irrep;
            unsigned int kappa12_max = L12_set[i12].tag;
            unsigned int L23 = L23_set[i23].irrep;
            unsigned int kappa23_max = L23_set[i23].tag;
            unsigned int L = L_set[0].irrep;
            for (unsigned int kappa1 = 1; kappa1 <= kappa1_max; ++kappa1)
              for (unsigned int kappa2 = 1; kappa2 <= kappa2_max; ++kappa2)
                for (unsigned int kappa3 = 1; kappa3 <= kappa3_max; ++kappa3)
                  for (unsigned int kappa12 = 1; kappa12 <= kappa12_max; ++kappa12)
                    for (unsigned int kappa23 = 1; kappa23 <= kappa23_max;
                         ++kappa23)
                    {
                      sum +=
                          u3::W(
                              x1, kappa1, L1, x2, kappa2, L2, x12, kappa12, L12, r12
                            )
                          * u3::W(x12, kappa12, L12, x3, kappa3, L3, x, 1, L, r12_3)
                          * u3::W(
                              x2, kappa2, L2, x3, kappa3, L3, x23, kappa23, L23, r23
                            )
                          * u3::W(x1, kappa1, L1, x23, kappa23, L23, x, 1, L, r1_23)
                          * am::Unitary6J(L1, L2, L12, L3, L, L23);
                    }
          }
  return sum;
}


void phi_caching_test()
// Test use of caching wrapper for U coefficients
{
  std::vector<u3::PhiCoefLabels> label_set;

  for (unsigned int lambda1 = 1; lambda1 < 3; ++lambda1)
    for (unsigned int mu1 = 1; mu1 < 3; ++mu1)
      for (unsigned int lambda2 = 1; lambda2 < 4; ++lambda2)
        for (unsigned int mu2 = 1; mu2 < 4; ++mu2)
        {
          u3::SU3 x1(lambda1, mu1);
          u3::SU3 x2(lambda2, mu2);
          MultiplicityTagged<u3::SU3>::vector x3_set =
              u3::KroneckerProduct(x1, x2);
          for (auto x3_tagged : x3_set)
          {
            u3::PhiCoefLabels labels(x1, x2, x3_tagged.irrep);
            label_set.push_back(labels);
          }
        }
  // cache coefficients
  std::cout << "Caching coefficients" << std::endl;
  u3::PhiCoefCache phi_coef_cache;
  for (const u3::PhiCoefLabels& labels : label_set)
  {
    if ((phi_coef_cache.size() % 100) == 0)
      std::cout << "  cache size " << phi_coef_cache.size() << "..." << std::endl;
    phi_coef_cache[labels] = u3::PhiCoefBlock(labels);
  }
  std::cout << "  cached " << phi_coef_cache.size() << std::endl;

  // retrieve labels and compare with on-the-fly values
  std::cout << "Checking cached values" << std::endl;
  for (const u3::PhiCoefLabels& labels : label_set)
  {
    // retrieve labels
    u3::SU3 x1, x2, x3;
    std::tie(x1, x2, x3) = labels.Key();
    // retrieve coefficient block
    // const u3::PhiCoefBlock& block = phi_coef_cache[labels];

    unsigned int rho_max = u3::OuterMultiplicity(x1, x2, x3);
    for (unsigned int rho1 = 1; rho1 <= rho_max; ++rho1)
      for (unsigned int rho2 = 1; rho2 <= rho_max; ++rho2)
      {
        double coef_direct = u3::Phi(x1, x2, x3, rho1, rho2);
        double coef_cached =
            u3::PhiCached(phi_coef_cache, x1, x2, x3, rho1, rho2);
        bool compare_ok = (coef_direct == coef_cached);
        if (!compare_ok)
          std::cout << " " << coef_direct << " " << coef_cached << " "
                    << compare_ok << std::endl;
      }
  }
  std::cout << "Done." << std::endl;
}


int main(int argc, char** argv)
{
  // initialize su3lib
  u3::U3CoefInit(39);

  // U coefficient test--comparison with escher formula for U in terms of
  // SU(3)\supset SO(3) coupling coefficients
  u3::SU3 x1(2u, 2u);
  u3::SU3 x2(2u, 0u);
  u3::SU3 x(2u, 1u);
  u3::SU3 x3(2u, 0u);
  u3::SU3 x12(2u, 0u);
  u3::SU3 x23(0u, 2u);
  std::cout << u3::U(x1, x2, x, x3, x12, 1, 1, x23, 1, 1) << std::endl;


  // W coefficient test--comparison with prototype
  x1 = u3::SU3(4u, 2u);
  x2 = u3::SU3(4u, 1u);
  x = u3::SU3(9u, 1u);
  unsigned int kappa1 = 2, kappa2 = 1, kappa = 1;
  unsigned int L1 = 2, L2 = 5, L = 7;
  unsigned int rho = 1;
  std::cout << "W test " << std::endl;
  std::cout << "kappa1_max " << u3::BranchingMultiplicitySO3(x1, L1) << std::endl;
  std::cout << "kappa2_max " << u3::BranchingMultiplicitySO3(x2, L2) << std::endl;
  std::cout << "kappa3_max " << u3::BranchingMultiplicitySO3(x, L) << std::endl;
  std::cout << "coef= " << u3::W(x1, kappa1, L1, x2, kappa2, L2, x, kappa, L, rho)
            << std::endl;
  // basic tests
  basic_test();

  phi_caching_test();
}
