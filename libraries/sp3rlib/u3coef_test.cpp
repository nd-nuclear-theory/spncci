/****************************************************************
  u3coef_test.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.
****************************************************************/
#include "cppformat/format.h"

#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "utilities/utilities.h"
#include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include <map>

void basic_test()
{
  u3::SU3 x1(2,0);
  u3::SU3 x2(2,0);
  u3::SU3 x12(4,0);
  u3::SU3 x3(2,0);
  u3::SU3 x4(1,0);
  u3::SU3 x23(4,0);
  u3::SU3 x34(3,0);
  u3::SU3 x13(4,0);
  u3::SU3 x24(3,0);
  u3::SU3 x(6,0);
  u3::SU3 xx(7,0);


  // U(x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23);
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << Z(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;

  std::cout << W(x1,1,2,x2,1,2,x12,1,4,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,2,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,0,1) << std::endl;

  std::cout << Unitary9LambdaMu(
                                x1,  x2,  x12, 1,
                                x3,  x4,  x34, 1,
                                x13, x24, xx,  1,
                                1,    1,   1
                                )<<std::endl;

}

void block_test()
{

  u3::SU3 x1(2,0);
  u3::SU3 x2(2,0);
  u3::SU3 x12(4,0);
  u3::SU3 x3(2,0);
  u3::SU3 x23(4,0);
  u3::SU3 x(6,0);

  // block access
  std::cout << "U block test" << std::endl;
  u3::UCoefLabels labels(x1, x2, x, x3, x12, x23);
  u3::UCoefBlock block(labels);
  int r12_max, r12_3_max, r23_max, r1_23_max;
  std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();
  std::cout << "multiplicities " << r12_max << " " << r12_3_max << " " << r23_max << " " << r1_23_max << std::endl;
  std::cout << block.GetCoef(1,1,1,1) << std::endl;
  std::cout << std::endl;
  
}

void iteration_test()
{

  // label set
  std::vector<u3::UCoefLabels> label_set;
  u3::SU3 x1(10,0);
  u3::SU3 x2(8,5);
  u3::SU3 x3(5,8);
  MultiplicityTagged<u3::SU3>::vector x12_vector=KroneckerProduct(x1,x2);
  MultiplicityTagged<u3::SU3>::vector x23_vector=KroneckerProduct(x2,x3);
  for(int i=0; i<x12_vector.size(); i++)
    {
      u3::SU3 x12=x12_vector[i].irrep;
      MultiplicityTagged<u3::SU3>::vector x_vector=KroneckerProduct(x12,x3);
      for(int j=0; j<x23_vector.size(); j++)
        {
          u3::SU3 x23=x23_vector[j].irrep;
          for(int k=0; k<x_vector.size(); k++)
            {
              u3::SU3 x=x_vector[k].irrep;
              if(u3::OuterMultiplicity(x1,x23,x)>0)
                label_set.push_back(u3::UCoefLabels(x1,x2,x,x3,x12,x23));
            }
        }
    }

  int countHash = 0;
  std::map<std::size_t,int> uniqueHash;
  for (int a=0; a<label_set.size(); a++)
    {
      std::size_t newHash = hash_value(label_set[a]);
      std::cout << label_set[a].Str() << " " << newHash << std::endl;
      countHash++;
      uniqueHash[newHash]++;

    }

  int collisionHash = 0;
  for (std::map<std::size_t,int>::const_iterator it = uniqueHash.begin();
       it != uniqueHash.end(); ++it)
    {
      //      if (it->second > 1)
      std::cout << it->first << "\t" << it->second << std::endl;
      //        collisionHash++;
    }

  std::cout << "Total number of hashes made: " << std::to_string(countHash) << std::endl;
  //std::cout << "Number of hash collision: " << std::to_string(collisionHash) << std::endl;
}


void caching_test()
// Test use of caching wrapper for U coefficients
{

  // generate label set for testing
  u3::SU3 x1(4,0);
  u3::SU3 x2(3,2);
  u3::SU3 x3(2,3);
  // u3::SU3 x1(0,0);
  // u3::SU3 x2(0,0);
  // u3::SU3 x3(0,0);
  std::vector<u3::UCoefLabels> label_set;
  MultiplicityTagged<u3::SU3>::vector x12_values=KroneckerProduct(x1,x2);
  MultiplicityTagged<u3::SU3>::vector x23_values=KroneckerProduct(x2,x3);
  for (auto it12 = x12_values.begin(); it12 != x12_values.end(); ++it12)
    for (auto it23 = x23_values.begin(); it23 != x23_values.end(); ++it23)
      {
        u3::SU3 x12 = it12->irrep;
        u3::SU3 x23 = it23->irrep;
        MultiplicityTagged<u3::SU3>::vector x_values=KroneckerProduct(x12,x3);
        for(auto it = x_values.begin(); it != x_values.end(); ++it)
          {
            u3::SU3 x = it->irrep;
            if(u3::OuterMultiplicity(x1,x23,x)>0)
              {
                u3::UCoefLabels labels(x1,x2,x,x3,x12,x23);

                label_set.push_back(labels);
              }
          }
      }

  // cache coefficients
  std::cout << "Caching coefficients" << std::endl;
  u3::UCoefCache u_coef_cache;
  for (const u3::UCoefLabels& labels : label_set)
    {
      //if (u_coef_cache.count(labels)>0)
      //  std::cout << "  duplicate " << labels.Str() << std::endl;
      if ((u_coef_cache.size()%100)==0)
        std::cout << "  cache size " << u_coef_cache.size() << "..." << std::endl;
      u_coef_cache[labels] = u3::UCoefBlock(labels);
    }
  std::cout << "  cached " << u_coef_cache.size() << std::endl;

  // retrieve labels and compare with on-the-fly values
  std::cout << "Checking cached values" << std::endl;
  for (const u3::UCoefLabels& labels : label_set)
    {
      // retrieve labels
      u3::SU3 x1, x2, x, x3, x12, x23;
      std::tie(x1,x2,x,x3,x12,x23) = labels.Key();

      // retrieve coefficient block
      const u3::UCoefBlock& block = u_coef_cache[labels];

      // retrieve multiplicities
      int r12_max, r12_3_max, r23_max, r1_23_max;
      std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();

      // loop over multiplicity indices
      // std::cout << "   " << labels.Str() << std::endl;
      for (int r12 = 1; r12 <= r12_max; ++r12)
        for (int r12_3 = 1; r12_3 <= r12_3_max; ++r12_3)
          for (int r23 = 1; r23 <= r23_max; ++r23)
            for (int r1_23 = 1; r1_23 <= r1_23_max; ++r1_23)
              {
                double coef_direct = u3::U(x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23);
                double coef_cached = u3::UCached(u_coef_cache,x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23);
                bool compare_ok = (coef_direct == coef_cached);
                if (!compare_ok)
                  std::cout << " " << coef_direct << " " << coef_cached << " " << compare_ok << std::endl;
              }
    }
  std::cout << "Done." << std::endl;

}

void OrthogonalitySum1(const u3::SU3& x1, const u3::SU3& x2)
{
  MultiplicityTagged<u3::SU3>::vector product=KroneckerProduct(x1,x2);
  MultiplicityTagged<int>::vector branch1=BranchingSO3(x1);
  MultiplicityTagged<int>::vector branch2=BranchingSO3(x2);
  // std::cout<<fmt::format("{} {}",x1.Str(),x2.Str())<<std::endl;
  for(int l1=0; l1<branch1.size(); ++l1)
    {
      int L1=branch1[l1].irrep;
      int kappa1_max=branch1[l1].tag;
      for (int l2=0; l2<branch2.size(); ++l2)
        {
          int L2=branch2[l2].irrep;
          int kappa2_max=branch2[l2].tag;
          for(int kappa1=1; kappa1<=kappa1_max; ++kappa1)
            for(int kappa2=1; kappa2<=kappa2_max; ++kappa2)
              {
                for(int L=abs(L1-L2); L<=(L1+L2); ++L)
                  {
                    // sum over x, rho, kappa
                    double coef=0;
                    for(int i=0; i<product.size(); ++i)
                      {
                        u3::SU3 x(product[i].irrep);
                        // std::cout<<x.Str()<<std::endl;
                        int rho_max=product[i].tag;
                        int kappa_max=u3::BranchingMultiplicitySO3(x,L);
                        for(int rho=1; rho<=rho_max; ++rho)  
                          for(int kappa=1; kappa<=kappa_max; ++kappa)
                              coef+=u3::W(x1,kappa1,L1,x2,kappa2,L2,x,kappa,L,rho)
                                                    *u3::W(x1,kappa1,L1,x2,kappa2,L2,x,kappa,L,rho);
                      }
                    if(fabs(coef-1)>10e-10)
                      std::cout<<fmt::format("test1: W({} {} {}; {} {} {}; {})  {}",
                                             x1.Str(),kappa1,L1,x2.Str(),kappa2,L2,L,coef)
                      <<std::endl;      
                  }
              }
                  
        }
    }
}


void OrthogonalitySum2(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, int rho)
  {
    MultiplicityTagged<int>::vector branch1=BranchingSO3(x1);
    MultiplicityTagged<int>::vector branch2=BranchingSO3(x2);
    MultiplicityTagged<int>::vector branch=BranchingSO3(x);
    for(int l=0; l<branch.size(); ++l)
      {
        int L=branch[l].irrep;
        int kappa_max=branch[l].tag;
        for(int kappa=1; kappa<=kappa_max; ++kappa)
          {
            double coef=0;
            //sum over kappa1, L1, kappa2, L2
            for(int l1=0; l1<branch1.size(); ++l1)
              {
                int L1=branch1[l1].irrep;
                int kappa1_max=branch1[l1].tag;
                for(int l2=0; l2<branch2.size(); ++l2)
                  {
                    int L2=branch2[l2].irrep;
                    int kappa2_max=branch2[l2].tag;
                    for(int kappa1=1; kappa1<=kappa1_max; ++kappa1)
                      for(int kappa2=1; kappa2<=kappa2_max; ++kappa2)
                        {
                          coef+=u3::W(x1,kappa1,L1,x2,kappa2,L2,x,kappa,L,rho)
                                *u3::W(x1,kappa1,L1,x2,kappa2,L2,x,kappa,L,rho);
                        }
                  }
              }
            if(fabs(coef-1)>10e-10)
              std::cout<<fmt::format("test2: W({}; {}| {} {} {}){}  {}",
                                     x1.Str(),x2.Str(),x.Str(), kappa,L,rho,coef)
              <<std::endl;      
          }
      }


  }


void TestOrthogonalityW(int lm_min, int lm_max, int mu_min, int mu_max)
  {

    for(int l1=lm_min; l1<=lm_max; l1++)
      for(int m1=mu_min; m1<=mu_max; m1++)
        {
          u3::SU3 x1(l1,m1);
          for(int l2=lm_min; l2<=lm_max; l2++)
            for(int m2=mu_min; m2<=mu_max; m2++)
              {
                u3::SU3 x2(l2,m2);
                MultiplicityTagged<u3::SU3>::vector product=KroneckerProduct(x1,x2);
                // std::cout<<"Testing Orthogonality sum over kappa1,L1, kappa2,L2"<<std::endl;
                OrthogonalitySum1(x1,x2);
                for(int i=0; i<product.size(); ++i)
                  {
                    u3::SU3 x(product[i].irrep);
                    int rho_max=product[i].tag;
                    // std::cout<<"Testing Orthogonality sum over (lambda,mu)rho kappa"<<std::endl;
                    for(int rho=1; rho<=rho_max; ++rho)
                      OrthogonalitySum2(x1, x2, x, rho);
                  }
              }
        }
  }

void TestWSymmetries13(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x3, int rho)
  {
    MultiplicityTagged<int>::vector branch1=BranchingSO3(x1);
    MultiplicityTagged<int>::vector branch2=BranchingSO3(x2);
    MultiplicityTagged<int>::vector branch3=BranchingSO3(x3);
    for(int l1=0; l1<branch1.size(); ++l1)
      {
        int L1=branch1[l1].irrep;
        int kappa1_max=branch1[l1].tag;
        for(int l2=0; l2<branch2.size(); ++l2)
          {
            int L2=branch2[l2].irrep;
            int kappa2_max=branch2[l2].tag;
            for(int l3=0; l3<branch3.size(); ++l3)
              {
                int L3=branch3[l3].irrep;
                int kappa3_max=branch3[l3].tag;
                for(int kappa1=1; kappa1<=kappa1_max; ++kappa1)
                  for(int kappa2=1; kappa2<=kappa2_max; ++kappa2)
                    for(int kappa3=1; kappa3<=kappa3_max; ++kappa3)
                      {
                        double coef1=u3::W(x1,kappa1,L1,x2,kappa2,L2,x3,kappa3,L3,rho);
                        double coef12=u3::W(x3,kappa3,L3,Conjugate(x2),kappa2,L2,x1,kappa1,L1,rho);
                        double coef2=parity(x1.lambda()+x1.mu()-x3.lambda()-x3.mu()+L1+L2-L3)
                                      *std::sqrt(u3::dim(x3)*(2*L1+1.)/(u3::dim(x1)*(2*L3+1.)))
                                      *u3::W(x3,kappa3,L3,Conjugate(x2),kappa2,L2,x1,kappa1,L1,rho);
                        
                        if(fabs(coef1-coef2)>10e-13)
                          {
                            std::cout<<fmt::format("w({} {} {}; {} {} {}| {} {} {}){}  {}  {}  {}",       
                            x1.Str(),kappa1,L1,x2.Str(),kappa2,L2,x3.Str(),kappa3,L3,rho,coef1,coef12,coef2)<<std::endl;
                          }
                      }
              }
          }
      }
  }

void TestWSymmetries(int lm_max)
{
for(int l1=0; l1<=lm_max; l1++)
    for(int m1=0; m1<=lm_max; m1++)
      {
        u3::SU3 x1(l1,m1);
        for(int l2=0; l2<=lm_max; l2++)
          for(int m2=0; m2<=lm_max; m2++)
            {
              u3::SU3 x2(l2,m2);
              MultiplicityTagged<u3::SU3>::vector product=KroneckerProduct(x1,x2);
              for(int i=0; i<product.size(); ++i)
                {
                  u3::SU3 x3(product[i].irrep);
                  int rho_max=product[i].tag;
                  for(int rho=1; rho<=rho_max; ++rho)
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
  u3::SU3 x1(3,2);
  u3::SU3 x2(2,2);
  MultiplicityTagged<u3::SU3>::vector x3_values=u3::KroneckerProduct(x1,x2);
  std::vector<u3::WCoefLabels> label_set;
  for(auto it=x3_values.begin(); it!=x3_values.end(); ++it)
    {
      u3::SU3 x3(it->irrep);
      MultiplicityTagged<int>::vector L1_kappa1_values=u3::BranchingSO3(x1);
      MultiplicityTagged<int>::vector L2_kappa2_values=u3::BranchingSO3(x2);
      MultiplicityTagged<int>::vector L3_kappa3_values=u3::BranchingSO3(x3);
      
      for (auto it1 = L1_kappa1_values.begin(); it1!= L1_kappa1_values.end(); ++it1)
        for (auto it2 = L2_kappa2_values.begin(); it2!= L2_kappa2_values.end(); ++it2)
          for (auto it3 = L3_kappa3_values.begin(); it3!= L3_kappa3_values.end(); ++it3)
            {
              int L1=it1->irrep;
              int kappa1_max=it1->tag;
              int L2=it2->irrep;
              int kappa2_max=it2->tag;
              int L3=it3->irrep;
              int kappa3_max=it3->tag;
              u3::WCoefLabels labels(x1,L1,x2,L2,x3,L3);
              label_set.push_back(labels);
             }
    }
  // cache coefficients
  std::cout << "Caching coefficients" << std::endl;
  u3::WCoefCache w_coef_cache;
  for (const u3::WCoefLabels& labels : label_set)
    {
      //if (u_coef_cache.count(labels)>0)
      //  std::cout << "  duplicate " << labels.Str() << std::endl;
      // if ((w_coef_cache.size()%1000)==0)
      //   std::cout << "  cache size " << w_coef_cache.size() << "..." << std::endl;
      w_coef_cache[labels] = u3::WCoefBlock(labels);
    }
  std::cout << "  cached " << w_coef_cache.size() << std::endl;

  // retrieve labels and compare with on-the-fly values
  std::cout << "Checking cached values" << std::endl;
  for (const u3::WCoefLabels& labels : label_set)
    {
      // retrieve labels
      u3::SU3 x1, x2, x3;
      int L1,L2,L3;
      std::tie(x1,L1,x2,L2,x3,L3) = labels.Key();

      // retrieve coefficient block
      const u3::WCoefBlock& block = w_coef_cache[labels];

      // retrieve multiplicities
      int rho_max, kappa1_max, kappa2_max, kappa3_max;
      std::tie(kappa1_max, kappa2_max, kappa3_max, rho_max) = block.Key();
      // loop over multiplicity indices
      for(int kappa1=1; kappa1<=kappa1_max; ++kappa1)
        for(int kappa2=1; kappa2<=kappa2_max; ++kappa2)
          for(int kappa3=1; kappa3<=kappa3_max; ++kappa3)
            for(int rho=1; rho<=rho_max; ++rho)
              {
                
                double coef_direct = u3::W(x1,kappa1,L1,x2,kappa2,L2,x3,kappa3,L3,rho);
                double coef_cached = u3::WCached(w_coef_cache,x1,kappa1,L1,x2,kappa2,L2,x3,kappa3,L3,rho);
                bool compare_ok = (coef_direct == coef_cached);
                if (!compare_ok)
                  std::cout << " " << coef_direct << " " << coef_cached << " " << compare_ok << std::endl;
              }
    }
  std::cout << "Done." << std::endl;

}



double 
UTest(
  const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3,
  const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23
  )
{
  double sum=0;
  MultiplicityTagged<int>::vector L1_set=u3::BranchingSO3(x1);
  MultiplicityTagged<int>::vector L2_set=u3::BranchingSO3(x2);
  MultiplicityTagged<int>::vector L3_set=u3::BranchingSO3(x3);
  MultiplicityTagged<int>::vector L12_set=u3::BranchingSO3(x12);
  MultiplicityTagged<int>::vector L23_set=u3::BranchingSO3(x23);
  MultiplicityTagged<int>::vector L_set=u3::BranchingSO3(x);
  for(int i1=0; i1<L1_set.size(); ++i1)
    for(int i2=0; i2<L2_set.size(); ++i2)
      for(int i3=0; i3<L3_set.size(); ++i3)
        for(int i12=0; i12<L12_set.size(); ++i12)
          for(int i23=0; i23<L23_set.size(); ++i23)
            // for(int i=0; i<L_set.size(); ++i)
              {
                int L1=L1_set[i1].irrep;
                int kappa1_max=L1_set[i1].tag;
                int L2=L2_set[i2].irrep;
                int kappa2_max=L2_set[i2].tag;
                int L3=L3_set[i3].irrep;
                int kappa3_max=L3_set[i3].tag;
                int L12=L12_set[i12].irrep;
                int kappa12_max=L12_set[i12].tag;
                int L23=L23_set[i23].irrep;
                int kappa23_max=L23_set[i23].tag;
                int L=L_set[0].irrep;
                for(int kappa1=1; kappa1<=kappa1_max; ++kappa1)
                  for(int kappa2=1; kappa2<=kappa2_max; ++kappa2)
                    for(int kappa3=1; kappa3<=kappa3_max; ++kappa3)
                      for(int kappa12=1; kappa12<=kappa12_max; ++kappa12)
                        for(int kappa23=1; kappa23<=kappa23_max; ++kappa23)
                            {
                              sum
                              +=u3::W(x1,kappa1,L1,x2,kappa2,L2,x12,kappa12,L12,r12)
                              *u3::W(x12, kappa12,L12,x3,kappa3,L3,x,1,L,r12_3)
                              *u3::W(x2,kappa2,L2,x3,kappa3,L3,x23,kappa23,L23,r23)
                              *u3::W(x1,kappa1,L1,x23,kappa23,L23,x,1,L,r1_23)
                              *am::Unitary6J(L1,L2,L12,L3,L,L23);
                            }


              }
  return sum;
}



int main(int argc, char **argv)
{
  // initialize su3lib
  u3::U3CoefInit();

  // // U coefficient test--comparison with escher formula for U in terms of SU(3)\supset SO(3) coupling coefficients
  // u3::SU3 x1(2,2);
  // u3::SU3 x2(2,0);
  // u3::SU3 x(2,1);
  // u3::SU3 x3(2,0);
  // u3::SU3 x12(2,0); 
  // u3::SU3 x23(0,2);
  // std::cout<<u3::U(x1, x2, x,x3,x12, 1,1,x23,1,1)<<std::endl; 
  // std::cout<<UTest(x1, x2, x,x3,x12, 1,1,x23,1,1)<<std::endl;
  
  // W coefficient test--comparison with prototype
  u3::SU3 x1(4,2);
  u3::SU3 x2(4,1);
  u3::SU3 x(9,1);
  int kappa1=2, kappa2=2, kappa=1;
  int L1=5, L2=5, L=10;
  int rho=1;
  std::cout<<"W test "<<std::endl;
  std::cout<<"coef= "<<u3::W(x1,kappa1,L1,x2,kappa2,L2,x,kappa,L,rho)<<std::endl;
  // // basic tests
  // basic_test();

  // // test of block storage
  // block_test();
  
  // // // test of hash key
  // // iteration_test();

  // // test cache storage and retrieval
  // caching_test();  

  //test symmetries of W coefficients 
  int lm_max=4;
  int lm_min=2;
  int mu_min=0;
  int mu_max=2;

  // Test prints out coefficients that fail the test.  
  // Test expected for fail in cases where for x2, lambda=mu.
  // TestWSymmetries(lm_max);
  
  //test orthogonality of W coefficients 
  TestOrthogonalityW(lm_min,lm_max, mu_min,mu_max);

  // caching_W_test();

  // for(int q1=0; q1<=20; ++q1)
  //   for(int q2=0; q2<=20; ++q2)
  //     for(int b=0; b<=20; b++)
  //       {
  //         double coef1=U(u3::SU3(q1+q2+b,0),u3::SU3(0,q2+b),u3::SU3(q1+q2,0),u3::SU3(q2,0),u3::SU3(q1,0),1,1,u3::SU3(0,b),1,1);
  //         double coef2=sqrt(1.*dim(u3::SU3(q1,0))*dim(u3::SU3(b,0))/(dim(u3::SU3(q1+q2,0))*dim(u3::SU3(q2+b,0))));
  //         if(fabs(coef1-coef2)>10e-10)
  //           std::cout<<fmt::format( "{} {} {}   {}  {}",q1,q2,b,coef1,coef2)<<std::endl;

  //       }

  // for(int q=0; q<=4; ++q)
  //   for(int rb=0; rb<=4; ++rb)
  //     for(int rbp=0; rbp<=4; rbp++)
  //       {
  //         if(rb<q)
  //           continue;
  //         MultiplicityTagged<u3::SU3>::vector w0_set=KroneckerProduct(u3::SU3(rbp,0),u3::SU3(0,rb));
  //         MultiplicityTagged<u3::SU3>::vector w0pp_set=KroneckerProduct(u3::SU3(rbp,0),u3::SU3(0,rb-q));
  //         for(int i=0; i<w0_set.size(); ++i) 
  //           for(int j=0; j<w0pp_set.size(); ++j)
  //             {
  //               u3::SU3 x0=w0_set[i].irrep;
  //               u3::SU3 x0pp=w0pp_set[j].irrep;
  //               if(
  //                 (u3::OuterMultiplicity(x0,u3::SU3(q,0),x0pp)>0)
  //                 &&(u3::OuterMultiplicity(u3::SU3(0,rb),u3::SU3(q,0),u3::SU3(0,rb-q))>0)
  //                 )
  //               {
  //                 std::cout<<x0.Str()<<"  "<<x0pp.Str()<<std::endl;

  //                 double coef1=U(u3::SU3(rbp,0),u3::SU3(0,rb),x0pp,u3::SU3(q,0),x0,1,1,u3::SU3(0,rb-q),1,1);
  //                 double coef2=sqrt(1.*dim(u3::SU3(rb-q,0))*dim(x0)/(dim(x0pp)*dim(u3::SU3(rb,0))));
  //                 if(fabs(coef1-coef2)>10e-10)
  //                   std::cout<<fmt::format( "{} {} {}   {}  {}",q,rbp,rb,coef1,coef2)<<std::endl;
  //               }
  //             }
  //       }
  // int Nmax=3;
  // for(int lambda1=2; lambda1<=2; ++lambda1)
  //   for(int mu1=2; mu1<=2; ++mu1)
  //     for(int lambda2=0; lambda2<=0; ++lambda2) //For comparison
  //       for(int mu2=0; mu2<=0; ++mu2)
  //         for(int lambda3=2; lambda3<=2; ++lambda3)
  //           for(int mu3=0; mu3<=0; ++mu3)
  //             for(int lambda4=0; lambda4<=Nmax; ++lambda4) //restrict to case when coef=1
  //               for(int mu4=0; mu4<=0; ++mu4)
  //                 {
  //                   u3::SU3 x1(lambda1,mu1), x2(lambda2,mu2), x3(lambda3,mu3), x4(lambda4,mu4);
  //                   MultiplicityTagged<u3::SU3>::vector x12_set=u3::KroneckerProduct(x1,x2);
  //                   MultiplicityTagged<u3::SU3>::vector x34_set=u3::KroneckerProduct(x3,x4);
  //                   MultiplicityTagged<u3::SU3>::vector x13_set=u3::KroneckerProduct(x1,x3);
  //                   MultiplicityTagged<u3::SU3>::vector x24_set=u3::KroneckerProduct(x2,x4);
  //                   for (auto a12 : x12_set)
  //                     for(auto a34 : x34_set)
  //                       for(auto a13: x13_set)
  //                         for(auto a24: x24_set)
  //                           {
  //                             if (a13.irrep.mu()!=0)
  //                               continue;
  //                             u3::SU3 x12(a12.irrep), x34(a34.irrep), x13(a13.irrep), x24(a24.irrep);
  //                             int rho12_max(a12.tag), rho34_max(a34.tag), rho13_max(a13.tag), rho24_max(a24.tag);
  //                             MultiplicityTagged<u3::SU3>::vector x_set=u3::KroneckerProduct(x12,x34);
  //                             for(auto a: x_set)
  //                               {
  //                                 u3::SU3 x(a.irrep);
  //                                 int rho12_34_max=a.tag;
  //                                 int rho13_24_max=u3::OuterMultiplicity(x13,x24,x);
  //                                 for(int rho13_24=1; rho13_24<=rho13_24_max; ++rho13_24)
  //                                   for(int rho12_34=1; rho12_34<=rho12_34_max; ++rho12_34)
  //                                     for(int rho12=1; rho12<=rho12_max; ++rho12)
  //                                       for(int rho34=1; rho34<=rho34_max; ++rho34)
  //                                         for(int rho13=1; rho13<=rho13_max; ++rho13)
  //                                           for(int rho24=1; rho24<=rho24_max; ++rho24)
  //                                             {
  //                                               std::cout
  //                                               <<
  //                                               u3::Unitary9LambdaMu(
  //                                                 x1,    x2,    x12,    rho12,
  //                                                 x3,    x4,    x34,    rho34,
  //                                                 x13,   x24,   x,      rho13_24,
  //                                                 rho13, rho24, rho12_34
  //                                                 )
  //                                               <<"  "<<fmt::format("{} {} {} {}; {} {} {}; {} {} {},  {}",
  //                                               x1.Str(),x3.Str(),x.Str(),x4.Str(),x13.Str(),rho13,1,x34.Str(),rho34,1,  
  //                                               u3::U(x1,x3,x,x4,x13,rho13,1,x34,rho34,1))
  //                                               <<std::endl;
  //                                             }                            
  //                               }
  //                           }
  //                 }

  // // u3::SU3 x1(9,3);
  // // u3::SU3 x2(1,2);
  // // u3::SU3 x12(8,3);
  // // u3::SU3 x3(2,1);
  // // u3::SU3 x23(2,2);
  // // u3::SU3 x(8,2);

  // u3::SU3 x1(2,2);
  // u3::SU3 x2(2,0);
  // u3::SU3 x12(2,0);
  // u3::SU3 x3(1,0);
  // u3::SU3 x23(3,0);
  // u3::SU3 x(1,1);

  // // block access
  // std::cout << "U block test" << std::endl;
  // u3::UCoefLabels labels(x1, x2, x, x3, x12, x23);
  // u3::UCoefBlock block(labels);
  // int r12_max, r12_3_max, r23_max, r1_23_max;
  // std::tie(r12_max, r12_3_max, r23_max, r1_23_max) = block.Key();
  // std::cout << "multiplicities " << r12_max << " " << r12_3_max << " " << r23_max << " " << r1_23_max << std::endl;
  // for(int r12=1; r12<=r12_max; ++r12)
  //   for(int r12_3=1; r12_3<=r12_3_max; ++r12_3)
  //     for(int r23=1; r23<=r23_max; ++r23)
  //       for(int r1_23=1; r1_23<=r1_23_max; ++r1_23)
  //         std::cout << fmt::format("{} {} {} {}   {}  {}", r12, r12_3,r23, r1_23,block.GetCoef(r12,r12_3,r23,r1_23),
  //           u3::U(x1,x2,x,x3,x12,r12,r12_3,x23,r23,r1_23)
  //           ) << std::endl;
  //         std::cout << std::endl;


}
