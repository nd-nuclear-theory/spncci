/****************************************************************
  check_Wcoef_orthogonality.cpp

  Checking orthogonality of SU(3) coupling coefficients  
                                  
  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/10/16 (aem): Created.
****************************************************************/
#include "fmt/format.h"

//#include "am/halfint.h"
// #include "am/wigner_gsl.h"
// #include "utilities/utilities.h"
// #include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
// #include <map>


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
                std::cout<<fmt::format("({},{}) x ({},{})",l1,m1,l2,m2)<<std::endl;
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

int main(int argc, char **argv)
{
  // initialize su3lib
  int max_lambda_plus_mu=39;
  u3::U3CoefInit(max_lambda_plus_mu);

  if(argc<4)
  {
  	std::cout<<"Syntax: lambda_min lambda_max mu_min mu_max"<<std::endl;
    return 0;
  }

  int lm_min=std::stoi(argv[1]);
  int lm_max=std::stoi(argv[2]);
  int mu_min=std::stoi(argv[3]);
  int mu_max=std::stoi(argv[4]);

  //test orthogonality of W coefficients 
  std::cout<<"Testing orthogonality"<<std::endl;
  TestOrthogonalityW(lm_min,lm_max, mu_min,mu_max);
  return 0;
}
