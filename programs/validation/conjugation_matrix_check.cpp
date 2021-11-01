/****************************************************************
  conjugation_matrix_check.cpp

  Checking definition of conjugation matrix as function of identity coupled coefficient                                  
  
  Anna E. McCoy
  TRIUMF
  
  SPDX-License-Identifier: MIT

  08/28/18 (aem): Created.
****************************************************************/
#include "fmt/format.h"

//#include "am/halfint.h"
// #include "am/wigner_gsl.h"
// #include "utilities/utilities.h"
// #include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
// #include <map>


void CheckIdentityCoupling(const u3::SU3& x)
  {
    MultiplicityTagged<int>::vector branch1=u3::BranchingSO3(x);
    for(auto& L_kappa : branch1)
      {
        int L=L_kappa.irrep;
        int kappa_max=u3::BranchingMultiplicitySO3(x,L);
        for(int kappa=1; kappa<=kappa_max; ++kappa)
          for(int kappap=1; kappap<=kappa_max; ++kappap)
            {
              double coef=u3::W(x, kappa,L,u3::SU3(0,0),1,0,x,kappap,L,1);
              std::cout<<fmt::format("Identity coupling with {} {} : {} {}     {:10f}",x.Str(),L,kappa,kappap, coef)<<std::endl;
            }
      }
  }

void CheckZeroCoupledClebsch(const u3::SU3& x)
  {
    MultiplicityTagged<int>::vector branch1=u3::BranchingSO3(x);
    for(auto& L_kappa : branch1)
      {
        int L=L_kappa.irrep;
        int kappa_max=L_kappa.tag;
        // int kappa_max=u3::BranchingMultiplicitySO3(x,L);
        for(int kappa=1; kappa<=kappa_max; ++kappa)
          for(int kappap=1; kappap<=kappa_max; ++kappap)
            {
              double coef=u3::W(x, kappa,L,u3::Conjugate(x),kappap,L,u3::SU3(0,0),1,0,1);
              double coef_expected=std::sqrt((2.*L+1)/u3::dim(x))*ParitySign(u3::ConjugationGrade(x));
              double diff=coef/coef_expected;
              std::cout<<fmt::format("{} {} {}:  {} {}     {:10f}      {:10f}     {:10f}",x.Str(),L,kappa_max,kappa,kappap,coef,coef_expected,diff)<<std::endl;
            }
      }

  }


void CheckReversalSymmetry(const u3::SU3& x1, const u3::SU3& x2)
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
                    for(int i=0; i<product.size(); ++i)
                      {
                        u3::SU3 x(product[i].irrep);
                        // std::cout<<x.Str()<<std::endl;
                        int rho_max=product[i].tag;
                        int kappa_max=u3::BranchingMultiplicitySO3(x,L);
                        for(int rho=1; rho<=rho_max; ++rho)  
                          for(int kappa=1; kappa<=kappa_max; ++kappa)
                          {
                            double coef=u3::W(x1,kappa1,L1,x2,kappa2,L2,x,kappa,L,rho);
                            // double coef_rev=u3::W(x,kappa,L,u3::Conjugate(x2),kappa2p,L2,x1,kappa1,L1,rho);

                            double coef_rev2=0.0;

                            for(int kappa2p=1; kappa2p<=kappa2_max; ++kappa2p)
                              {
                                coef_rev2+=u3::W(u3::Conjugate(x2),kappa2p,L2,x2,kappa2,L2,u3::SU3(0,0),1,0,1)
                                              *ParitySign(x1.lambda()+x1.mu()+x2.lambda()+x2.mu()+x.lambda()+x.mu())
                                              *ParitySign(L1+L2+L)
                                              *std::sqrt(u3::dim(x2)/(2.*L2+1))
                                              *std::sqrt((2.*L1+1)*u3::dim(x)/(2.*L+1)/u3::dim(x1))
                                              *u3::W(x,kappa,L,u3::Conjugate(x2),kappa2p,L2,x1,kappa1,L1,rho);
                              }

                            double coef_rev=ParitySign(u3::ConjugationGrade(x2))
                                              *ParitySign(x1.lambda()+x1.mu()+x2.lambda()+x2.mu()+x.lambda()+x.mu())
                                              *ParitySign(L1+L2+L)
                                              *std::sqrt((2.*L1+1)*u3::dim(x)/(2.*L+1)/u3::dim(x1))
                                              *u3::W(x,kappa,L,u3::Conjugate(x2),kappa2,L2,x1,kappa1,L1,rho);

                            if(fabs(coef-coef_rev2)>1e-6)
                            {
                            std::cout<<fmt::format("{} {} {}; {} {} {} | {} {} {}  {}    {:10f}   {}",
                                          x1.Str(),kappa1,L1,x2.Str(),kappa2,L2,x.Str(),kappa,L,rho, coef, kappa2_max)<<std::endl;
                            // std::cout<<fmt::format("{} {} {}; {} {} {} | {} {} {}  {}    {:10f}",
                            //               x.Str(),kappa,L,u3::Conjugate(x2).Str(),kappa2,L2,x1.Str(),kappa1,L1,rho,coef_rev)<<std::endl;
                            std::cout<<fmt::format("{} {} {}; {} {} {} | {} {} {}  {}    {:10f}",
                                          x.Str(),kappa,L,u3::Conjugate(x2).Str(),kappa2,L2,x1.Str(),kappa1,L1,rho,coef_rev2)<<std::endl<<std::endl;
                            }


                            // double coef_rev=ParitySign(u3::ConjugationGrade(x2))*std::sqrt((2.*L1+1)*u3::dim(x)/(2.*L+1)/u3::dim(x1))
                            //                   *u3::W(x,kappa,L,u3::Conjugate(x2),kappa2p,L2,x1,kappa1,L1,rho);

                          }
                      }  
                  }
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
  for(int lambda=lm_min; lambda<=lm_max; ++lambda)
    for(int mu=mu_min; mu<=mu_max; ++mu)
      {
        u3::SU3 x(lambda,mu);
        CheckIdentityCoupling(x);
      }

  for(int lambda=lm_min; lambda<=lm_max; ++lambda)
    for(int mu=mu_min; mu<=mu_max; ++mu)
      {
        u3::SU3 x(lambda,mu);
        CheckZeroCoupledClebsch(x);
      }


  for(int lambda1=lm_min; lambda1<=lm_max; ++lambda1)
    for(int mu1=mu_min; mu1<=mu_max; ++mu1)
      for(int lambda2=lm_min; lambda2<=lm_max; ++lambda2)
        for(int mu2=mu_min; mu2<=mu_max; ++mu2)
          {
            // if(lambda2==mu2)
            //   continue;
            u3::SU3 x1(lambda1,mu1);
            u3::SU3 x2(lambda2,mu2);
            CheckReversalSymmetry(x1, x2);
          }

  return 0;
}
