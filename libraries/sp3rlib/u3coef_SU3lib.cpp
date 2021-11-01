/****************************************************************
  u3coef.cpp

  SU(3) coupling coefficient wrappers for Dytrych and Langr
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "sp3rlib/u3coef.h"
#include <algorithm>
#include <cassert>
#include "fmt/format.h"
#include "su3.h"

namespace u3
{

  void U3CoefInit(int max_lambda_plus_mu)
  {
    int jjmax = 2*max_lambda_plus_mu;
    su3::init_thread(max_lambda_plus_mu, jjmax);
  }

  int WCoefBlock::CoefIndex(
    const int kappa1, const int kappa2, const int kappa3, const int rho) const
  {
    // SU3lib convention is mulitplicity index is 0<=index<multiplicity
    int index = (kappa3-1);
    index = index*kappa2_max_ +(kappa2-1);
    index = index*kappa1_max_ +(kappa1-1);
    index = index*rho_max_+(rho-1);
    return index;
  }
  

  WCoefBlock::WCoefBlock(const u3::WCoefLabels& labels)
  {
    // calculate multiplicities
    const auto&[x1,L1,x2,L2,x3,L3] = labels.Key();
    const auto&[lambda1,mu1]=x1.Key();
    const auto&[lambda2,mu2]=x2.Key();
    const auto&[lambda3,mu3]=x3.Key();
    std::tie(kappa1_max_,kappa2_max_,kappa3_max_,rho_max_) = WMultiplicity(x1,L1,x2,L2,x3,L3);
    
    coefs_.reserve(std::max(lambda1,mu1) * std::max(lambda2,mu2) * std::max(lambda3,mu3));
    
    su3::wu3r3w(lambda1,mu1,lambda2,mu2,lambda3,mu3,
      L1,L2,L3,rho_max_,kappa1_max_,kappa2_max_,kappa3_max_,coefs_
      );
  }



  // ////////////////////////////////////////////////////////////////

  // Same for U and Z
  int RecouplingCoefBlock::CoefIndex(int r12, int r12_3, int r23, int r1_23) const
  {
    int index = (r1_23-1);
    index = index * r23_max_ + (r23-1);
    index = index * r12_3_max_ + (r12_3-1);
    index = index * r12_max_ + (r12-1);
    return index;
  }

 RecouplingCoefBlock::RecouplingCoefBlock(const u3::UCoefLabels& labels, const u3::RecouplingMode mode)
  {
    // calculate multiplicities
    const auto& [x1,x2,x,x3,x12,x23] = labels.Key();
    std::tie(r12_max_,r12_3_max_,r23_max_,r1_23_max_) = UMultiplicity(x1,x2,x,x3,x12,x23);

    mode_=mode;

    // pre-size vector
    int r_max = r12_max_*r12_3_max_*r23_max_*r1_23_max_;
    assert(r_max > 0);

    // populate vector
    if (mode ==RecouplingMode::kZ)
    {
      su3::wzu3(
        x1.lambda(), x1.mu(), x2.lambda(), x2.mu(), x.lambda(), x.mu(),
        x3.lambda(), x3.mu(), x12.lambda(), x12.mu(), x23.lambda(), x23.mu(),
        r12_max_, r12_3_max_, r23_max_, r1_23_max_,coefs_
        );
    }
    // mode is kU
    else
    {
      su3::wru3(
        x1.lambda(), x1.mu(), x2.lambda(), x2.mu(), x.lambda(), x.mu(),
        x3.lambda(), x3.mu(), x12.lambda(), x12.mu(), x23.lambda(), x23.mu(),
        r12_max_, r12_3_max_, r23_max_, r1_23_max_,coefs_
      );

    }
  }


  PhiCoefBlock::PhiCoefBlock(const u3::PhiCoefLabels& labels)
  {
    // calculate multiplicities
    u3::SU3 x1,x2,x3;
    std::tie(x1,x2,x3) = labels.Key();
    rho_max_=u3::OuterMultiplicity(x1,x2,x3);
    // zero initialize array
    su3::wru3(x1.lambda(), x1.mu(), 0, 0, x3.lambda(), x3.mu(), x2.lambda(),
              x2.mu(), x1.lambda(), x1.mu(), x2.lambda(), x2.mu(), 1, rho_max_,
              1, rho_max_, coefs_);
  }




} // namespace 
