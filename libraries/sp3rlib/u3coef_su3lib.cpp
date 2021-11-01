/****************************************************************
  u3coef.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include "sp3rlib/u3coef.h"

#include <algorithm>
#include <cassert>

#include "fmt/format.h"

#define WRU3_FUNCTION su3lib::wru3optimized_

namespace u3
{

  ////////////////////////////////////////////////////////////////
  // direct access to su3lib FORTRAN subroutines
  ////////////////////////////////////////////////////////////////

  namespace su3lib
  {

    const size_t MAX_K = 9;

    // Subroutines of original Fortran SU(3) library
    extern "C"
    {
      extern void wu3r3w_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[MAX_K][MAX_K][MAX_K][MAX_K]);
      extern void wru3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
      extern void wzu3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
      extern void wu39lm_(const int&, const int& , const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, double[], const int&);
      extern void blocks_(void);
    }

  } //namespace
  
  void U3CoefInit(int max_lambda_plus_mu)
  {
    if(max_lambda_plus_mu>=40)
      throw std::runtime_error("max (lambda + mu) must be < 40");

    su3lib::blocks_();
  }

  int WCoefBlock::CoefIndex(const int kappa1, const int kappa2, const int kappa3, const int rho) const
  {
    // equivalent to looking up w_array[k3-1][k2-1][k1-1][r0-1]
    // build up index, dimension by dimension
    int index = (rho-1);
    index = index * kappa2_max_ + (kappa2-1);
    index = index * kappa1_max_ + (kappa1-1);
    index = index * kappa3_max_ + (kappa3-1);
    return index;
  }


  WCoefBlock::WCoefBlock(const u3::WCoefLabels& labels)
  {
    // calculate multiplicities
    const auto&[x1,L1,x2,L2,x3,L3] = labels.Key();
    std::tie(kappa1_max_,kappa2_max_,kappa3_max_,rho_max_) = WMultiplicity(x1,L1,x2,L2,x3,L3);

    // zero initialize array
    double w_array[su3lib::MAX_K][su3lib::MAX_K][su3lib::MAX_K][su3lib::MAX_K];
    memset(w_array,0,sizeof(w_array));
    su3lib::wu3r3w_(x1.lambda(), x1.mu(), x2.lambda(), x2.mu(), x3.lambda(), x3.mu(), L1 , L2, L3, 1,1,1,1, w_array);
    int size=rho_max_*kappa1_max_*kappa2_max_*kappa3_max_;

    
    // coefs are in column-major Fortran array
    auto position=coefs_.begin();
    for(int rho=1; rho<=rho_max_; ++rho)
      for(int kappa2=1; kappa2<=kappa2_max_; ++kappa2)
        for(int kappa1=1; kappa1<=kappa1_max_; ++kappa1)
          for(int kappa3=1; kappa3<=kappa3_max_; ++kappa3)
            //Using row-major C to access column-major Fortran array
            coefs_.push_back(w_array[kappa3-1][kappa2-1][kappa1-1][rho-1]);
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
    coefs_.resize(r_max);
    
    // populate vector
    if (mode ==RecouplingMode::kZ)
    {
      // populate vector
      su3lib::wzu3optimized_(
        x1.lambda(), x1.mu(), x2.lambda(), x2.mu(), x.lambda(), x.mu(), 
        x3.lambda(), x3.mu(), x12.lambda(), x12.mu(), x23.lambda(), x23.mu(),
        r12_max_, r12_3_max_, r23_max_, r1_23_max_,
        coefs_.data(), r_max
      );
    }
    // mode is kU
    else
    {
      WRU3_FUNCTION(
        x1.lambda(), x1.mu(), x2.lambda(), x2.mu(), x.lambda(), x.mu(), 
        x3.lambda(), x3.mu(), x12.lambda(), x12.mu(), x23.lambda(), x23.mu(),
        r12_max_, r12_3_max_, r23_max_, r1_23_max_,
        coefs_.data(), r_max
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
    int rho_dummy=1;
    int dim=rho_max_*rho_max_;
    coefs_.resize(dim);
    su3lib::wzu3optimized_(
     x1.lambda(), x1.mu(), 0, 0, x3.lambda(), x3.mu(), x2.lambda(), x2.mu(),
     x1.lambda(), x1.mu(), x2.lambda(), x2.mu(),rho_dummy, rho_max_, rho_dummy, rho_max_,
     coefs_.data(), dim
     );
  }

////////////////////////////////////////////////////////////////
  double Unitary9LambdaMu(
                          const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x12, int r12,
                          const u3::SU3& x3,  const u3::SU3& x4,  const u3::SU3& x34, int r34,
                          const u3::SU3& x13, const u3::SU3& x24, const u3::SU3& x,   int r13_24,
                          int r13,     int r24,     int r12_34)
  {
    int r12_max=u3::OuterMultiplicity(x1,x2,x12);
    int r34_max=u3::OuterMultiplicity(x3,x4,x34);
    int r13_max=u3::OuterMultiplicity(x1,x3,x13);
    int r24_max=u3::OuterMultiplicity(x2,x4,x24);
    int r13_24_max=u3::OuterMultiplicity(x13,x24,x);
    int r12_34_max=u3::OuterMultiplicity(x12,x34,x);
    //std::vector<double> NLM_array;
    int index=(
               r12_max*(r34-1)
               +r12_max*r34_max*(r12_34_max-1)
               +r12_max*r34_max*r12_34_max*(r13_max-1)
               +r12_max*r34_max*r12_34_max*r13_max*(r24_max-1)
               +r12_max*r34_max*r12_34_max*r13_max*r24_max*(r13_24_max-1)
               );
    int r_max=r12_max*r34_max*r13_max*r24_max*r12_34_max*r13_24_max;
    std::vector<double> NLM_array(r_max);
    // double NLM_array[r_max];
    //NLM_array.resize(r_max);
    su3lib::wu39lm_(
                    x1.lambda(), x1.mu(), x2.lambda(), x2.mu(), x12.lambda(), x12.mu(),
                    x3.lambda(), x3.mu(), x4.lambda(), x4.mu(), x34.lambda(), x34.mu(),
                    x13.lambda(), x13.mu(), x24.lambda(), x24.mu(), x.lambda(), x.mu(),
                    NLM_array.data(),r_max
                    );
    return NLM_array[index];
  }



} // namespace
