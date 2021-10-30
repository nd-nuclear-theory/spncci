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
#include "su3.h"

namespace u3
{
  
  int WCoefIndex(const int k1, const int k2, const int k3, const int rho, const WMultiplicity& multiplicities)
  {
    const auto&[k1max,k2max,k3max,rhomax] = multiplicities;
    // int index = rho-1+(k1-1)*rhomax+(k2-1)*k1max*rhomax+(k3-1)*k2max*k1max*rhomax;
    // SU3lib convention is mulitplicity index is 0<=index<multiplicity
    index = (k3-1);
    index = index*k2max +(k2-1);
    index = index*k1max +(k1-1);
    index = index*rhomax+(rho-1);
    return index;

  }
  
  // double W(const u3::SU3& x1, int k1, int L1, const u3::SU3& x2, int k2, int L2, const u3::SU3& x3, int k3, int L3, int rho)
  //   {
      
  //     auto multiplicities = WMultiplicity(x1, L1, x2, L2,x3, L3);
  //     const auto&[k1max,k2max,k3max,rhomax] = multiplicities;
      
  //     std::vector<double> su3cgs;
  //     su3::wu3r3w(x1.lambda(),x1.mu(),x2.lambda(),x2.mu(),x3.lambda(),x3.mu(),L1,L2,L3,rhomax-1,k1max-1,k2max-1,k3max-1,su3cgs);
  //     // Coeffs returned in order, by k3, by k2, by k1,rho
  //     int index = WCoefIndex(k1,k2,k3,rho,multiplicities);

  //     return su3cgs[index];
  //   }


  // std::string WCoefLabels::Str() const
  // {
  //   return fmt::format("[{} {}; {} {}| {} {}]",x1_.Str(),L1_,x2.Str(),L2_,x3_.Str(),L3_);
  // }

  WCoefBlock::WCoefBlock(const u3::WCoefLabels& labels)
  {
    // calculate multiplicities
    const auto&[x1,L1,x2,L2,x3,L3] = labels.Key();
    std::tie(kappa1_max_,kappa2_max_,kappa3_max_,rho_max_) = WMultiplicity(x1,L1,x2,L2,x3,L3);
    su3::wu3r3w(x1.lambda(),x1.mu(),x2.lambda(),x2.mu(),x3.lambda(),x3.mu(),
      L1,L2,L3,rho_max_-1,kappa1_max_-1,kappa2_max_-1,kappa3_max_-1,coefs_
      );
  }


  // double WCoefBlock::GetCoef(int kappa1, int kappa2, int kappa3, int rho) const
  // {
  //   // validate multiplicity indices
  //   assert(
  //          (rho <= rho_max_)
  //          &&(kappa1 <= kappa1_max_)
  //          &&(kappa2 <= kappa2_max_)
  //          &&(kappa3 <= kappa3_max_)
  //          );
  //   // Lookup desired coefficient
  //   int index = WCoefIndex(kappa1, kappa2, kappa3, rho,{kappa1_max_,kappa2_max_,kappa3_max_,rho_max_});
  //   return coefs_[index];
  // }


  // double WCached(
  //   u3::WCoefCache& cache,
  //   const u3::SU3& x1, const int kappa1, const int L1,
  //   const u3::SU3& x2, const int kappa2, const int L2,
  //   const u3::SU3& x3, const int kappa3, const int L3, 
  //   int rho
  //   )
  // {
  //   double value;
  //   const u3::WCoefLabels labels(x1,L1,x2,L2,x3,L3);
  //   {
  //     if (cache.count(labels)==0)
  //       cache[labels]=u3::WCoefBlock(labels);
  //   }
  //   const u3::WCoefBlock& block = cache.at(labels); 
  //   return block.GetCoef(kappa1, kappa2, kappa3, rho);
  // }


} // namespace 
