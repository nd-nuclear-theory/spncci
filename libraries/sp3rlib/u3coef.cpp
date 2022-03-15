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

// #define WRU3_FUNCTION su3lib::wru3optimized_

namespace u3
{

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // W coefs
  ////////////////////////////////////////////////////////////////////////////////////////////////
  std::string WCoefLabels::Str() const
  {
    return fmt::format("({} {}; {} {}| {} {})",x1_.Str(),L1_,x2_.Str(),L2_,x3_.Str(),L3_);
  }

  double WCached(
                 u3::WCoefCache& cache,
                 const u3::SU3& x1, int kappa1, int L1, const u3::SU3& x2, int kappa2, int L2,
                 const u3::SU3& x3, int kappa3, int L3, int rho
                 )
  {
    const u3::WCoefLabels labels(x1,L1,x2,L2,x3,L3);

    // If block not cached.  Compute and cache. 
    if (cache.count(labels)==0)
      cache[labels]=u3::WCoefBlock(labels);

    const u3::WCoefBlock& block = cache.at(labels); 
    return block.GetCoef(kappa1, kappa2, kappa3, rho);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // U and Z coefficients
  ////////////////////////////////////////////////////////////////////////////////////////////////



  std::string UCoefLabels::Str() const
  {
    std::string label_str 
    = fmt::format("[{} {} {} {} {} {}",
      x1_.Str(),x2_.Str(),x_.Str(),x3_.Str(),x12_.Str(),x23_.Str()
      );
    
    return label_str;
  }



  double UCached(
                 u3::UCoefCache& cache,
                 const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, 
                 const u3::SU3& x3, const u3::SU3& x12,
                 int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23
                 )
  {
    const u3::UCoefLabels labels(x1,x2,x,x3,x12,x23);
    
    //If  coef block not found in cache, compute
    if (cache.count(labels)==0)
      {cache[labels]=u3::RecouplingCoefBlock(labels,RecouplingMode::kU);}

    const auto& block = cache.at(labels);

    return block.GetCoef(r12,r12_3,r23,r1_23);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////

  double ZCached(
                 u3::ZCoefCache& cache,
                 const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, 
                 const u3::SU3& x3, const u3::SU3& x12,
                 int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23
                 )
  {
    const u3::ZCoefLabels labels(x1,x2,x,x3,x12,x23);
    
    //If  coef block not found in cache, compute
    if (cache.count(labels)==0)
      {cache[labels]=u3::RecouplingCoefBlock(labels,RecouplingMode::kZ);}

    const auto& block = cache.at(labels);

    return block.GetCoef(r12,r12_3,r23,r1_23);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Interchange matrices and phases
  ////////////////////////////////////////////////////////////////////////////////////////////////

  std::string PhiCoefLabels::Str() const
  {return fmt::format("[{} {} {}]",x1_.Str(),x2_.Str(),x3_.Str());}

  double PhiCached(
         u3::PhiCoefCache& cache,
         const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x3, int rho1, int rho2
        )
  {
    const u3::PhiCoefLabels labels(x1,x2,x3);
    // #pragma omp critical (phi)
    {
      if (cache.count(labels)==0)
        cache[labels]=u3::PhiCoefBlock(labels);
    }
    const u3::PhiCoefBlock& block = cache.at(labels);
    return block.GetCoef(rho1, rho2);
  }






} // namespace
