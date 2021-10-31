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
// namespace su3lib
//   {

//     const size_t MAX_K = 9;

//     // Subroutines of original Fortran SU(3) library
//     extern "C"
//     {
//       // extern void wu3r3w_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[MAX_K][MAX_K][MAX_K][MAX_K]);
//       extern void wru3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
//       extern void wzu3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
//       extern void wu39lm_(const int&, const int& , const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, double[], const int&);
//       extern void blocks_(void);
//     }

//   } //namespace

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // W coefs
  ////////////////////////////////////////////////////////////////////////////////////////////////
  std::string WCoefLabels::Str() const
  {
    return fmt::format("({} {}; {} {}| {} {})",x1_.Str(),L1_,x2_.Str(),L2_,x3_.Str(),L3_);
  }

  double WCoefBlock::GetCoef(const int kappa1, const int kappa2, const int kappa3, const int rho) const
  {
    
    // validate multiplicity indices
    assert(
           (rho <= rho_max_)
           &&(kappa1 <= kappa1_max_)
           &&(kappa2 <= kappa2_max_)
           &&(kappa3 <= kappa3_max_)
           );

    // calculate index into block
    int index = CoefIndex(kappa1,kappa2,kappa3,rho);
    return coefs_[index];
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

  double W(const u3::SU3& x1, int k1, int L1, const u3::SU3& x2, int k2, int L2, const u3::SU3& x3, int k3, int L3, int r0)
  {
    WCoefLabels labels(x1,L1,x2,L2,x3,L3);
    auto block = WCoefBlock(labels);
    return block.GetCoef(k1,k2,k3,r0);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // U and Z coefficients
  ////////////////////////////////////////////////////////////////////////////////////////////////

  u3::UMultiplicityTuple UMultiplicity(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x,
                                       const u3::SU3& x3, const u3::SU3& x12, const u3::SU3& x23)
  {
    int r12_max = u3::OuterMultiplicity(x1,x2,x12);
    int r12_3_max = u3::OuterMultiplicity(x12,x3,x);
    int r23_max = u3::OuterMultiplicity(x2,x3,x23);
    int r1_23_max = u3::OuterMultiplicity(x1,x23,x);
    return UMultiplicityTuple(r12_max,r12_3_max,r23_max,r1_23_max);
  }

  u3::UMultiplicityTuple ZMultiplicity(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x,
                                       const u3::SU3& x3, const u3::SU3& x12, const u3::SU3& x23)
  {
    int r12_max = u3::OuterMultiplicity(x1,x2,x12);
    int r12_3_max = u3::OuterMultiplicity(x12,x3,x);
    int r23_max = u3::OuterMultiplicity(x2,x3,x23);
    int r1_23_max = u3::OuterMultiplicity(x1,x23,x);
    return ZMultiplicityTuple(r12_max,r12_3_max,r23_max,r1_23_max);
  }


  u3::WMultiplicityTuple WMultiplicity(
    const u3::SU3& x1, const int L1, 
    const u3::SU3& x2, const int L2,
    const u3::SU3& x3, const int L3
    )
  {
    int rho_max=u3::OuterMultiplicity(x1,x2,x3);
    int kappa1_max=u3::BranchingMultiplicitySO3(x1,L1);
    int kappa2_max=u3::BranchingMultiplicitySO3(x2,L2);
    int kappa3_max=u3::BranchingMultiplicitySO3(x3,L3);

    return WMultiplicityTuple(kappa1_max, kappa2_max, kappa3_max, rho_max);
  }



  std::string UCoefLabels::Str() const
  {
    std::string label_str 
    = fmt::format("[{} {} {} {} {} {}",
      x1_.Str(),x2_.Str(),x_.Str(),x3_.Str(),x12_.Str(),x23_.Str()
      );
    
    return label_str;
  }

  double UCoefBlock::GetCoef(int r12, int r12_3, int r23, int r1_23) const
  {
    // validate multiplicity indices
    assert(
           (r12 <= r12_max_)
           &&(r12_3 <= r12_3_max_)
           &&(r23 <= r23_max_)
           &&(r1_23 <= r1_23_max_)
           );
    int index = CoefIndex(r12,r12_3,r23,r1_23);
    return coefs_[index];
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
      {cache[labels]=u3::UCoefBlock(labels);}

    const u3::UCoefBlock& block = cache.at(labels);

    return block.GetCoef(r12,r12_3,r23,r1_23);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////

  std::string ZCoefLabels::Str() const
  {
    std::string label_str 
    = fmt::format("[{} {} {} {} {} {}",
      x1_.Str(),x2_.Str(),x_.Str(),x3_.Str(),x12_.Str(),x23_.Str()
      );
    
    return label_str;
  }

  double ZCoefBlock::GetCoef(int r12, int r12_3, int r23, int r1_23) const
  {
    // validate multiplicity indices
    assert(
           (r12 <= r12_max_)
           &&(r12_3 <= r12_3_max_)
           &&(r23 <= r23_max_)
           &&(r1_23 <= r1_23_max_)
           );
    int index = CoefIndex(r12,r12_3,r23,r1_23);
    return coefs_[index];
  }

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
      {cache[labels]=u3::ZCoefBlock(labels);}

    const u3::ZCoefBlock& block = cache.at(labels);

    return block.GetCoef(r12,r12_3,r23,r1_23);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Interchange matrices and phases
  ////////////////////////////////////////////////////////////////////////////////////////////////

  std::string PhiCoefLabels::Str() const
  {
    std::ostringstream ss;

    ss << "[" << x1_.Str()<< x2_.Str() << x3_.Str() << "]";
    return ss.str();
  }



  double PhiCoefBlock::GetCoef(int rho1, int rho2) const
  {
    // validate multiplicity indices
    assert((rho1 <= rho_max_)&&(rho2<=rho_max_));

    int index = (rho2-1);
    index = index * rho_max_ + (rho1-1);
    // retrieve entry
    double value = cache_[index];

    return value;
  }

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


  double Phi(const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x3, int r, int rp)
  // Phi phase factor that arrises in chainging the coupling order of SU(3) irreps
  {
    
    return PhiCoefBlock({x1,x2,x3}).GetCoef(r,rp);
  }



} // namespace
