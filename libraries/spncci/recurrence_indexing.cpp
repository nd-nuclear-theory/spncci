/****************************************************************
  recurrence_indexing.cpp

  Anna E. McCoy and Patrick J. Fasano
  Institute for Nuclear Theory  
  and University of Notre Dame

****************************************************************/

#include "spncci/recurrence_indexing.h"
#include "basis/basis.h"
#include <fstream>
#include <iostream>
#include <algorithm>
// #include "mcutils/parsing.h"
#include "fmt/format.h"
// #include "am/halfint_fmt.h"

// #include "sp3rlib/vcs.h"

namespace spncci
{
  SpNCCISpinSubspace::SpNCCISpinSubspace(
    const HalfInt& S, 
    const MultiplicityTagged<SpSn>::vector& spin_vector
  )
  {
    labels_ = {S};
    for(const auto& [SpSn,gamma_max] : spin_vector)
      PushStateLabels(SpSn,gamma_max); 
  }

  SpNCCISpinSpace::SpNCCISpinSpace(
    const u3::U3& sigma, 
    const std::map<HalfInt,MultiplicityTagged<SpSn>::vector>& spin_map
    )
  {
    labels_ = {sigma};
    for(const auto& [S,spin_vector] : spin_map)
      PushSubspace(SpNCCISpinSubspace(S,spin_vector));
  }
  
  SpNCCILGISpace::SpNCCILGISpace(const lgi::MultiplicityTaggedLGIVector& lgi_vector)
  {
    std::map<u3::U3,std::map<HalfInt,MultiplicityTagged<SpSn>::vector>> sigma_spins_map;
    for(const auto& [lgi,gamma_max] : lgi_vector)
      {
        const auto& [Nex,sigma,Sp,Sn,S] = lgi.Key();
        sigma_spins_map[sigma][S].emplace_back(SpSn{Sp,Sn},gamma_max);
      }
 
    for(const auto& [sigma,spin_map] : sigma_spins_map)
        PushSubspace(SpNCCISpinSpace(sigma,spin_map));
  }


}

