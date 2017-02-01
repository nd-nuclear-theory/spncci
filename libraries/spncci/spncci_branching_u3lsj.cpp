/****************************************************************
  spncci_branching_u3lsj.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/spncci_branching_u3lsj.h"

#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"
#include "cppformat/format.h"

namespace spncci
{

  SubspaceLS::SubspaceLS(const int& L, const HalfInt& S,const SpaceU3S& u3s_space)
  {
    labels_=std::tuple<int,HalfInt>(L,S);
    int index=0;
    // iterate over U(3)xSU(2) irreps
    for(int subspace_index=0; subspace_index<u3s_space.size(); ++subspace_index)
      {
        SubspaceU3S subspace=u3s_space.GetSubspace(subspace_index);
        u3::U3 omega(subspace.omega());
        int kappa_max=u3::BranchingMultiplicitySO3(omega.SU3(),L);
        // if space contains S and omega can branch to L
        if(kappa_max>0 && subspace.S()==S)
        {
          //Construct subspace
          int dim=subspace.sector_dim();
          PushStateLabels(StateLabelsType(omega,kappa_max,index));
          // increment index 
          index+=kappa_max*dim;
        }
      }

    sector_size_=index;
  }

  std::string SubspaceLS::Str() const
  {
    return fmt::format("[{} {}]",L(),S());
  }

  std::string StateLS::Str() const
  {
    return fmt::format("[{} {}]",omega().Str(),kappa_max());
  }


  SpaceLS::SpaceLS(const SpaceU3S& u3s_space, HalfInt J)
  {
    // iterate over U(3)xSU(2) irreps
    for(int subspace_index=0; subspace_index<u3s_space.size(); ++subspace_index)
    // for(auto u3s_subspace : u3s_space)
      {
        const SubspaceU3S& u3s_subspace=u3s_space.GetSubspace(subspace_index);
        HalfInt S(u3s_subspace.S());
        // iterate through omega space
        u3::U3 omega(u3s_subspace.omega());
        // interate over possible L values
        for(int L=int(abs(S-J)); L<=(S+J); ++L)
          {
            if(lookup_.count(std::tuple<int,HalfInt>(L,S)))
              continue;
            SubspaceLS ls_subspace(L,S,u3s_space);
            PushSubspace(ls_subspace);
          }
       }
    // get dimension and starting index of last subspace
    const SubspaceType& subspace=GetSubspace(subspaces_.size()-1);
    HalfInt S=subspace.S();
    u3::U3 omega=std::get<0>(subspace.GetStateLabels(subspace.size()-1));
    int index=std::get<2>(subspace.GetStateLabels(subspace.size()-1));
    u3::U3S omegaS(omega,S);
    int dim=u3s_space.LookUpSubspace(omegaS).sector_dim();
    dimension_=dim+index;
  }


}  // namespace
