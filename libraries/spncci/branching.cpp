/****************************************************************
  branching.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching.h"

#include <sstream>
#include <iostream>

#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "spncci/unit_tensor.h"

namespace spncci
{

  SubspaceSpU3S::SubspaceSpU3S(const u3::U3S& omegaS,const BabySpNCCISpace& baby_spncci_space)
  {

    // save labels
    labels_ = omegaS;

    // scan BabySpNCCISpace for states to accumulate
    int substate_offset = 0;  // accumulated offset
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace = baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // short circuit if subspace not relevant
        if (!(omegaS==baby_spncci_subspace.omegaS()))
          continue;

        // push state
        PushStateLabels(baby_spncci_subspace.sigmaSPN(),baby_spncci_subspace.size());

        // record auxiliary state information
        state_gamma_max_.push_back(baby_spncci_subspace.gamma_max());
        state_baby_spncci_subspace_index_.push_back(baby_spncci_subspace_index);
      }

  }

  std::string SubspaceSpU3S::LabelStr() const
  {
    return labels().Str();
  }

  std::string SubspaceSpU3S::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpU3S state(*this,state_index);

        os << fmt::format(
            "  index {} omegaS {} sigmaSPN {} multiplicity {} offset {}",
            state_index,
            state.omegaS().Str(),state.sigmaSPN().Str(),
            state.multiplicity(),state.offset()
          ) << std::endl;
      }

    return os.str();
  }

  SpaceSpU3S::SpaceSpU3S(const BabySpNCCISpace& baby_spncci_space)
  {
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // create new subspace -- only if not already constructed for this (omega,S)
        u3::U3S omegaS = u3::U3S(baby_spncci_subspace.omega(),baby_spncci_subspace.S());
        if(ContainsSubspace(omegaS))
          continue;
        PushSubspace(SubspaceSpU3S(omegaS,baby_spncci_space));
      }
  }

  std::string SpaceSpU3S::DebugStr(bool show_subspaces) const
  {
    std::ostringstream os;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // set up alias
        const SubspaceSpU3S& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "index {} omegaS {} full_dimension {}",
            subspace_index,
            subspace.omegaS().Str(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }

}  // namespace
