/****************************************************************
  branching.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching.h"

#include <set>
#include <sstream>
#include <iostream>

#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "spncci/unit_tensor.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to U3S level
  ////////////////////////////////////////////////////////////////  

  SubspaceSpU3S::SubspaceSpU3S(const u3::U3S& omegaS,const BabySpNCCISpace& baby_spncci_space)
  {

    // save labels
    labels_ = omegaS;

    // scan BabySpNCCISpace for states to accumulate
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
    return omegaS().Str();
  }

  std::string SubspaceSpU3S::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpU3S state(*this,state_index);

        os << fmt::format(
            "  index {} omegaS {} sigmaSPN {} degeneracy {} offset {}",
            state_index,
            state.omegaS().Str(),state.sigmaSPN().Str(),
            state.degeneracy(),state.offset()
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
            "index {} omegaS {} size {} full_dimension {}",
            subspace_index,
            subspace.omegaS().Str(),subspace.size(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }

  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to LS level
  ////////////////////////////////////////////////////////////////  

  SubspaceSpLS::SubspaceSpLS(const spncci::LSLabels& ls_labels, const SpaceSpU3S& spu3s_space)
  {

    // save labels
    labels_ = ls_labels;

    // scan SpU3S space for states to accumulate
    for(int spu3s_subspace_index=0; spu3s_subspace_index<spu3s_space.size(); ++spu3s_subspace_index)
      {

        // set up alias
        const SubspaceSpU3S& spu3s_subspace = spu3s_space.GetSubspace(spu3s_subspace_index);

        // determine branching multiplicity to the specified L
        int kappa_max=u3::BranchingMultiplicitySO3(spu3s_subspace.x(),L());

        // short circuit if U3S subspace not relevant to current LS subspace
        if ((kappa_max==0)||(S()!=spu3s_subspace.S()))
          continue;

        // accumulate LS-branched states
        //
        // To obtain the mandated state ordering, we must pass through
        // the U3S subspace repeatedly, with kappa as the outermost
        // index.
        for (int kappa=1; kappa<=kappa_max; ++kappa)
          for (int spu3s_state_index=0; spu3s_state_index<spu3s_subspace.size(); ++spu3s_state_index)
            {
              // push state
              StateSpU3S spu3s_state(spu3s_subspace,spu3s_state_index);
              StateLabelsSpLS spls_state_labels = StateLabelsSpLS(
                  spu3s_state.omega(),
                  kappa,
                  spu3s_state.sigmaSPN()
                );
              PushStateLabels(spls_state_labels,spu3s_state.degeneracy());

              // record auxiliary state information
              state_gamma_max_.push_back(spu3s_state.gamma_max());
              state_baby_spncci_subspace_index_.push_back(spu3s_state.baby_spncci_subspace_index());
              state_spu3s_subspace_index_.push_back(spu3s_subspace_index);
            }

      }
  }

  std::string SubspaceSpLS::LabelStr() const
  {
    return fmt::format("({},{})",L(),S().Str());
  }

  std::string SubspaceSpLS::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpLS state(*this,state_index);

        os << fmt::format(
            "  index {} omegaS {} kappa {} sigmaSPN {} degeneracy {} offset {}",
            state_index,
            state.omegaS().Str(),state.kappa(),state.sigmaSPN().Str(),
            state.degeneracy(),state.offset()
          ) << std::endl;
      }

    return os.str();
  }


  SpaceSpLS::SpaceSpLS(const SpaceSpU3S& spu3s_space)
  {

    // collect (L,S) branchings
    std::set<LSLabels> ls_labels_set;
    for(int spu3s_subspace_index=0; spu3s_subspace_index<spu3s_space.size(); ++spu3s_subspace_index)
      {

        // set up alias
        const SubspaceSpU3S& spu3s_subspace = spu3s_space.GetSubspace(spu3s_subspace_index);

        // find branching L values
        u3::SU3 x = spu3s_subspace.omega().SU3();
        HalfInt S = spu3s_subspace.S();
        MultiplicityTagged<int>::vector branching = u3::BranchingSO3(x);

        // accumulate (L,S) pairs from U3S subspace
        for (const MultiplicityTagged<int>& l_kappa_max : branching)
          {
            int L = l_kappa_max.irrep;
            ls_labels_set.insert(LSLabels(L,S));
          }
      }
    
    // create subspaces
    for (const LSLabels& ls_labels : ls_labels_set)
      PushSubspace(SubspaceSpLS(ls_labels,spu3s_space));
  }

  SpaceSpLS::SpaceSpLS(const SpaceSpU3S& spu3s_space, HalfInt J)
  {

    // collect (L,S) branchings
    std::set<LSLabels> ls_labels_set;
    for(int spu3s_subspace_index=0; spu3s_subspace_index<spu3s_space.size(); ++spu3s_subspace_index)
      {

        // set up alias
        const SubspaceSpU3S& spu3s_subspace = spu3s_space.GetSubspace(spu3s_subspace_index);

        // find branching L values
        u3::SU3 x = spu3s_subspace.omega().SU3();
        HalfInt S = spu3s_subspace.S();
        HalfInt::pair l_range = am::ProductAngularMomentumRange(J,S);
        MultiplicityTagged<int>::vector branching = u3::BranchingSO3Constrained(x,l_range);

        // accumulate (L,S) pairs from U3S subspace
        for (const MultiplicityTagged<int>& l_kappa_max : branching)
          {
            int L = l_kappa_max.irrep;
            ls_labels_set.insert(LSLabels(L,S));
          }
      }
    
    // create subspaces
    for (const LSLabels& ls_labels : ls_labels_set)
      PushSubspace(SubspaceSpLS(ls_labels,spu3s_space));
  }

  std::string SpaceSpLS::DebugStr(bool show_subspaces) const
  {
    std::ostringstream os;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // set up alias
        const SubspaceSpLS& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "index {} (L,S) ({},{}) size {} full_dimension {}",
            subspace_index,
            subspace.L(),
            subspace.S().Str(),subspace.size(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }


}  // namespace
