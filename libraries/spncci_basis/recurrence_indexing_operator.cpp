/****************************************************************
  recurrence_indexing_relative_operator.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "spncci_basis/recurrence_indexing_operator.h"

#include "fmt/format.h"
#include <cppitertools/imap.hpp>
#include <cppitertools/unique_everseen.hpp>



namespace relative
{
namespace spatial
{
  std::unordered_map<u3::SU3,std::vector<std::tuple<unsigned int,unsigned int>>>
  GenerateSpatialOperators(
    const int N0,
    const uint8_t parity_bar,
    const OperatorParameters& operator_parameters
  )
  {
    int Nbar_min=parity_bar;
    int Nbar_max=operator_parameters.Nbar_max;
    const auto& Allowed_w0_values=operator_parameters.Allowed_w0_values;
    bool restrict_x0 = Allowed_w0_values.size()>0;
    //////////////////////////////////////////////////////////////////////////////
    // int Nbar_min = unit_tensor_constraints.parity_bar;
    std::unordered_map<u3::SU3,std::vector<std::tuple<unsigned int,unsigned int>>> x0_Nbar_pairs;

    for (int Nbar = Nbar_min; Nbar <= Nbar_max; Nbar += 2)
      {
        int Nbar_p = N0 + Nbar;
        if (Nbar_p <= Nbar_max & Nbar_p >= 0)
        {
          MultiplicityTagged<u3::SU3>::vector possible_x0 =
              u3::KroneckerProduct({Nbar_p, 0}, {0, Nbar});

          for (const auto& [x0, rho] : possible_x0)
            {
              if(!restrict_x0  || Allowed_w0_values.count({N0,x0}))
                x0_Nbar_pairs[x0].push_back({Nbar, Nbar_p});

            }
        }
      }
    return x0_Nbar_pairs;
  }


  OperatorSpace::OperatorSpace(const OperatorParameters& operator_parameters)
  : BaseDegenerateSpace{0}
  {
    // For a two-body operator, the maximum value of S0 possible is 2.
    // S0 the minimum possible L0 that will couple to with S0 to J0 is
    L0min_=std::max(0,int(operator_parameters.J0)-2);

    std::vector<HalfInt> Allowed_N0_values;
    const auto& Allowed_w0_values = operator_parameters.Allowed_w0_values;
    if(operator_parameters.Allowed_w0_values.size()>0)
      {
        // Get list of possible N0 values
        for( const HalfInt& N0 : iter::unique_everseen(iter::imap([] (u3::U3 w) {return w.N();}, Allowed_w0_values)))
          Allowed_N0_values.push_back(N0);
      }

    for(uint8_t parity_bar : {0,1})
      for(int N0=0; N0<=operator_parameters.Nbar_max; N0+=2)
        {
          auto x0_Nbar_pairs
            = GenerateSpatialOperators(N0,parity_bar,operator_parameters);

          for(const auto& [x0,Nbar_pairs] : x0_Nbar_pairs)
            {

              int degeneracy=0;

              std::array<std::size_t,5> L0_offsets;
              const auto& Allowed_L0_values = operator_parameters.Allowed_L0_values;
              for(int L0=0; L0<L0min_+5; ++L0)
                {
                  // Save offset information
                  L0_offsets[L0-L0min_]=degeneracy;

                  if(Allowed_L0_values.size()==0 || Allowed_L0_values.count(L0)!=0)
                    {
                      auto kappa0_max = u3::BranchingMultiplicitySO3(x0,L0);

                      // increment degeneracy
                      degeneracy += kappa0_max;
                    }
                }

              if(degeneracy==0)
                continue;

              auto subspace = OperatorU3Subspace(parity_bar,N0,x0,Nbar_pairs);
              if(subspace.dimension()==0)
                continue;

              PushSubspace(std::move(subspace),degeneracy);
              L0_offsets_.push_back(L0_offsets);
            }
        }
  }





  std::string OperatorSpace::DebugStr() const
  {
    std::string debug_str="";

    for (int i_subspace=0; i_subspace<size(); ++i_subspace)
    {
      const auto& subspace = GetSubspace(i_subspace);
      debug_str+=fmt::format("Subspace: {}\n",i_subspace);
      debug_str+=fmt::format("  ----------\n");
      for(int l=0; l<5; ++l)
        {
          unsigned int L0 = L0min_+l;
          debug_str+=fmt::format("   L0 {}: kappa0_max: {}\n",
          L0,
          Getkappa0max(i_subspace,L0)
          );
        }
      debug_str+=fmt::format("  ----------\n");
      debug_str+=fmt::format("  labels: {} degeneracy: {} size: {}  dimensions: {}\n",
        subspace.LabelStr(),
        GetSubspaceDegeneracy(i_subspace),
        subspace.size(),
        subspace.dimension()
        );

      for(int i_state=0; i_state<subspace.size(); ++i_state)
        {
          debug_str+=fmt::format("    {}\n",subspace.GetState(i_state).LabelStr());
        }
    }

    return debug_str;
  }

}//end spatial namespace

namespace spin
{

  std::string OperatorState::LabelStr() const
    {
      return fmt::format("[{} {} {} {}]\n",Sbar(),Sbarp(),Tbar(),Tbarp());
    }


  std::string OperatorT0Subspace::DebugStr() const
  {
    std::string debug_str
      =fmt::format("    T0subspace: {} size: {}  dimension: {}\n",
        T0(),size(),dimension()
      );

    for (int i_state=0; i_state<size(); ++i_state)
    {
      const auto& state = GetState(i_state);
      debug_str+=fmt::format("     {}",state.LabelStr());
    }

    return debug_str;
  }

  std::string OperatorS0Subspace::DebugStr() const
  {
    std::string debug_str
      =fmt::format("  S0Subspace: S0: {} parity_bar: {} size: {}  dimension: {}\n",
        S0(),parity_bar(),size(),dimension()
      );

    for (int i_subspace=0; i_subspace<size(); ++i_subspace)
    {
      const auto& subspace = GetSubspace(i_subspace);
      debug_str+=subspace.DebugStr();
    }

    return debug_str;
  }


  std::string OperatorSpace::DebugStr() const
  {
    std::string debug_str="";

    for (int i_subspace=0; i_subspace<size(); ++i_subspace)
    {
      const auto& subspace = GetSubspace(i_subspace);
      debug_str+=subspace.DebugStr();
    }

    return debug_str;
  }


}//end spin namespace



  OperatorSectors::OperatorSectors(
    const relative::spatial::OperatorSpace& spatial_operator_space,
    const relative::spin::OperatorSpace& spin_operator_space)
  {
    for(int i_spatial=0; i_spatial<spatial_operator_space.size(); ++i_spatial)
      for(int i_spin=0; i_spin<spin_operator_space.size(); ++i_spin)
        {
          const auto& U3_subspace = spatial_operator_space.GetSubspace(i_spatial);
          const auto& S0_subspace = spin_operator_space.GetSubspace(i_spin);

          auto  spatial_ptr = spatial_operator_space.GetSubspacePtr(i_spatial);
          auto spin_ptr = spin_operator_space.GetSubspacePtr(i_spin);

          unsigned int spatial_degeneracy
            = spatial_operator_space.GetSubspaceDegeneracy(i_spatial);
          unsigned int spin_degeneracy = 1;

          if(U3_subspace.parity_bar()==S0_subspace.parity_bar())
            PushSector(
              OperatorSector(
                i_spatial,
                i_spin,
                spatial_degeneracy,
                spin_degeneracy,
                spatial_ptr,
                spin_ptr
              )
              // i_spatial,i_spin
            );
        }
  }


}// end relative namespace





