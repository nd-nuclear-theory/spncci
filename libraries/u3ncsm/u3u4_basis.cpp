/****************************************************************
  u3u4_basis.cpp

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3u4_basis.h"

// #include <unordered_map>
// #include <map>
#include "fmt/format.h"
#include "UNtoU3.h"
#include "cppitertools/unique_everseen.hpp"
// debugging 
#include "utilities/utilities.h"

namespace u4
{
  
  std::string U4Subspace::DebugStr(const std::string offset) const
  {
    std::string debug_str=fmt::format("{}{}\n",offset,f().Str());
    for (std::size_t i_state = 0; i_state < size(); ++i_state)
    {

      debug_str+=fmt::format("    {}[{} {}] {}\n",
          offset,
          GetState(i_state).S().Str(),
          GetState(i_state).T().Str(),
          GetState(i_state).beta_max()
        );
    }
    return debug_str;

  }


  std::string U4Space::DebugStr(const std::string offset) const
  {
    std::string debug_str="\n";
    for (std::size_t i = 0; i < size(); ++i)
    {
      debug_str+=GetSubspace(i).DebugStr(offset);
    }
    return debug_str;

  }
} //end namespace u4

namespace un
{
  std::vector<MultiplicityTagged<u3::U3>>
  BranchUNtoU3(const u4::U4&f, const unsigned int Nshell)
  {    

    UNtoU3<> gen;          
    gen.generateXYZ(Nshell);

    // Get U(N) conjugate tableau
    std::vector<MultiplicityTagged<u3::U3>> u3_irreps;

    u4::U4Conjugate conjugate_tableau(f, un::OmegaUN(Nshell));
    const auto& [a4,a3,a2,a1,a0]=conjugate_tableau.Key();  
    gen.generateU3Weights(a4,a3,a2,a1,a0);
    // iteration over generated U(3) weights
    for (const auto & pair : gen.multMap()) 
      {
        // get U(3) cartesian lables
        const auto & u3_tableau = pair.first;
        // get U(N)->U(3) branching multiplicity
        int alpha = gen.getLevelDimensionality(u3_tableau); 
        if (alpha)
          {
            const auto& [w1,w2,w3]=u3_tableau;
            u3_irreps.emplace_back(u3::U3(int(w1),int(w2),int(w3)),alpha); 
          }  
      }
    return u3_irreps;   
  }


  UNU4Subspace::UNU4Subspace(
      unsigned int Nshell,  
      const std::shared_ptr<const u4::U4Subspace>& u4subspace_ptr //pointer to corresponding U(4) subspace
    )
      : BaseDegenerateSubspace{{Nshell,u4subspace_ptr}}//space.GetSubspacePtr()
  {
    
    const u4::U4& f = (*u4subspace_ptr).f();
    assert(f.HasSnConjugate(un::OmegaUN(Nshell)));
    auto u3_irreps = un::BranchUNtoU3(f,Nshell);
    for(const auto& [w,alpha] : u3_irreps)
    {
     PushStateLabels({w}, alpha);
    }
  }


  U3U4IrrepTableType GetU3U4ForConfigurations(
    const shell::Configuration& configuration,
    const un::UNU4Space& unu4space,
    // const unsigned int Nex,
    const unsigned int Nshell_max
    )
  { 

    U3U4IrrepTableType u3u4_irreps_map;
    std::vector<std::vector<std::size_t>> Nshell_unu4subspace_list(Nshell_max+1);

    // For a given distribution of particles, iterate over shells and idenitfy [f] irreps (UNU4Subspaces) which
    // can occur for a given number of particles in a given shell.
    // Then, iterate through specific wSTshell which occur in that unu4_subspace and create master listunsigned int total_num_particles = 0;
    
    for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
    {
      const auto& num_particles = configuration[Nshell];
      // total_num_particles+=num_particles;

      auto& unu4subspace_list = Nshell_unu4subspace_list[Nshell];
      
      for(std::size_t subspace_index=0; subspace_index<unu4space.size(); ++subspace_index)
      {
        const auto& unu4subspace = unu4space.GetSubspace(subspace_index);
        if (unu4subspace.num_particles()==num_particles && unu4subspace.Nshell()==Nshell)
          unu4subspace_list.push_back(subspace_index);
      }
    }
    // std::cout<<"defined config"<<std::endl;
    int num_unu4_configurations=1;
    for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
    {
      if (Nshell_unu4subspace_list[Nshell].size()!=0)
        num_unu4_configurations*=Nshell_unu4subspace_list[Nshell].size();
    }

    // For the given distirbution over particles, construct list of possible unu4configurations 
    // which are represented as vector of length Nshell_max+1, where each element of vector is
    // the corresponding index of the subspace index of a particular unu4subpace
    std::vector<std::vector<std::size_t>> unu4_configurations(num_unu4_configurations, std::vector<std::size_t>(Nshell_max+1,basis::kNone));
    
    std::size_t index_repeat = num_unu4_configurations;
    std::size_t pattern_repeat = 1;
    for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
    {
      // number of times a given index is repeated
      index_repeat = index_repeat/Nshell_unu4subspace_list[Nshell].size();
      
      std::size_t index=0;
      // iterate over number of times list of indices is repeated 
      for(std::size_t p=0; p<pattern_repeat; ++p)
        for(std::size_t subspace_index=0; subspace_index<Nshell_unu4subspace_list[Nshell].size(); ++subspace_index)
          for(std::size_t s=0; s<index_repeat; ++s)
          {
            unu4_configurations[index][Nshell]=Nshell_unu4subspace_list[Nshell][subspace_index];

            ++index;
          }

      pattern_repeat*=Nshell_unu4subspace_list[Nshell].size();
    }
    
    
    // std::cout<<"finished setup"<<std::endl;
    // iterate over configurations and get coupled products 
    for(const auto& config : unu4_configurations)
      {
        std::unordered_map<u4::U4,unsigned int> u4_irreps_map = {{{0,0,0,0},1}};
        std::unordered_map<u3::U3,unsigned int> u3_irreps_map = {{{0,0,0},1}};
        for(unsigned int Nshell=0; Nshell<=Nshell_max; ++Nshell)
        {
          std::unordered_map<u4::U4,unsigned int> target_u4_irreps_map;
          std::unordered_map<u3::U3,unsigned int> target_u3_irreps_map;
          auto subspace_index  =config[Nshell];
          {
            const auto& subspace = unu4space.GetSubspace(subspace_index);
            for(const auto& [f_source,mult_source] : u4_irreps_map)
            {
              auto u4_product = u4::KroneckerProduct(f_source, subspace.f());
              for(const auto f_target_tagged : u4_product)
                target_u4_irreps_map[f_target_tagged.irrep]+=f_target_tagged.tag*mult_source;
            }
            // get U3 products 
            for(const auto& state : subspace)
            {
              
              for(const auto& [w_source,mult_source] : u3_irreps_map)
              {
                auto u3_product = u3::KroneckerProduct(w_source,state.w());
                for(const auto w_target_tagged : u3_product)
                  target_u3_irreps_map[w_target_tagged.irrep]+=w_target_tagged.tag*mult_source;
              }
            }
          }

          // Update map with new irreps 
          u4_irreps_map = target_u4_irreps_map;
          u3_irreps_map = target_u3_irreps_map;

        } //Nshell 
      
      //copy irreps into master lists
      for(const auto& [w,mult_w] : u3_irreps_map)
        for(const auto& [f,mult_f] : u4_irreps_map)
          u3u4_irreps_map[{w,f}]+=mult_w*mult_f;
      }//end config loop

return u3u4_irreps_map;
}//end function




}