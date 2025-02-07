/****************************************************************
  u4_test.cpp

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT

  + 01/15/25 (aem): Created.
****************************************************************/

#include "su4lib/u4.h"
#include "u3ncsm/u3u4_basis.h"
#include "u3ncsm/configurations.h"
#include "utilities/utilities.h"
#include "u3ncsm/u3st_basis.h"
int main()
{
if(false)
{
  u4::U4Subspace u4subspace(u4::U4(8u,3u,3u,0u));
  std::cout<<u4subspace.f().Str()<<std::endl;
  for(const auto& state : u4subspace)
    std::cout<<state.S()<<" "<<state.T()<<"  "<<state.beta_max()<<std::endl;
}
if(false)
{
  // std::cout<<u4subspace.DebugStr("")<<std::endl;

  unsigned int N_particles_max = 4;
  u4::U4Space u4space(N_particles_max);
  std::cout<<u4space.DebugStr()<<std::endl;

  for(unsigned int Nshell=0; Nshell<3; Nshell++)
  {
    unsigned int N = un::OmegaUN(Nshell);
    
    for(std::size_t index=0; index<u4space.size(); ++index)
    {
      
      const auto& u4subspace = u4space.GetSubspace(index);
      if(not u4subspace.f().HasSnConjugate(N)) continue;
      un::UNU4Subspace unu4subspace(Nshell,u4space.GetSubspacePtr(index));  

      std::cout<<Nshell<<"  "<<N<<"  "<<unu4subspace.N()<<"  "<<unu4subspace.f().Str()<<std::endl;
      std::cout<<"" ""<<unu4subspace.size()<<std::endl;
    }
  }
}
if (false)
{
  unsigned int A=10; 
  unsigned int Nshell_max=2;
  unsigned int Nmax=2;

  // std::vector<shell::ShellConfigurations> configurations = shell::generate_configurations(A,Nmax);

  un::wSTConfigurations wST_configurations(A,Nshell_max,Nmax);
  
  std::cout<<"num configs "<<wST_configurations.configurations_table().size()<<std::endl;
  for(const auto& config : wST_configurations.configurations_table())
  {
    utils::print_vector(config);
  }

  std::cout<<"print labels"<<std::endl;
  for(auto index=0; index<wST_configurations.configurations_table().size(); ++index)
  {
    const auto& labels = wST_configurations.configuration_labels(index);
    utils::print_vector(labels);
  }

  un::U3STSpace u3st_space(A, Nshell_max,Nmax);
  std::cout<<u3st_space.size()<<std::endl;
  for(const auto& subspace : u3st_space)
  {
    std::cout<<subspace.LabelStr()<<"  "<<subspace.dimension()<<std::endl;
  }

  std::cout<<"For specific nuclide with Tz=1"<<std::endl;
  un::U3STSpace u3st_space2(A, Nshell_max,Nmax,1);
  std::cout<<u3st_space2.size()<<std::endl;
  for(const auto& subspace : u3st_space2)
  {
    std::cout<<subspace.LabelStr()<<"  "<<subspace.dimension()<<std::endl;
  }
}

if(true)
{
  unsigned int A=6; 
  unsigned int Nshell_max=5;
  unsigned int Nmax=4;

  un::U3U4Configurations u3u4_configs(A,Nshell_max,Nmax);
  std::map<std::tuple<unsigned int, u3::SU3,HalfInt,HalfInt>,unsigned int> test_map1;
  for(const auto& [Nex,irreps] : u3u4_configs.u3_u4_irreps())
  {
    for(const auto& [labels,mult] :irreps)
    {
      const auto&[w,f]=labels;

      // fmt::print("{}  {} {}  {}\n",Nex, w,f,mult);
      auto weights = u4::GenerateU4toSU2SU2Weights(f);
      for(const auto& [S,T,multiplicity] : weights)
        test_map1[{Nex,w.SU3(),S,T}]+=mult*multiplicity;
        // fmt::print("   {} {}  {}\n",S,T,multiplicity);

    }
  }
  fmt::print("----------------------\n");
  std::map<std::tuple<unsigned int, u3::SU3,HalfInt,HalfInt>,unsigned int> test_map2;
  un::U3STSpace u3st_space(A, Nshell_max,Nmax);
  std::cout<<u3st_space.size()<<std::endl;
  for(const auto& subspace : u3st_space)
  {
    test_map2[{subspace.Nex(),subspace.SU3(), subspace.S(), subspace.T()}]+=subspace.dimension();
    // std::cout<<subspace.LabelStr()<<"  "<<subspace.dimension()<<std::endl;
  }

  // std::cout<<test_map1.size()<<"  "<<test_map2.size()<<std::endl;
  for(const auto& [labels,mult]: test_map1)
  {
    if (test_map2.count(labels)==0)
    {
      const auto& [Nex,w,S,T]=labels;
      fmt::print("irrep {} {} {} {}  not found\n",Nex,w,S,T);
      continue;
    }
    if(test_map2[labels]!=mult)
    {
      const auto& [Nex,w,S,T]=labels;
      fmt::print("irrep {} {} {} {}  has different multiplicity\n",Nex,w,S,T);
      fmt::print("   mult: {} {}\n",test_map2[labels],mult);

    }

  }

}
}