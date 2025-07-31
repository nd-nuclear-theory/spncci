/****************************************************************
  u3u4_basis.cpp

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3st_basis.h"

#include "fmt/format.h"
#include "UNtoU3.h"
#include "cppitertools/unique_everseen.hpp"
#include "u3ncsm/u3u4_basis.h"
// debugging 
#include "utilities/utilities.h"


namespace un {

  std::vector<un::wSTshell> 
  generate_wSTshell_list(const un::UNU4Subspace& unu4_subspace)
  {
    // declare container 
    std::vector<un::wSTshell> wSTshell_list;    
    // make shared pointer to unu4 subspace for use in wSTshell constructor
    std::shared_ptr<un::UNU4Subspace> unu4_ptr = std::make_shared<un::UNU4Subspace>(unu4_subspace);
    // Iterate over all spins in [f] -> st and all U3 in [f]->w
    for(std::size_t index_st=0; index_st<unu4_subspace.getU4Subspace().size(); ++index_st)
      for(std::size_t index_w=0; index_w<unu4_subspace.size(); ++index_w )
      {
        wSTshell_list.emplace_back(unu4_ptr,index_w,index_st);
      }
    return wSTshell_list;
  }


  template<typename U>
  typename std::enable_if<std::is_same<U, u3::U3>::value, std::vector<MultiplicityTagged<u3::U3>>>::type
  KroneckerProduct(const U& u, const MultiplicityTagged<U>& v) {
      return u3::KroneckerProduct(u, v.irrep);
  }

  template<typename U>
  typename std::enable_if<std::is_same<U, HalfInt>::value, std::vector<HalfInt>>::type
  KroneckerProduct(const U& u, const U& v) {
      return am::ProductAngularMomenta(u, v);
  }

    template<typename T, typename U>
    void updateIntermediateCouplings(std::vector<std::vector<T>>& intermediate_couplings_w, const U& w)
    {
      std::vector<std::vector<T>> temp_intermediate_w;

      for(const auto& intermediate_w : intermediate_couplings_w)
      {
        for(const auto& wi : intermediate_w)
          // std::cout<<"im "<<w.Str()<<"  "<<wi.Str()<<std::endl;
        const auto& last = intermediate_w.back();
        // std::cout<<last.Str()<<std::endl;
        for(const auto& w_rho : un::KroneckerProduct(w,intermediate_w.back()))  
        {
          // Copy previous intermediate couplings and add new intermediate coupling.
          auto intermediate_w_new(intermediate_w);
          intermediate_w_new.push_back(w_rho);
          temp_intermediate_w.push_back(intermediate_w_new);
        }
      }
      // update intermediate w coupling container
      intermediate_couplings_w = temp_intermediate_w;                
    }


  std::map<u3::U3,std::vector<MultiplicityTagged<std::size_t>>>
  getIntermediateCouplingIndexed(
    const std::vector<std::vector<MultiplicityTagged<u3::U3>>>& intermediate_couplings_w,
    std::unordered_map<std::vector<u3::U3>,std::size_t>& intermediate_couplings_w_map
    )
  {
      std::map<u3::U3,std::vector<MultiplicityTagged<std::size_t>>> intermediate_w_indexed;
      
      for(const auto& intermediate_w : intermediate_couplings_w)
      {
        std::vector<u3::U3> temp_w_vec(intermediate_w.size());
        unsigned int degeneracy = 1; 
        for(int j=0; j<intermediate_w.size(); ++j)
        {
          temp_w_vec[j]=intermediate_w[j].irrep;
          degeneracy*=intermediate_w[j].tag;
        }

        // If not in master list, add
        // if(not intermediate_couplings_w_map.contains(temp_w_vec)) // need c++20, but c++20 not working with clang. 
        // Has issues with fmt::format.  
        if(intermediate_couplings_w_map.count(temp_w_vec)==0)
          intermediate_couplings_w_map[temp_w_vec]=intermediate_couplings_w_map.size();

        intermediate_w_indexed[intermediate_w.back().irrep].push_back({intermediate_couplings_w_map[temp_w_vec],degeneracy});

      }
      return intermediate_w_indexed;
  }

  std::map<HalfInt,std::vector<std::size_t>>
  getIntermediateCouplingIndexed(
    const std::vector<std::vector<HalfInt>>& intermediate_couplings,
    std::unordered_map<std::vector<HalfInt>,std::size_t>& intermediate_couplings_map
    )
  {
      std::map<HalfInt,std::vector<std::size_t>> intermediate_couplings_indexed;
      for(const auto& intermediate : intermediate_couplings)
      {     

        // If not in master list, add
        if(intermediate_couplings_map.count(intermediate)==0)
          intermediate_couplings_map[intermediate]=intermediate_couplings_map.size();

        intermediate_couplings_indexed[intermediate.back()].push_back(intermediate_couplings_map[intermediate]);

      }
      return intermediate_couplings_indexed;
  }

  wSTConfigurations::wSTConfigurations(unsigned int A, const unsigned int Nshell_max, unsigned int Nmax)
  {
    // Generate and save U4Space for given number of particles A
    u4space_ = u4::U4Space(A);
    // maximum number of particles in a given shell for given Nmax

    // generate all possible U(N) irreps for a given shell with quantum number Nshell with
    // U(4) conjugate [under Sn] that exists in u4space
    for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
    {
      unu4subspace_table_.resize(Nshell_max+1);
      unsigned int N = un::OmegaUN(Nshell); // N of U(N) corresponding to the shell.
      for(std::size_t index=0; index<u4space_.size(); ++index)
      {
        const auto& u4subspace = u4space_.GetSubspace(index);
        
        if(not u4subspace.f().HasSnConjugate(N)) continue;

        // if conjugate existis, add to list for given Nshell
        unu4subspace_table_[Nshell].push_back(UNU4Subspace(Nshell,u4space_.GetSubspacePtr(index)));
        

      }
    }


    // Distribute particles over Nshell_max+1 shells restricted to Nex<=Nmax
    // shell::ShellConfigurations is an ordered set of "typedef std::vector<int> Configuration"
    std::vector<shell::ShellConfigurations> configurations_by_Nex = shell::generate_configurations(A,Nmax,Nshell_max);
    
    // Get distributions over wSTShell for a given distribution of particles over Nshells
    
    // Lookup table to from specific wSTshell to it's index within the wSTshell_list_by_shell table 
    std::vector<std::unordered_map<un::wSTshell,std::size_t>> Nshell_wSTShell_map(Nshell_max+1);

    for(unsigned int Nex=Nmax%2; Nex<=Nmax; Nex+=2)
      {
        auto Nex_start=wSTshell_configurations_table_.size();
        
        for(const shell::Configuration& configuration : configurations_by_Nex[Nex])
        {
          // For a given distribution of particles, iterate over shells and idenitfy [f] irreps (UNU4Subspaces) which
          // can occur for a given number of particles in a given shell.
          // Then, iterate through specific wSTshell which occur in that unu4_subspace and create master list
          std::vector<std::vector<std::size_t>> Nshell_wSTShell_configuration_list(Nshell_max+1);
          unsigned int total_num_particles = 0;
          // utils::print_vector(configuration);
          for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
          {
            
            const auto& num_particles = configuration[Nshell];
            total_num_particles+=num_particles;
            auto& Nshell_wSTshell_list = Nshell_wSTShell_configuration_list[Nshell];
            for(const auto& unu4_subspace: unu4subspace_table_[Nshell])
            {
              if (unu4_subspace.num_particles()!=num_particles)
                continue;
              

              // get list of ([f],w,S,T for a given unu4_subspace)
              const auto temp = generate_wSTshell_list(unu4_subspace);
              for(const auto& wstshell : temp)
              {
                if (Nshell_wSTShell_map[Nshell].count(wstshell)==0)
                  Nshell_wSTShell_map[Nshell][wstshell]=Nshell_wSTShell_map[Nshell].size();

                // accumulate list of indices corresponding to wstshell which can appear for a given Nshell 
                // for the particular configuration.
                Nshell_wSTShell_configuration_list[Nshell].push_back(Nshell_wSTShell_map[Nshell][wstshell]);
              }
            }
          }
          
          assert(total_num_particles==A);
          
          // Compute number of wSTConfigurations which appear for a given particle distribution.
          // This is just the product of all possible combinations of different allowed wSTShell in each shell.
          int num_wst_configurations=1;
          for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
          {
            if (Nshell_wSTShell_configuration_list[Nshell].size()!=0)
              num_wst_configurations*=Nshell_wSTShell_configuration_list[Nshell].size();
          }

          
          // For the given distirbution over particles, construct list of possible wSTshell configurations 
          // which are represented as vector of length Nshell_max+1, where each element of vector is
          // the corresponding index of the particular wSTshell in the master list wSTshell_list_by_shell;
          std::vector<std::vector<std::size_t>> wSTshell_configurations(num_wst_configurations, std::vector<std::size_t>(Nshell_max+1,basis::kNone));
          
          std::size_t index_repeat = num_wst_configurations;
          std::size_t pattern_repeat = 1;
          for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
          {
            // number of times a given index is repeated
            index_repeat = index_repeat/Nshell_wSTShell_configuration_list[Nshell].size();
            
            std::size_t index=0;
            // iterate over number of times list of indices is repeated 
            for(std::size_t p=0; p<pattern_repeat; ++p)
              for(std::size_t wst_index=0; wst_index<Nshell_wSTShell_configuration_list[Nshell].size(); ++wst_index)
                for(std::size_t s=0; s<index_repeat; ++s)
                {
                  wSTshell_configurations[index][Nshell]=Nshell_wSTShell_configuration_list[Nshell][wst_index];

                  ++index;
                }

            pattern_repeat*=Nshell_wSTShell_configuration_list[Nshell].size();
            
          }


          // Add the list of configurations to the master list 
          wSTshell_configurations_table_.insert(wSTshell_configurations_table_.end(),wSTshell_configurations.begin(),wSTshell_configurations.end());
          
        } //end loop over Nex configurations
      auto Nex_end = wSTshell_configurations_table_.size();
      Nex_limits_[Nex] = {Nex_start,Nex_end}; 
      } // end loop over Nex

      // transfer wSTshell list from map into final container 
      wSTshell_list_by_shell_.resize(Nshell_max+1);
      for(unsigned int Nshell=0; Nshell<=Nshell_max; Nshell++)
      {
        // std::vector<std::unordered_map<un::wSTshell,std::size_t>> Nshell_wSTShell_map(Nshell_max+1);
        auto& wSTshell_list=wSTshell_list_by_shell_[Nshell];
        wSTshell_list.resize(Nshell_wSTShell_map[Nshell].size());
        for(const auto& [wst, index] : Nshell_wSTShell_map[Nshell])
          wSTshell_list[index]=wst;
      }

  } // end constructor

  


  void GenerateIntermediateCouplingsU3STSpace(
    const un::wSTConfigurations& wst_configurations,
    const unsigned int Nshell_max,
    const std::size_t config_index,
    std::vector<std::vector<MultiplicityTagged<u3::U3>>>& intermediate_couplings_w,
    std::vector<std::vector<HalfInt>>& intermediate_couplings_s,
    std::vector<std::vector<HalfInt>>& intermediate_couplings_t
    )
  {
      // Get quantum numbers for first two shells        
    const auto& wSTshell1 = wst_configurations.getwSTshell(config_index,0);// wSTshell_lists_by_shell[0][wSTshell_configuration[0]];
    const auto& wSTshell2 = wst_configurations.getwSTshell(config_index,1);//wSTshell_lists_by_shell[1][wSTshell_configuration[1]];
    
    auto u3_product = u3::KroneckerProduct(wSTshell1.w(),wSTshell2.w());
    
    for(const auto& w_rho : u3_product)
      intermediate_couplings_w.push_back({w_rho});
        
      
    const auto spin_product = am::ProductAngularMomenta(wSTshell1.S(),wSTshell2.S());
    for(const HalfInt& spin : spin_product) {intermediate_couplings_s.push_back({spin});}
    const auto isospin_product = am::ProductAngularMomenta(wSTshell1.T(),wSTshell2.T());
    for(const HalfInt& isospin : isospin_product) {intermediate_couplings_t.push_back({isospin});}

    // Loop over remaining shells and get complete set of all intermediate couplings. 
    for(unsigned int Nshell=2; Nshell<=Nshell_max; ++Nshell)
    {
      
      //Look up quantum numbers for Nshell
      const auto& wSTshell = wst_configurations.getwSTshell(config_index,Nshell);//wSTshell_lists_by_shell[Nshell][wSTshell_configuration[Nshell]];
      // Update intermediate couplings for u3, spin and isospin

      un::updateIntermediateCouplings(intermediate_couplings_w, wSTshell.w());
      un::updateIntermediateCouplings(intermediate_couplings_s, wSTshell.S());
      un::updateIntermediateCouplings(intermediate_couplings_t, wSTshell.T());
    }
  }


  U3STSpace::U3STSpace( const unsigned int A, const unsigned int Nshell_max, const unsigned int Nmax, HalfInt Tz)
  : BaseSpace{}
  {

  wst_configurations_=wSTConfigurations(A,Nshell_max,Nmax);

  std::unordered_map<std::vector<u3::U3>,std::size_t> intermediate_couplings_w_map;
  std::unordered_map<std::vector<HalfInt>,std::size_t> intermediate_couplings_st_map;

  // Loop over all configs and get intermediate couplings 
  for(unsigned int Nex=Nmax%2; Nex<=Nmax; Nex+=2)
  {
    const auto& Nex_limits = wst_configurations_.Nex_limits(Nex);  
    std::map<
        std::tuple<unsigned int, u3::SU3,HalfInt,HalfInt>, // omega,S,T
        std::map<
          std::size_t, // index for configuration in table table[Nex][index]
          std::vector<std::tuple<std::size_t,std::size_t,std::size_t,unsigned int>> //inter_w_index, inter_s_index,inter_t_index,degeneracy
          > 
    > U3STSpace_map;


    
    for(std::size_t config_index=Nex_limits.first; config_index<Nex_limits.second; ++config_index)
    {
      // inialize intermediate coupling list
      std::vector<std::vector<MultiplicityTagged<u3::U3>>> intermediate_couplings_w;
      std::vector<std::vector<HalfInt>> intermediate_couplings_s;
      std::vector<std::vector<HalfInt>> intermediate_couplings_t;

      GenerateIntermediateCouplingsU3STSpace(
        wst_configurations_,
        Nshell_max,
        config_index,
        intermediate_couplings_w,
        intermediate_couplings_s,
        intermediate_couplings_t
      );

      // Accumuate intermedaite couplings for in master containers and get vectors 
      // with indices of coupling labels in master container 
      auto intermediate_w_indexed
        =  getIntermediateCouplingIndexed(intermediate_couplings_w,intermediate_couplings_w_map);
      auto intermediate_s_indexed 
        = getIntermediateCouplingIndexed(intermediate_couplings_s,intermediate_couplings_st_map);
      auto intermediate_t_indexed 
        = getIntermediateCouplingIndexed(intermediate_couplings_t,intermediate_couplings_st_map);


      // create list of all combinations of wST with multiplicity for a given configuration and same 
      // in subsapce map container for use in constructing the u3u4subspaces.
      for(const auto& omega_list: intermediate_w_indexed)
        for(const auto& S_list : intermediate_s_indexed)
          for(const auto& T_list : intermediate_t_indexed)
            {
              if(Tz!=-1 and T_list.first<Tz) continue;

              auto& subspace_map = U3STSpace_map[{Nex,omega_list.first.SU3(),S_list.first,T_list.first}][config_index];

              for(const auto& [w_index,mult] : omega_list.second)
                for(const auto s_index : S_list.second)
                  for(const auto t_index : T_list.second)
                  {
                    subspace_map.emplace_back(w_index,s_index,t_index,mult);
                  }
            }
    }

    for(const auto& [label,subspace_map] : U3STSpace_map)
    {
      const auto&[Nex,lm,S,T] = label;
      PushSubspace(wSTSubspace(label,wst_configurations_,subspace_map));
    }

  }//end Nex

  // Transfer intermediat coupling information to saved containers 
  intermediate_couplings_w_.resize(intermediate_couplings_w_map.size());
  for(const auto& [w_list,index] : intermediate_couplings_w_map)
    intermediate_couplings_w_[index]=w_list;

  intermediate_couplings_st_.resize(intermediate_couplings_st_map.size());
  for(const auto& [st_list,index] : intermediate_couplings_st_map)
    intermediate_couplings_st_[index]=st_list;

  }
  // end constructor

}