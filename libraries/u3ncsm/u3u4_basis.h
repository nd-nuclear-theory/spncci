/****************************************************************
  u4.h

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT

  + 01/15/25 (aem): Created.
****************************************************************/
#ifndef U3U4_BASIS_H_
#define U3U4_BASIS_H_

#include <unordered_map>
#include <map>
#include <functional>
#include <omp.h>

#include "basis/basis.h"
#include "basis/degenerate.h" 
#include "sp3rlib/u3.h"
#include "su4lib/u4.h"
#include "su4lib/u4_branching.h"
#include "am/halfint.h"
#include "utilities/nuclide.h"
#include "u3ncsm/configurations.h"

// U4Space
// Generate list of all possible U(4) configurations is N=0...A particles.  
//  -> U4Subpsace [f], has accessor giving num particles 
//     -> SpinSubspace: Generate all possible beta[ST] in each U(4) irrep.
//  

  // UNU4Space(U4Space, [N_shell]) 
  // For each N_shell, get corresponding U(N) 
  //    ->UNU4subspace 
  //      For each U4 in U4Space, if has valid Sn_conjugate for U(N), add subspace labeled by
  //        u4_subspace_index and has pointer to U4Subspace. 
  //        -> U3Subspa
  // Subspace constructor


namespace u4{
  // Funcation to generate list of U(4) tableau with A particles 
  inline std::vector<u4::U4>
  GenerateU4Tableaux(unsigned int A)
    {
      std::vector<u4::U4> tableaux;
      // Generate set of tableau for A particles 
      // Since f1>=f2>=f3>=f4, f4<= floor(A/4)
      // Once f4 is fix, then f3<=floor((A-f4)/3), etc.
      for(int d=0; d<=int(A/4);++d)
        for(int c=d; c<=int((A-d)/3); ++c)
          for(int b=c; b<=int((A-c-d)/2); ++b)
            {
              int a=A-b-c-d;        
              tableaux.push_back({a,b,c,d});
            }
      return tableaux;
    }


  class U4Subspace;
  class U4State;
  class U4BosonSpace;
  // Class for U(4) space, which is 
  // u4::U4Space (N particles)
  // -> u4::U4Subspace [f]
  //  -> STSubspace [S,T] beta
  
  // Subspace constructor
  class U4Subspace
      : public basis::
            BaseDegenerateSubspace<U4Subspace, std::tuple<u4::U4>, U4State, std::tuple<HalfInt,HalfInt>>
  {
   public:
    U4Subspace() = default;

    U4Subspace(
        const u4::U4& f
      )
        : BaseDegenerateSubspace{{f}}
    {
      auto weights = u4::GenerateU4toSU2SU2Weights(f);
      for(const auto& [S,T,multiplicity] : weights)
      {
       PushStateLabels({S,T}, multiplicity);
      }
    }

    u4::U4 f() const { return std::get<0>(labels()); }
    unsigned int Num_particles() const {return f().N();}
    std::string DebugStr(const std::string offset="") const;
    std::string LabelStr() const {return f().Str();}

   private:
  };

  // State constructor
  class U4State
      : public basis::BaseDegenerateState<U4Subspace>
  {
   public:
    // pass-through constructors

    U4State(const SubspaceType& subspace, std::size_t index)
    // Construct state by index.
        : basis::BaseDegenerateState<U4Subspace>(subspace, index)
    {}

    U4State(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
    // Construct state by reverse lookup on labels.
        : basis::BaseDegenerateState<U4Subspace>(subspace, state_labels)
    {}

    // pass-through accessors for subspace labels
    HalfInt S() const { return std::get<0>(labels()); }
    HalfInt T() const { return std::get<1>(labels()); }
    unsigned int beta_max() const { return subspace().GetStateDegeneracy(index()); }

    // private:
  };

  // Space constructor
  class U4Space
      : public basis::BaseSpace<U4Space, U4Subspace>
  {
    public:
      U4Space() = default;
      U4Space( const unsigned int N_particles_max)
      : BaseSpace{}
      {
        for(unsigned int N_particles=0; N_particles<=N_particles_max; ++N_particles)
        {
          auto tableaux = u4::GenerateU4Tableaux(N_particles);
          for (const auto& f : tableaux)
            {
              PushSubspace(U4Subspace(f));
            }
        }
      }
      
      std::string DebugStr(const std::string offset="") const;

    private:

  };

}


namespace un
{

  //Get N of U(Omega)->U(N)xU(4)
  inline unsigned int OmegaUN(unsigned int N) {
    return int((N+1)*(N+2)/2);
  }



  // Forward declarations
  class UNU4Subspace;
  class U3State;

  // UNU4subspace 
  //   For each U4 in U4Space, if has valid Sn_conjugate for U(N), add subspace labeled by
  //        u4_subspace_index and has pointer to U4Subspace. 
  //        -> U3State


  // Subspace constructor
  class UNU4Subspace
      : public basis::
            BaseDegenerateSubspace<UNU4Subspace, std::tuple<unsigned int,std::shared_ptr<const u4::U4Subspace>>, U3State, std::tuple<u3::U3>>
  {
   public:
    UNU4Subspace() = default;

    UNU4Subspace(
        unsigned int Nshell,  // corresponding to U(N) 
        const std::shared_ptr<const u4::U4Subspace>& u4subspace_ptr //pointer to corresponding U(4) subspace
      );

    u4::U4Subspace getU4Subspace() const {
      // add check for null pointer before dereferencing
      return (*std::get<1>(labels()));
    }
    u4::U4 f() const { return getU4Subspace().f();}

    // u4::U4 f() const { return (*std::get<1>(labels())).f(); }
    unsigned int num_particles() const {return f().N();}
    unsigned int Nshell() const {return std::get<0>(labels());}
    unsigned int N() const {return un::OmegaUN(Nshell());}


    std::string DebugStr(const std::string offset="") const;
    std::string LabelStr() const {return f().Str();}

   private:
  };

  // State constructor
  class U3State
      : public basis::BaseDegenerateState<UNU4Subspace>
  {
   public:
    // pass-through constructors

    U3State(const SubspaceType& subspace, std::size_t index)
    // Construct state by index.
        : basis::BaseDegenerateState<UNU4Subspace>(subspace, index)
    {}

    U3State(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
    // Construct state by reverse lookup on labels.
        : basis::BaseDegenerateState<UNU4Subspace>(subspace, state_labels)
    {}

    // pass-through accessors for subspace labels
    u3::U3 w() const { return std::get<0>(labels()); }
    unsigned int alpha_max() const { return subspace().GetStateDegeneracy(index()); }

    // private:
  };




  // Space constructor
  class UNU4Space
      : public basis::
            BaseSpace<UNU4Space,UNU4Subspace, std::tuple<unsigned int,std::shared_ptr<const u4::U4Subspace>> >
  {
   public:
    UNU4Space() = default;

    UNU4Space(
        unsigned int A, 
        unsigned int Nshell_max
      ): Nshell_max_(Nshell_max)
    {
      
      u4space_ = u4::U4Space(A);
      
      for(unsigned int Nshell=0; Nshell<=Nshell_max; ++Nshell)
      {
        unsigned int N = un::OmegaUN(Nshell); // N of U(N) corresponding to the shell.
        for(const auto& subspace : u4space_)
        {
          if(not subspace.f().HasSnConjugate(N)) continue;
            PushSubspace(UNU4Subspace(Nshell,std::make_shared<u4::U4Subspace>(subspace)));
        }
      }
    }

    inline const u4::U4Space& u4space() const {return u4space_;}

    // std::string DebugStr(const std::string offset="") const;
    // std::string LabelStr() const {return f().Str();}

   private:

    u4::U4Space u4space_;
    unsigned int Nshell_max_;

  };



  using U3U4ConfigurationType = std::vector<std::size_t>;
  using U3U4ConfigurationTableType = std::vector<U3U4ConfigurationType>;

  using U3U4IrrepTableType = std::unordered_map<
        std::tuple<u3::U3,u4::U4>,
        unsigned int,
        boost::hash<std::tuple<u3::U3,u4::U4>>
        >;
  using NexU3U4IrrepTableType = std::map<
      unsigned int, 
      U3U4IrrepTableType
      >;


  U3U4IrrepTableType GetU3U4ForConfigurations(
    const shell::Configuration& configuration,
    const un::UNU4Space& unu4space,
    // const unsigned int Nex,
    const unsigned int Nshell_max
    );


  struct U3U4Configurations
  {
    // Members 
    un::UNU4Space unu4space_;
    NexU3U4IrrepTableType  master_u3u4_irreps_map_;
    // std::vector<std::vector<un::UNU4Subspace>> unu4subspace_table_; //By Nshell, by subspace

    // std::vector<std::vector<un::wSTshell>> wSTshell_list_by_shell_;
    // wSTConfigurationTableType u3u4_configurations_table_;
    // std::map<unsigned int,std::pair<std::size_t,std::size_t>> Nex_limits_;

    

    // Null constructor
    U3U4Configurations():unu4space_() {}

    // Constructor

    /// Not going to work.  Need to branch to U3State first....
    U3U4Configurations(unsigned int A, const unsigned int Nshell_max, unsigned int Nmax)
    {
      // Generate and save U4Space for given number of particles A
      unu4space_ = un::UNU4Space(A,Nshell_max);

      // Distribute particles over Nshell_max+1 shells restricted to Nex<=Nmax
      std::vector<shell::ShellConfigurations> configurations_by_Nex = shell::generate_configurations(A,Nmax,Nshell_max);
      
      // {Nex,{<w,f>:mult}}
      
      
      for(unsigned int Nex=Nmax%2; Nex<=Nmax; Nex+=2)
      {
        
        std::vector<shell::Configuration> Nex_configurations(configurations_by_Nex[Nex].begin(),configurations_by_Nex[Nex].end());
        #pragma omp for schedule(dynamic) nowait
        for(auto i = 0; i<Nex_configurations.size(); ++i)
        // for(const shell::Configuration& configuration : configurations_by_Nex[Nex])
        {
          const shell::Configuration& configuration = Nex_configurations[i];
          auto u3u4_irreps_map = GetU3U4ForConfigurations(
              configuration,
              unu4space_,
              // Nex,
              Nshell_max
              );

          #pragma omp critical
          {
          for(const auto& [irrep,multiplicity] : u3u4_irreps_map)
            master_u3u4_irreps_map_[Nex][irrep]+=multiplicity;
          }
            
        }//end configuration
      }//Nex
    }//Constructor     

    inline const NexU3U4IrrepTableType & u3_u4_irreps(){return master_u3u4_irreps_map_;}

  };




}// end un namespace




#endif