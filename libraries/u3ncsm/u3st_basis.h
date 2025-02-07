/****************************************************************
  u4.h

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT

  + 2/6/25 (aem): Created.
****************************************************************/
#ifndef U3ST_BASIS_H_
#define U3ST_BASIS_H_

#include <unordered_map>
#include <map>
#include <functional>
#include "basis/basis.h"
#include "basis/degenerate.h" 
#include "sp3rlib/u3.h"
#include "su4lib/u4.h"
#include "su4lib/u4_branching.h"
#include "am/halfint.h"
#include "utilities/nuclide.h"
#include "u3ncsm/configurations.h"
#include "u3ncsm/u3u4_basis.h"

namespace un
{
  // Structure containing pointer to corresponding unu4subspace
  // w_index gives state index within unu4subspace
  // st_index gives state index within u4subspace associated with unu4subspace
  struct wSTshell
  {
    // Members 
    std::shared_ptr<const un::UNU4Subspace> unu4_ptr;
    std::size_t w_index;
    std::size_t st_index;

    wSTshell():unu4_ptr(nullptr),w_index(basis::kNone),st_index(basis::kNone){} //null constructor for use in, e.g., map

    // Constructor 
    wSTshell(
      const std::shared_ptr<const un::UNU4Subspace>& unu4_ptr_, 
      const std::size_t w_index_, 
      const std::size_t st_index_
      ): unu4_ptr(unu4_ptr_), w_index(w_index_), st_index(st_index_){}

    wSTshell(
      const un::UNU4Subspace& unu4_subspace, 
      const std::size_t w_index_, 
      const std::size_t st_index_
      )
    {
      std::shared_ptr<un::UNU4Subspace> unu4_ptr_ = std::make_shared<un::UNU4Subspace>(unu4_subspace);
      wSTshell(unu4_ptr_, w_index_,st_index_);
    }

    // Accessor
    inline un::UNU4Subspace getUNU4Subspace() const {return *unu4_ptr;}
    inline u4::U4 f() const {return getUNU4Subspace().f();}
    inline HalfInt S() const {return getUNU4Subspace().getU4Subspace().GetState(st_index).S();}
    inline HalfInt T() const {return getUNU4Subspace().getU4Subspace().GetState(st_index).T();}
    inline unsigned int beta_max() const {return getUNU4Subspace().getU4Subspace().GetState(st_index).beta_max();}
    inline u3::U3 w() const {return getUNU4Subspace().GetState(w_index).w();}
    inline unsigned int alpha_max() const {return getUNU4Subspace().GetState(w_index).alpha_max();}
    inline unsigned int multiplicity() const {return beta_max()*alpha_max();}

    typedef std::tuple<u4::U4,u3:: U3,HalfInt,HalfInt> KeyType;
    inline KeyType Key() const
    {
      return {f(),w(),S(),T()};
    }


    inline friend bool operator == (const wSTshell& a, const wSTshell& b)
    {
      return a.Key() == b.Key();
    }

    inline friend bool operator < (const wSTshell& a, const wSTshell& b)
    {
      return a.Key() < b.Key();
    }
    inline friend std::size_t hash_value(const wSTshell& v)
    {
      boost::hash<wSTshell::KeyType> hasher;
      return hasher(v.Key());
    }


  };


  using wSTConfigurationType = std::vector<std::size_t>;
  using wSTConfigurationTableType = std::vector<wSTConfigurationType>;


  struct wSTConfigurations
  {
    // Members 
    u4::U4Space u4space_;
    std::vector<std::vector<un::UNU4Subspace>> unu4subspace_table_; //By Nshell, by subspace

    std::vector<std::vector<un::wSTshell>> wSTshell_list_by_shell_;
    wSTConfigurationTableType wSTshell_configurations_table_;
    std::map<unsigned int,std::pair<std::size_t,std::size_t>> Nex_limits_;

    

    // Null constructor
    wSTConfigurations():  
      wSTshell_list_by_shell_(), 
      wSTshell_configurations_table_(),
      u4space_(),
      unu4subspace_table_(),
      Nex_limits_()
    {}

    // Constructor
    wSTConfigurations(unsigned int A, const unsigned int Nshell_max, unsigned int Nmax);

    // accessors
    inline const std::vector<std::vector<wSTshell>>& wSTshell_list_by_shell() const
    {return wSTshell_list_by_shell_;}

    inline const wSTConfigurationTableType& configurations_table() const
      {return wSTshell_configurations_table_;}

    // inline std::vector<std::reference_wrapper<wSTConfigurationType>> 
    //   configurations_table(unsigned int  Nex) const
    //   {
    //     std::vector<std::reference_wrapper<wSTConfigurationType>> subvector;
    //     for(std::size_t i = Nex_limits_.at(Nex).first; i<Nex_limits_.at(Nex).second; ++i)
    //       subvector.push_back(configurations_table()[i]);
        
    //     return subvector;
    //   }

    inline const u4::U4Space& U4Space() const {return u4space_;}

    inline const std::vector<std::vector<un::UNU4Subspace>>& unu4subspace_table()const
      {return unu4subspace_table_;}


    inline std::pair<std::size_t,std::size_t>
      Nex_limits(unsigned int Nex) const {return Nex_limits_.at(Nex);}

    
    inline const un::wSTshell& getwSTshell(std::size_t config_index, unsigned int Nshell) const 
      { 
          return wSTshell_list_by_shell()[Nshell][configurations_table()[config_index][Nshell]];
      }


    inline std::vector<un::wSTshell> configuration_labels(const std::size_t index)
    {
      
      const auto& config = configurations_table()[index];
      std::vector<un::wSTshell> labels(config.size());
      for(auto i=0; i<config.size(); ++i)
      {
          labels[i] = wSTshell_list_by_shell()[i][config[i]];
        
      }
      
      return labels;
    }


  };



 // Forward declarations
  class U3STSpace; // space containing all wST
  class wSTSubSpace; // subspace labeled by [wST] containing all coupled configurations which results in total wST quantum numbers 
  class wSTConfigurationSubspace; //Subspace for each wSTshell_configuration, stored as vector of indices into wSTshell_lists_by_shell_
  class wSTIntermediateCouplings; //State // indices into tables corresponding to su3 coupling and s/t coupling. Degeneracy coming from w couplings


  using WSTConfigurationSubspaceLabelType = std::shared_ptr<std::vector<wSTshell>>;
  // Subspace constructor
  class wSTConfigurationSubspace
      : public basis::
            BaseDegenerateSubspace<
              wSTConfigurationSubspace,
              std::tuple<unsigned int>, 
              wSTIntermediateCouplings, 
              std::tuple<std::size_t,std::size_t,std::size_t>
            >
  {
   public:
    // wSTConfigurationSubspace() = default;

    // TODO
    wSTConfigurationSubspace(
        const unsigned int config_index,
        const std::vector<std::tuple<std::size_t,std::size_t,std::size_t,unsigned int>>& intermediate_couplings
      ): BaseDegenerateSubspace{{config_index}}
    {
      for(const auto& intermediate_coupling : intermediate_couplings)
      {
        const auto& [w_index,s_index,t_index,multiplicity] = intermediate_coupling;
         PushStateLabels({w_index,s_index,t_index}, multiplicity);
      }
    }

    // un::wSTshell wSTshell(const unsigned int Nshell) const {return (*(std::get<0>(labels())))[Nshell];}
    std::string DebugStr(const std::string offset="") const;
    // std::string LabelStr() const {return f().Str();}

   private:
    

  };



  // Subspace constructor
  class wSTSubspace
      : public basis::BaseDegenerateSpace<
              wSTSubspace,
              wSTConfigurationSubspace, 
              std::tuple<unsigned int,u3::SU3,HalfInt,HalfInt>
              // std::tuple<unsigned int> // index of wSTConfiguration in wst_configurations_.configurations_table()[Nex]
            >
  {
   
   public:
    
    wSTSubspace(
      const std::tuple<unsigned int,u3::SU3,HalfInt,HalfInt>& labels,
      const un::wSTConfigurations& wst_configurations,
      //key: index for configuration in table table
      //value: vector of <inter_w_index, inter_s_index,inter_t_index,degeneracy>
      const std::map<
          std::size_t, 
          std::vector<std::tuple<std::size_t,std::size_t,std::size_t,unsigned int>> 
          >& subspace_map
      ):BaseDegenerateSpace(labels)
    {
      
      //Config table is empty...why?
      const auto& config_table = wst_configurations.configurations_table();

      const auto& shell_table = wst_configurations.wSTshell_list_by_shell();
      for(const auto& [config_index, subspace_info] : subspace_map)
      {
        unsigned int multiplicity=1; 
        const auto& config = config_table[config_index];

        for (auto shell_index=0; shell_index<config.size(); ++ shell_index)
        {
          multiplicity *= shell_table[shell_index][config[shell_index]].multiplicity();
        }
        
        PushSubspace(wSTConfigurationSubspace(config_index,subspace_info),multiplicity);
        
      }

    }

    

    
    inline unsigned int Nex() const {return std::get<0>(labels());}
    inline u3::SU3 SU3() const {return std::get<1>(labels());}
    inline HalfInt S() const {return std::get<2>(labels());}
    inline HalfInt T() const {return std::get<3>(labels());}

    inline std::string DebugStr(const std::string offset="") const;
    std::string LabelStr() const {return fmt::format("{}{} {} {}",Nex(),SU3(),S(),T());}

   private:
  };



  // Space constructor
  class U3STSpace
      : public basis::BaseSpace<
      U3STSpace, 
      wSTSubspace,
      std::tuple<unsigned int,u3::SU3,HalfInt,HalfInt> // Nex(lambda,mu)ST
      > 
  {
    public:
      // wSTSpace() = default;

      U3STSpace( const unsigned int A, const unsigned int Nshell_max, const unsigned int Nmax,HalfInt Tz=-1);
      
      // U3STSpace(const nuclide::NuclideType& nuclide, const unsigned int Nshell_max, const unsigned int Nmax)
      // {
      //   U3STSpace(nuclide[0]+nuclide[1],Nshell_max,Nmax,HalfInt(abs(nuclide[0]-nuclide[1]),2));
      // }
      

    private:
      std::vector<std::vector<u3::U3>> intermediate_couplings_w_;
      std::vector<std::vector<HalfInt>> intermediate_couplings_st_;
      un::wSTConfigurations wst_configurations_;




  };

}

namespace std
{
  template<> struct hash<un::wSTshell>
  {
    inline std::size_t operator()(const un::wSTshell& h) const
    {
      return hash_value(h);
    }
  };

template <>
  struct hash<std::vector<un::wSTshell>> {
      inline size_t operator()(const std::vector<un::wSTshell>& vec) const {
          size_t result = 0;
          for (const auto& obj : vec) {
              result ^= std::hash<un::wSTshell>{}(obj); 
          }
          return result;
      }
  };
}

namespace fmt
{
  template<> struct formatter<un::wSTshell>
  {
    char presentation = 'g';

    template<typename ParseContext>
    FMT_CONSTEXPR auto parse(ParseContext& ctx) -> decltype(ctx.begin())
    {
      auto it = ctx.begin(), end = ctx.end();
      if (it != end && (*it == 'd' || *it == 'g' || *it == 'f'))
        presentation = *it++;

      // Check if reached the end of the range:
      if (it != end && *it != '}')
        throw format_error("invalid format");

      // Return an iterator past the end of the parsed range:
      return it;
    }


// return {f(),w(),S(),T()};
    template<typename FormatContext>
    FMT_CONSTEXPR auto format(const un::wSTshell& labels, FormatContext& ctx)
        -> decltype(ctx.out())
    {
      // if (presentation == 'f')
        return fmt::format_to(ctx.out(), "{} {} {} {}", labels.f(), labels.w(), labels.S(), labels.T());
      // else
      //   return fmt::format_to(ctx.out(), "{:g}{:d}", w.N(), w.SU3());
    }
  };
}  // namespace fmt
#endif