/****************************************************************
  operator_indexing_spatial.h

  Indexing for operators defined by matrix elements in a harmonic
  oscillator basis.
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/23/22 (aem): Created.
****************************************************************/

#ifndef OPERATOR_INDEXING_SPATIAL_H_
#define OPERATOR_INDEXING_SPATIAL_H_

#include <unordered_set>
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "sp3rlib/u3.h"
#include "u3shell/operator_parameters.h"

namespace spatial
{

namespace onecoord
{

  std::map<u3::SU3,std::vector<unsigned int>>
  GenerateSpatialOperators(
    const int N0,
    const uint8_t parity_bar,
    const u3shell::relative::OperatorParameters& operator_parameters
  );

  ////////////////////////////////////////////////////////////////
  // By N0, by L0, by kappa0, by rho0,omega0, by Nbar,Nbarp
  //
  // States
  ////////////////////////////////////////////////////////////////
  class OperatorNbarStates;
  class OperatorSU3Subspace;
  class OperatorL0Space;
  class OperatorN0Space;
  class OperatorSpace;

  ////////////////////////////////////////////////////////////////
  // x0 subspaces
  // States within subspace are Nbar,Nbarp pairs,
  // but containing spaces has fixed N0, so only Nbar is stored.
  ////////////////////////////////////////////////////////////////
  class OperatorU3Subspace
    : public basis::BaseSubspace<
          OperatorU3Subspace,
          std::tuple<u3::SU3>,
          OperatorNbarStates,
          std::tuple<unsigned int>
        >
  {
  public:

    // constructors
    OperatorU3Subspace() = default;

    OperatorU3Subspace(
      const int N0,
      const u3::SU3& x0,
      const std::vector<unsigned int>& Nbar_vector
    )
    : BaseSubspace{x0}, N0_{N0}
    {
      for(const auto& Nbar : Nbar_vector){PushStateLabels(Nbar);}
    }

    // accessors
    u3::SU3 x0() const {return std::get<0>(labels());}
    int N0() const {return N0_;}
    std::string LabelStr() const {return x0().Str();}

    private:
      int N0_;
  };


  ////////////////////////////////////////////////////////////////
  // States given by Nbar,Nbar' pairs
  ////////////////////////////////////////////////////////////////
  class OperatorNbarStates
      : public basis::BaseState<OperatorU3Subspace>
    {
     public:

      OperatorNbarStates() = default;
      // pass-through constructors

      OperatorNbarStates(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
          : basis::BaseState<OperatorU3Subspace>(subspace, index)
      {}

      OperatorNbarStates(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
      // Construct state by reverse lookup on labels.
          : basis::BaseState<OperatorU3Subspace>(subspace, state_labels)
      {}

      unsigned int Nbar()  const {return std::get<0>(labels());}

      //Reconstruct Nbarp from Nbar and N0;
      inline unsigned int Nbarp() const
      {
        return Nbar()+subspace().N0();
      }

      std::string LabelStr() const
        {return fmt::format("[{} {}]",Nbar(),Nbarp());}

      // private:
    };


  //////////////////////////////////////////////////////////////////////////////////////
  /// OperatorN0Space: operators grouped by N0. Contains OperatorU3Subspaces with
  /// degeneracy given by (L0,kappa0)
  /////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // L0 subspace: constains x0 subspaces with degeneracy kappa0_max
  ////////////////////////////////////////////////////////////////
  class OperatorL0Space
    : public basis::BaseDegenerateSpace<
          OperatorL0Space,
          OperatorU3Subspace,
          std::tuple<unsigned int>
        >
  {
  public:

    // constructors
    OperatorL0Space() = default;

    OperatorL0Space(
      const int N0, const unsigned int L0,
      const std::map<u3::SU3,std::vector<unsigned int>>& x0_Nbar_vector
    )
    : BaseDegenerateSpace{L0}
    {
      for(const auto&[x0,Nbar_vector] : x0_Nbar_vector)
        {
          unsigned int kappa0_max = u3::BranchingMultiplicitySO3(x0,L0);
          auto subspace = OperatorU3Subspace(N0,x0,Nbar_vector);
          if(subspace.size()>0)
            PushSubspace(std::move(subspace),kappa0_max);
        }

    }

    // accessors
    unsigned int L0() const {return std::get<0>(labels());}
    unsigned int kappa0_max(std::size_t i) const {return GetSubspaceDegeneracy(i);}
    std::string LabelStr() const {return fmt::format("{}",L0());}

    // private:

  };

  //////////////////////////////////////////////////////////////////////////////////////
  /// OperatorN0Space: operators grouped by N0. Contains OperatorL0Spaces
  /////////////////////////////////////////////////////////////////////////////////////
  class OperatorN0Space
  : public basis::BaseSpace<OperatorN0Space, OperatorL0Space,std::tuple<int>>
  {

   public:
    OperatorN0Space() = default;

    OperatorN0Space(
      const uint8_t parity_bar, const int N0,
      const u3shell::relative::OperatorParameters& operator_parameters
    ): BaseSpace{N0}
    {

      // Extract S0max from operator_parameters.
      int S0max = 2;
      const auto& allowed_S0_values = operator_parameters.Allowed_S0_values;
      // allowed_S0_values are in ordered set with largest value being the right-most
      if(allowed_S0_values.size()>0)
        S0max= *(allowed_S0_values.rbegin());

      // Extract allowed values for L0 from operator_parameters.
      std::set<unsigned int> allowed_L0_values;
      if(operator_parameters.Allowed_L0_values.size()>0)
         allowed_L0_values = operator_parameters.Allowed_L0_values;
      else
        {
          unsigned int L0min=std::max(0,int(operator_parameters.J0)-S0max);

          for(unsigned int L0=L0min; L0<=operator_parameters.J0+S0max; ++L0)
            {
              allowed_L0_values.insert(L0);
            }
        }

      std::map<u3::SU3,std::vector<unsigned int>> x0_Nbar_vector
        = GenerateSpatialOperators(N0,parity_bar,operator_parameters);

      // Iterate over L0 values and push L0Subspaces
      for(const unsigned int L0 : allowed_L0_values)
        {
          auto subspace = OperatorL0Space(N0,L0,x0_Nbar_vector);

          if(subspace.size()>0)
            PushSubspace(std::move(subspace));
        }
    }

    int N0() const {return std::get<0>(labels());}
    std::string DebugStr() const;

    // private:

  };


  //////////////////////////////////////////////////////////////////////////////////////
  /// Space
  ///   OperatorSpace [parity_bar]
  ///     -> OperatorN0Space [N0]
  ///       -> OperatorL0Space [L0]
  ///          -> OperatorSU3Subspace (kappa0)[omega0]
  ///              ->OperatorStates [Nbar]  Nbarp fixed by N0.
  ///
  /// Takes as an argument an u3shell::relative::OperatorParameters struct
  /// which has methods:
  ///   Nbar_max
  ///   J0
  ///   Allowed_L0_values
  ///   Allowed_S0_values
  ///   Allowed_T0_values
  ///   Allowed_w0_values
  ///
  /// If the Allowed_{qn}_values.size()==0, then all possible values of qn are considered.
  /// For e.g., ab initio Hamiltonin, operator_parameters would have
  ///
  ///    Nbar_max=Nmax+2*N1v
  ///    J0=0
  ///    Allowed_L0_values={} -> L0min=max(0,J0-2), L0=L0min,L0min+1,...,L0min+4;
  ///    Allowed_S0_values={} -> T0={0,1,2}
  ///    Allowed_T0_values={} -> T0={0,1,2}
  ///    Allowed_w0_values={} -> Given by GenerateSpatialOperators
  ///
  /////////////////////////////////////////////////////////////////////////////////////
  class OperatorU3Space
  : public basis::BaseSpace<OperatorU3Space, OperatorN0Space,std::tuple<uint8_t>>
  {

   public:
    OperatorU3Space() = default;

    OperatorU3Space(
      const uint8_t parity_bar,
      const u3shell::relative::OperatorParameters& operator_parameters
    ) : BaseSpace{parity_bar}
    {
      int N0_max=operator_parameters.Nbar_max-parity_bar;
      for(int N0 = -N0_max; N0<=N0_max; N0+=2)
        {
          auto subspace = OperatorN0Space(parity_bar,N0,operator_parameters);
          if(subspace.size()>0)
            PushSubspace(std::move(subspace));
        }

    }

    uint8_t parity_bar() const {return std::get<0>(labels());}

    std::string DebugStr() const;

    // private:

  };



  // [TODO: Rewrite with new U3Space indexing in mind]
  ////////////////////////////////////////////////////////////////
  // For upcoupling, define spatial indexing for relative
  ////////////////////////////////////////////////////////////////

  class OperatorSO3Space;
  class OperatorSO3Subspace;
  class OperatorSO3State;

  ////////////////////////////////////////////////////////////////
  // L0 subspaces with degeneracy given by (L0,kappa0) pairs
  ////////////////////////////////////////////////////////////////
  class OperatorSO3Subspace
    : public basis::BaseDegenerateSubspace<
          OperatorSO3Subspace,
          std::tuple<uint8_t, unsigned int>,
          OperatorSO3State,
          std::tuple<unsigned int, unsigned int>
        >
  {
  public:

    // constructors
    OperatorSO3Subspace() = default;

    OperatorSO3Subspace(
      const uint8_t parity_bar,
      const unsigned int Nmax,
      const unsigned int L0
    )
    : BaseDegenerateSubspace{{parity_bar,L0}}
    {
      for(int Lbar = parity_bar; Lbar<=Nmax; Lbar+=2)
        for(int Lbarp = std::abs(int(Lbar)-int(L0)); Lbarp<(Lbar+L0)&&Lbarp<=Nmax; Lbarp+=2)
          {
            int np_max = (Nmax-Lbarp-parity_bar)/2;
            int n_max = (Nmax-Lbar-parity_bar)/2;
            int degeneracy = np_max*n_max;
            PushStateLabels({Lbar,Lbarp},degeneracy);
            np_max_values_.push_back(np_max);
          }
    }

    // accessors
    uint8_t parity_bar() const {return std::get<0>(labels());}
    unsigned int L0() const {return std::get<1>(labels());}

    inline std::size_t GetStateOffset(const std::size_t i, const int Nbar, const int Nbarp) const
      {
        const auto& [Lbar,Lbarp] = GetStateLabels(i);
        int n = (Nbar-Lbar)/2;
        int np = (Nbarp-Lbarp)/2;
        int degeneracy_index = n*np_max_values_[i]+np;
        return BaseDegenerateSubspace::GetStateOffset(i, degeneracy_index);
      }

    std::string LabelStr() const
      {return fmt::format("[{} {}]",L0(),parity_bar());}

    private:
      std::vector<unsigned int> np_max_values_;
  };


  ////////////////////////////////////////////////////////////////
  // States given by Nbar,Nbar' pairs
  ////////////////////////////////////////////////////////////////
  class OperatorSO3State
      : public basis::BaseState<OperatorSO3Subspace>
    {
     public:

      OperatorSO3State() = default;
      // pass-through constructors

      OperatorSO3State(const SubspaceType& subspace, const std::size_t index)
      // Construct state by index.
          : basis::BaseState<OperatorSO3Subspace>(subspace, index)
      {}

      OperatorSO3State(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
      // Construct state by reverse lookup on labels.
          : basis::BaseState<OperatorSO3Subspace>(subspace, state_labels)
      {}

      unsigned int Lbar()  const {return std::get<0>(labels());}
      unsigned int Lbarp() const {return std::get<1>(labels());}

      std::string LabelStr() const
        {return fmt::format("[{} {}]",Lbar(),Lbarp());}

      // private:
    };

  ////////////////////////////////////////////////////////////////
  //  OperatorSO3Space
  //   -> OperatorSO3Subspace [L0] {Nbar,Nbar'}
  //       -> OperatorSO3State [Lbar,Lbarp]
  ////////////////////////////////////////////////////////////////
  // QUERRY PATRICK: Can we get rid of mandatory SpaceLabelType for BaseDegenerateSpace
  class OperatorSO3Space
  : public basis::BaseSpace<OperatorSO3Space, OperatorSO3Subspace>
  {

   public:
    OperatorSO3Space() = default;

    OperatorSO3Space(const unsigned int J0, const unsigned int Nbar_max)
    {
      int L0_min = std::max(0,int(J0)-2);
      for(int parity_bar : {0,1})
        for(unsigned int L0=L0_min; L0<=J0+2; L0++)
          {
            auto subspace = OperatorSO3Subspace(parity_bar,Nbar_max,L0);
            if(subspace.size()>0)
              PushSubspace(std::move(subspace));
          }
    }

    OperatorSO3Space(const u3shell::relative::OperatorParameters& operator_parameters)
    {
      const auto& Allowed_L0_values = operator_parameters.Allowed_L0_values;
      const unsigned int Nbar_max = operator_parameters.Nbar_max;
      for(int parity_bar : {0,1})
        for(const unsigned int L0 : Allowed_L0_values)
          {
            auto subspace = OperatorSO3Subspace(parity_bar,Nbar_max,L0);
            if(subspace.size()>0)
              PushSubspace(std::move(subspace));
          }
    }

    std::string DebugStr() const;
    // private:
  };


}// spatial namespace


}//relative namespace



#endif
