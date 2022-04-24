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

namespace u3shell
{

namespace spatial::onecoord
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
  class OperatorParitySpace;
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
    std::string DebugStr(const std::string& indent="") const;

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
        // Casting shouldn't be necessary, but I'm paranoid (aem).
        return static_cast<unsigned int>(static_cast<int>(Nbar())+subspace().N0());
      }

      std::string LabelStr() const
        {return fmt::format("[{} {}]",Nbar(),Nbarp());}

      std::string DebugStr(const std::string& indent="") const;
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
          if(subspace.size()>0 && kappa0_max>0)
            PushSubspace(std::move(subspace),kappa0_max);
        }

    }

    // accessors
    unsigned int L0() const {return std::get<0>(labels());}
    unsigned int kappa0_max(std::size_t i) const {return GetSubspaceDegeneracy(i);}
    std::string LabelStr() const {return fmt::format("{}",L0());}
    std::string DebugStr(const std::string& indent="") const;
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
      // if no S0 values given, then S0max is default value
      int S0max = 2;
      const auto& allowed_S0_values = operator_parameters.Allowed_S0_values;
      if(allowed_S0_values.size()>0)
      // allowed_S0_values are in ordered set with largest value being the right-most
        S0max= *(allowed_S0_values.rbegin());

      // Extract allowed values for L0 from operator_parameters.
      std::set<unsigned int> allowed_L0_values;
      if(operator_parameters.Allowed_L0_values.size()>0)
         allowed_L0_values = operator_parameters.Allowed_L0_values;
      else
        {
          unsigned int L0min=static_cast<int>(std::max(0,int(operator_parameters.J0)-S0max));

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
    std::string LabelStr() const {return fmt::format("{}",N0());}
    std::string DebugStr(const std::string& indent="") const;

    // private:

  };


  //////////////////////////////////////////////////////////////////////////////////////
  /// Space
  ///   OperatorSpace
  ///     -> OperatorParitySpace [parity_bar]
  ///       -> OperatorN0Space [N0]
  ///         -> OperatorL0Space [L0]
  ///           -> OperatorSU3Subspace (kappa0)[omega0]
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
  class OperatorParitySpace
  : public basis::BaseSpace<OperatorParitySpace, OperatorN0Space,std::tuple<uint8_t>>
  {

   public:
    OperatorParitySpace() = default;

    OperatorParitySpace(
      const uint8_t parity_bar,
      const u3shell::relative::OperatorParameters& operator_parameters
    ) : BaseSpace{parity_bar}
    {
      // For parity conserving Nbar and Nbarp
      // If Nmax even,
      //    If parity_bar=1, then max(Nbarp-Nbar)=Nbar_max-2
      //    If parity_bar=0, then max(Nbarp-Nbar)=Nbar_max
      // If Nmax odd
      //    then max(Nbarp-Nbar) = Nbar_max-1
      int N0_max=operator_parameters.Nbar_max-2*parity_bar+(operator_parameters.Nbar_max%2);
      for(int N0 = -N0_max; N0<=N0_max; N0+=2)
        {
          auto subspace = OperatorN0Space(parity_bar,N0,operator_parameters);
          if(subspace.size()>0)
            PushSubspace(std::move(subspace));
        }

    }

    uint8_t parity_bar() const {return std::get<0>(labels());}

    std::string DebugStr(const std::string& indent) const;
    std::string LabelStr() const {return fmt::format("{}",parity_bar());}
    // private:

  };

  class OperatorSpace
  : public basis::BaseSpace<OperatorSpace,OperatorParitySpace>
  {

   public:
    OperatorSpace() = default;

    OperatorSpace(
      const u3shell::relative::OperatorParameters& operator_parameters
    ) : BaseSpace{}
    {
      for(uint8_t parity_bar : {0,1})
        {
          auto subspace = OperatorParitySpace(parity_bar,operator_parameters);
          if(subspace.size()>0)
            PushSubspace(std::move(subspace));
        }

    }

    std::string DebugStr() const;

    // private:

  };

  }  // namespace spatial::onecoord
  }  // namespace u3shell

#endif
