/****************************************************************
  recurrence_indexing_relative_operator.h

  Relative operator enumeration for spncci recurrence.
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/23/22 (aem): Created.
****************************************************************/

#ifndef RECURRENCE_OPERATOR_H_
#define RECURRENCE_OPERATOR_H_

#include <unordered_set>
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "sp3rlib/u3.h"

namespace relative
{
  struct OperatorParameters
  // Provides information about the relative two-body operator
  // Includes SU(3), angular momentum, spin, and isospin
  // quantum numbers for selection rules.
  {
    OperatorParameters(
      const int N1v_,
      const int Nmax_,
      const unsigned int J0_,
      const std::unordered_set<u3::U3>& Allowed_w0_values_,
      const std::set<unsigned int>& Allowed_L0_values_,
      const std::set<unsigned int>& Allowed_S0_values_,
      const std::set<unsigned int>& Allowed_T0_values_
    )
        :
        Nbar_max{Nmax_+2*N1v_},
        J0{J0_},
        Allowed_w0_values{Allowed_w0_values_},
        Allowed_L0_values{Allowed_L0_values_},
        Allowed_S0_values{Allowed_S0_values_},
        Allowed_T0_values{Allowed_T0_values_}

    {}

    const int Nbar_max;
    const unsigned int J0;
    const std::set<unsigned int> Allowed_L0_values;
    const std::set<unsigned int> Allowed_S0_values;
    const std::set<unsigned int> Allowed_T0_values;
    const std::unordered_set<u3::U3> Allowed_w0_values;
  };

namespace spatial
{
  class OperatorSpace;
  class OperatorU3Subspace;
  class OperatorState;

  ////////////////////////////////////////////////////////////////
  // x0 subspaces with degeneracy given by (L0,kappa0) pairs
  ////////////////////////////////////////////////////////////////
  class OperatorU3Subspace
    : public basis::BaseSubspace<
          OperatorU3Subspace,
          std::tuple<uint8_t, unsigned int, u3::SU3>,
          OperatorState,
          std::tuple<unsigned int, unsigned int>
        >
  {
  public:

    // constructors
    OperatorU3Subspace() = default;

    OperatorU3Subspace(
      const  uint8_t parity_bar,
      const unsigned int N0,
      const u3::SU3& x0,
      const std::vector<std::tuple<unsigned int,unsigned int>>& Nbar_pairs
    )
    : BaseSubspace{{parity_bar,N0,x0}}
    {
      for(const auto& Nbar_pair : Nbar_pairs){PushStateLabels(Nbar_pair);}
    }

    // accessors
    const uint8_t parity_bar() const {return std::get<0>(labels());}
    const unsigned int N0() const {return std::get<1>(labels());}
    const u3::SU3 x0() const {return std::get<2>(labels());}

    std::string LabelStr() const
      {return fmt::format("[{}{} {}]",N0(),x0(),parity_bar());}

    // private:
  };

  ////////////////////////////////////////////////////////////////
  // States given by Nbar,Nbar' pairs
  ////////////////////////////////////////////////////////////////
  class OperatorState
      : public basis::BaseState<OperatorU3Subspace>
    {
     public:

      OperatorState() = default;
      // pass-through constructors

      OperatorState(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
          : basis::BaseState<OperatorU3Subspace>(subspace, index)
      {}

      OperatorState(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
      // Construct state by reverse lookup on labels.
          : basis::BaseState<OperatorU3Subspace>(subspace, state_labels)
      {}

      unsigned int Nbar()  const {return std::get<0>(labels());}
      unsigned int Nbarp() const {return std::get<1>(labels());}

      std::string LabelStr() const
        {return fmt::format("[{} {}]",Nbar(),Nbarp());}

      // private:
    };

  ////////////////////////////////////////////////////////////////
  // Space
  // Takes as an argument an relative::OperatorParameters struct
  // which has methods
  // Nbar_max
  // J0
  // Allowed_L0_values
  // Allowed_S0_values
  // Allowed_T0_values
  // Allowed_x0_values
  //
  // If any the Allowed_{qn}_values sets are empty,
  // then all possible values of {qn} are considered.
  ////////////////////////////////////////////////////////////////
  // TODO: get rid of annoying mandatory SpaceLabelType for BaseDegenerateSpace
  // PATRICK ^^
  class OperatorSpace
  : public basis::BaseDegenerateSpace<OperatorSpace, OperatorU3Subspace, char>
  {

   public:
    OperatorSpace() = default;

    OperatorSpace(
      const relative::OperatorParameters& operator_parameters
    );

    // accessors
    //
    // Gets subspace offset for particular values of L0 and kappa0
    inline int GetSubspaceOffset(std::size_t i, int L0, int kappa0) const
      {
        int degeneracy_index = L0_offsets_[i][L0-L0min_]+kappa0;
        return BaseDegenerateSpace::GetSubspaceOffset(i, degeneracy_index);
      }

    // Extracts kappa0_max from stored L0 offsets.  Primarily for debugging.
    inline unsigned int Getkappa0max(std::size_t i, unsigned int L0) const
      {
        int l = L0-L0min_;
        if(l==4) {return BaseDegenerateSpace::GetSubspaceDegeneracy(i)-L0_offsets_[i][l];}
        else {return L0_offsets_[i][l+1]-L0_offsets_[i][l];}
      }

    std::string DebugStr() const;

    private:

      // Because S0 is at most 2, then the only possible values of L0 for
      // an operator that has good J0 are the 5 values |J0-2|, |J0-2|+1,...,J0+2
      // Alternatively, L0_min=max(0,J0-2), L0=L0_min+0, L0_min+1,...,L0_min+4.
      // The ith entry of L0_offsets gives the offset for L0+i.
      // The different between the offsets of subsequent L0 values (L0_min+i+1-L0_min+1)
      // corresponds to kappa0_max.
      std::vector<std::array<std::size_t,5>> L0_offsets_;
      unsigned int L0min_;
  };

}// spatial namespace

namespace spin
{

  class OperatorSpace;
  class OperatorS0Subspace; // S0
  class OperatorT0Subspace; // T0
  class OperatorState;      // Sbar,Sbar',Tbar,Tbar'


 class OperatorT0Subspace
  : public basis::BaseSubspace<
    OperatorT0Subspace,
    std::tuple<unsigned int>,
    OperatorState,
    std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>
  >
  {
   public:
    OperatorT0Subspace() = default;


    // TODO: Switch to storing unit tensor quantum numbers as bytes
    OperatorT0Subspace(const unsigned int T0, const unsigned int S0, const uint8_t parity_bar)
    : BaseSubspace{T0}
    {
      for(unsigned int Sbar : {0,1})
        for(unsigned int Sbarp : {0,1})
          for(unsigned int Tbar : {0,1})
            for(unsigned int Tbarp : {0,1})
              {
                if(!am::AllowedTriangle(Sbar,Sbarp,S0)) continue;
                if(!am::AllowedTriangle(Tbar,Tbarp,T0)) continue;
                if((Sbar+Tbar+parity_bar)%2!=1) continue;
                if((Sbarp+Tbarp+parity_bar)%2!=1) continue;

                PushStateLabels({Sbar,Sbarp,Tbar,Tbarp});
              }
    }

    //accessor
    const unsigned int T0() const {return std::get<0>(labels());}
    // // diagnostic output
    std::string DebugStr() const;

  // private:
  };


  ////////////////////////////////////////////////////////////////
  //
  ////////////////////////////////////////////////////////////////
  class OperatorState
      : public basis::BaseState<OperatorT0Subspace>
    {
     public:

      OperatorState() = default;
      // pass-through constructors

      OperatorState(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
          : basis::BaseState<OperatorT0Subspace>(subspace, index)
      {}

      OperatorState(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
      // Construct state by reverse lookup on labels.
          : basis::BaseState<OperatorT0Subspace>(subspace, state_labels)
      {}

      unsigned int Sbar() const {return std::get<0>(labels());}
      unsigned int Sbarp()const {return std::get<1>(labels());}
      unsigned int Tbar() const {return std::get<2>(labels());}
      unsigned int Tbarp()const {return std::get<3>(labels());}

      std::string LabelStr() const;

      // private:
    };


 class OperatorS0Subspace
  : public basis::BaseSpace<
    OperatorS0Subspace,
    OperatorT0Subspace,
    std::tuple<unsigned int,uint8_t>
  >
  {
   public:
    OperatorS0Subspace() = default;

    OperatorS0Subspace(
      const unsigned int S0,
      const uint8_t parity_bar,
      const OperatorParameters& operator_parameters
    )
    : BaseSpace{{S0,parity_bar}}
    {

      auto Allowed_T0_values = operator_parameters.Allowed_T0_values.size()?
        operator_parameters.Allowed_T0_values : std::set<unsigned int>{0,1,2};

      for(const auto T0 : Allowed_T0_values)
        PushSubspace(OperatorT0Subspace(T0,S0,parity_bar));

    }

    // accessors
    const unsigned int S0() const {return std::get<0>(labels());}
    const uint8_t parity_bar() const {return std::get<1>(labels());}

    // // diagnostic output
    std::string LabelStr() const {return fmt::format("[{} {}]",S0(),parity_bar());}
    std::string DebugStr() const;

  // private:
  };


 class OperatorSpace
  : public basis::BaseSpace<OperatorSpace, OperatorS0Subspace>
  {
   public:
    OperatorSpace() = default;

    OperatorSpace(
      const OperatorParameters& operator_parameters
    )
    {
      auto Allowed_S0_values = operator_parameters.Allowed_S0_values.size()?
          operator_parameters.Allowed_S0_values : std::set<unsigned int>{0,1,2};


      for(uint8_t parity_bar : {0,1})
        for(unsigned int S0 : Allowed_S0_values)
          {

            auto subspace = OperatorS0Subspace(S0,parity_bar,operator_parameters);
            if(subspace.dimension()==0)
              continue;

            PushSubspace(std::move(subspace));
          }
    }

    // // diagnostic output
    std::string DebugStr() const;

  // private:
  };

}// spin namespace



using OperatorSector
  = basis::BaseDegenerateSector<
      relative::spatial::OperatorU3Subspace,
      relative::spin::OperatorS0Subspace
      >;


class OperatorSectors
  : public basis::BaseSectors<
      relative::spatial::OperatorSpace,
      relative::spin::OperatorSpace,
      relative::OperatorSector,
      false
    >
{
  public:

  // Default constructor
  OperatorSectors() = default;

  // Constructor
  OperatorSectors(
    const relative::spatial::OperatorSpace& spatial_operator_space,
    const relative::spin::OperatorSpace& spin_operator_space
  );
};

}//relative namespace



#endif
