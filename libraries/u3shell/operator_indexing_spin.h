/****************************************************************
  operator_indexing_spin.h

  Operator enumeration for one-body and two-body spin operators
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/23/22 (aem): Created.
****************************************************************/

#ifndef OPERATOR_INDEXING_SPIN_H_
#define OPERATOR_INDEXING_SPIN_H_

#include <unordered_set>
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "u3shell/operator_parameters.h"
#include "mcutils/bit_tuple.h"

namespace u3shell::spin::twobody
{
////////////////////////////////////////////////////////////////
// labels for good-ST unit tensors
////////////////////////////////////////////////////////////////
class OperatorLabelsST;
}
////////////////////////////////////////////////////////////////
// Define standard hash
////////////////////////////////////////////////////////////////
namespace std
{
template<> struct hash<u3shell::spin::twobody::OperatorLabelsST>;
}

namespace u3shell::spin::twobody
{

class OperatorLabelsST
{
 private:
  constexpr explicit OperatorLabelsST(mcutils::bit_tuple<uint8_t,2,2,1,1,1,1> bitrep)
      : bitrep_{bitrep}{}


  friend struct std::hash<OperatorLabelsST>;

public:
  constexpr OperatorLabelsST(int S0, int T0, int Sbar, int Sbarp, int Tbar, int Tbarp)
      : bitrep_{S0,T0,Sbar,Sbarp,Tbar,Tbarp}{}

  static std::array<OperatorLabelsST, 36> const ALLOWED_LABELS;

  inline auto bitrep() const {return bitrep_;}

  inline bool Allowed() const
  {
    if (Sbar() < 0 || Sbar() > 1)
      return false;
    if (Sbarp() < 0 || Sbarp() > 1)
      return false;
    if (!am::AllowedTriangle(Sbar(), S0(), Sbarp()))
      return false;
    if (Tbar() < 0 || Tbar() > 1)
      return false;
    if (Tbarp() < 0 || Tbarp() > 1)
      return false;
    if (!am::AllowedTriangle(Tbar(), T0(), Tbarp()))
      return false;
    return true;
  }

  inline unsigned int S0() const
  {return static_cast<unsigned int>(mcutils::get<0>(bitrep_));}

  inline unsigned int T0() const
  {return static_cast<unsigned int>(mcutils::get<1>(bitrep_));}

  inline unsigned int Sbar() const
  {return static_cast<unsigned int>(mcutils::get<2>(bitrep_));}

  inline unsigned int Sbarp() const
  {return static_cast<unsigned int>(mcutils::get<3>(bitrep_));}

  inline unsigned int Tbar() const
  {return static_cast<unsigned int>(mcutils::get<4>(bitrep_));}

  inline unsigned int Tbarp() const
  {return static_cast<unsigned int>(mcutils::get<5>(bitrep_));}

  inline std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> state_labels() const
  {
    return {Sbar(),Sbarp(),Tbar(),Tbarp()};
  }

  inline std::tuple<unsigned int,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int>
  labels() const
  {
    return {S0(),T0(),Sbar(),Sbarp(),Tbar(),Tbarp()};
  }


  inline uint8_t exchange_symm_bar() const { return (Sbar() + Tbar()) % 2; }
  inline uint8_t exchange_symm_barp() const { return (Sbarp() + Tbarp()) % 2; }

  inline bool operator==(const OperatorLabelsST& rhs) const
  {
    return bitrep() == rhs.bitrep();
  }

  std::string LabelStr() const
  {
    return fmt::format("U[{} {}]({} {}, {} {})",S0(),T0(),Sbarp(),Tbarp(),Sbar(),Tbar());
  }

 private:
  mcutils::bit_tuple<uint8_t,2,2,1,1,1,1> bitrep_;
};

  // 36 possible two-body spin and isospin labels
  inline std::array<OperatorLabelsST, 36> const OperatorLabelsST::ALLOWED_LABELS{{
    // (S0, T0, Sbar, Sbarp, Tbar, Tbarp)
    OperatorLabelsST{0, 0, 0, 0, 0, 0},
    OperatorLabelsST{0, 0, 0, 0, 1, 1},
    OperatorLabelsST{0, 0, 1, 1, 0, 0},
    OperatorLabelsST{0, 0, 1, 1, 1, 1},
    OperatorLabelsST{0, 1, 0, 0, 0, 1},
    OperatorLabelsST{0, 1, 0, 0, 1, 0},
    OperatorLabelsST{0, 1, 0, 0, 1, 1},
    OperatorLabelsST{0, 1, 1, 1, 0, 1},
    OperatorLabelsST{0, 1, 1, 1, 1, 0},
    OperatorLabelsST{0, 1, 1, 1, 1, 1},
    OperatorLabelsST{0, 2, 0, 0, 1, 1},
    OperatorLabelsST{0, 2, 1, 1, 1, 1},
    OperatorLabelsST{1, 0, 0, 1, 0, 0},
    OperatorLabelsST{1, 0, 0, 1, 1, 1},
    OperatorLabelsST{1, 0, 1, 0, 0, 0},
    OperatorLabelsST{1, 0, 1, 0, 1, 1},
    OperatorLabelsST{1, 0, 1, 1, 0, 0},
    OperatorLabelsST{1, 0, 1, 1, 1, 1},
    OperatorLabelsST{1, 1, 0, 1, 0, 1},
    OperatorLabelsST{1, 1, 0, 1, 1, 0},
    OperatorLabelsST{1, 1, 0, 1, 1, 1},
    OperatorLabelsST{1, 1, 1, 0, 0, 1},
    OperatorLabelsST{1, 1, 1, 0, 1, 0},
    OperatorLabelsST{1, 1, 1, 0, 1, 1},
    OperatorLabelsST{1, 1, 1, 1, 0, 1},
    OperatorLabelsST{1, 1, 1, 1, 1, 0},
    OperatorLabelsST{1, 1, 1, 1, 1, 1},
    OperatorLabelsST{1, 2, 0, 1, 1, 1},
    OperatorLabelsST{1, 2, 1, 0, 1, 1},
    OperatorLabelsST{1, 2, 1, 1, 1, 1},
    OperatorLabelsST{2, 0, 1, 1, 0, 0},
    OperatorLabelsST{2, 0, 1, 1, 1, 1},
    OperatorLabelsST{2, 1, 1, 1, 0, 1},
    OperatorLabelsST{2, 1, 1, 1, 1, 0},
    OperatorLabelsST{2, 1, 1, 1, 1, 1},
    OperatorLabelsST{2, 2, 1, 1, 1, 1}
  }};


  inline std::size_t GetSpinOperatorOffset(const unsigned int S0)
      {
        for(std::size_t index=0; index<36; index++)
          {
            if(OperatorLabelsST::ALLOWED_LABELS[index].S0()==S0)
              return index;
          }
        fmt::print("S0 = {} is not an allowed two-body spin");
        exit(EXIT_FAILURE);
      }

  inline std::size_t GetSpinOperatorOffset(const unsigned int S0, const unsigned int T0)
      {
        for(std::size_t index=0; index<36; index++)
          {
            if(OperatorLabelsST::ALLOWED_LABELS[index].S0()==S0)
              if(OperatorLabelsST::ALLOWED_LABELS[index].T0()==T0)
                return index;
          }

        exit(EXIT_FAILURE);
      }
  ////////////////////////////////////////////////////////////////
  class OperatorState;
  class OperatorSubspace;
  class OperatorSpace;
  ////////////////////////////////////////////////////////////////
  // OperatorSubspace [exchange_syn_bar]
  //  -> OperatorStates [i]
  // States within subspace are index into
  // OperatorLabelsST::ALLOWED_LABELS
  // corresponding to S0 T0 Sbar Sbarp Tbar Tbarp
  ////////////////////////////////////////////////////////////////
  class OperatorSubspace
      : public basis::BaseSubspace<
      OperatorSubspace,
      std::tuple<uint8_t,uint8_t>,
      OperatorState,
      std::tuple<uint8_t>
    >
  {
     public:
      // constructors
      OperatorSubspace() = default;

      OperatorSubspace(
        const uint8_t exchange_symm_bar,
        const uint8_t S0,
        const u3shell::relative::OperatorParameters& operator_parameters
        ): BaseSubspace{{exchange_symm_bar,S0}}
      {

        auto allowed_T0_values =
            (operator_parameters.Allowed_T0_values.size() == 0)?
            std::set<uint8_t>{0, 1, 2}
            :operator_parameters.Allowed_T0_values;

        for (int i = 0; i < 36; ++i)
        {
          const auto& labels = OperatorLabelsST::ALLOWED_LABELS[i];
          bool allowed_state = S0==labels.S0();
          allowed_state &= allowed_T0_values.count(labels.T0());
          allowed_state &= (labels.Sbar()+labels.Tbar())%2 == exchange_symm_bar;
          allowed_state &= (labels.Sbarp()+labels.Tbarp())%2 == exchange_symm_bar;

          if(allowed_state)
            PushStateLabels({i});
        }
      }

      // accessors
      uint8_t exchange_symm_bar() const { return std::get<0>(labels()); }
      uint8_t S0() const { return std::get<1>(labels()); }
      std::string DebugStr(const std::string& indent = "") const;
      std::string LabelStr() const {return fmt::format("[{} {}]",exchange_symm_bar(),S0());}
     // private:

  };


  ////////////////////////////////////////////////////////////////
  // States: index into OperatorLabelsST::ALLOWED_LABELS
  ////////////////////////////////////////////////////////////////
  class OperatorState
      : public basis::BaseState<OperatorSubspace>
    {
     public:

      OperatorState() = default;
      OperatorState(const SubspaceType& subspace, std::size_t index)
      // Construct state by index.
          : basis::BaseState<OperatorSubspace>(subspace, index)
      {}

      OperatorState(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
      // Construct state by reverse lookup on labels.
          : basis::BaseState<OperatorSubspace>(subspace, state_labels)
      {}

      inline int index()  const {return static_cast<int>(std::get<0>(BaseState::labels()));}
      inline std::tuple<unsigned int,unsigned int,unsigned int,unsigned int,unsigned int,unsigned int>
        labels() const {return OperatorLabelsST::ALLOWED_LABELS[index()].labels();}

      std::string LabelStr() const
        {return OperatorLabelsST::ALLOWED_LABELS[index()].LabelStr();}
    };


class OperatorSpace
  : public basis::BaseSpace<
    OperatorSpace, OperatorSubspace
    >
  {

   public:
    OperatorSpace() = default;

    OperatorSpace(
      const u3shell::relative::OperatorParameters& operator_parameters
    ) : BaseSpace{}
    {

      auto allowed_S0_values =
          (operator_parameters.Allowed_S0_values.size() == 0)?
          std::set<uint8_t>{0, 1, 2}
          :operator_parameters.Allowed_S0_values;


      for(uint8_t exchange_symm_bar : {0,1})
      {
        for (uint8_t S0 : allowed_S0_values)
        {
          auto subspace =
              OperatorSubspace(exchange_symm_bar, S0, operator_parameters);

          if(subspace.size()>0)
            PushSubspace(std::move(subspace));
        }
      }

    }

    std::string DebugStr() const;

    // private:

  };

}  // namespace u3shell::spin::twobody


namespace std
{
template<> struct hash<u3shell::spin::twobody::OperatorLabelsST>
{
  inline std::size_t operator()(const u3shell::spin::twobody::OperatorLabelsST& h) const noexcept
  {
    return std::hash<mcutils::bit_tuple<uint8_t,2,2,1,1,1,1>>()(h.bitrep());
  }
};
}  // namespace std

namespace u3shell::spin::twobody
{
inline std::size_t hash_value(const OperatorLabelsST& l)
{
  return std::hash<OperatorLabelsST>{}(l);
}
}


#endif
