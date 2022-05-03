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

namespace u3shell::spin::twobody{

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


  std::size_t GetSpinOperatorOffset(const unsigned int S0)
    {
      for(std::size_t index=0; index<36; index++)
        {
          if(OperatorLabelsST::ALLOWED_LABELS[index].S0()==S0)
            return index;
        }

      exit(EXIT_FAILURE);
    }


  //[TODO: determine if we actually need this space.]
  ////////////////////////////////////////////////////////////////
  /// Spin operator space
  ///
  ///  spin::twobody::OperatorSpace
  ///  -> spin::twobody::OperatorS0Subspace [S0,parity_bar]
  // /    -> spin::twobody::OperatorT0Subspapce [T0] remove
  ///      -> spin::twobody::OperatorState [T0,[Sbar,Sbarp,Tbar,Tbarp]]
  ////////////////////////////////////////////////////////////////

  class OperatorSpace;
  class OperatorS0Subspace;
  class OperatorT0Subspace;
  class OperatorState;

 class OperatorT0Subspace
  : public basis::BaseSubspace<
    OperatorT0Subspace,
    std::tuple<uint8_t>,
    OperatorState,
    mcutils::bit_tuple<uint8_t,1,1,1,1>
  >
  {
   public:
    OperatorT0Subspace() = default;

    OperatorT0Subspace(const uint8_t& T0, const std::vector<mcutils::bit_tuple<uint8_t,1,1,1,1>>& state_label_vector)
    : BaseSubspace{T0}
    {
      for(const auto& state_labels : state_label_vector)
        PushStateLabels(state_labels);
    }

    //accessor
    unsigned int T0() const {return std::get<0>(labels());}
    // // diagnostic output
    std::string DebugStr() const;

  // private:
  };


  class OperatorState
      : public basis::BaseState<OperatorT0Subspace>
    {
     public:

      OperatorState() = default;
      // pass-through constructors

      // Construct state by index.
      OperatorState(const SubspaceType& subspace, std::size_t index)
          : basis::BaseState<OperatorT0Subspace>(subspace, index){}

      // Construct state by reverse lookup on labels.
      OperatorState(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        ): basis::BaseState<OperatorT0Subspace>(subspace, state_labels){}

      unsigned int Sbar() const {return static_cast<unsigned int>(mcutils::get<0>(labels()));}
      unsigned int Sbarp()const {return static_cast<unsigned int>(mcutils::get<1>(labels()));}
      unsigned int Tbar() const {return static_cast<unsigned int>(mcutils::get<2>(labels()));}
      unsigned int Tbarp()const {return static_cast<unsigned int>(mcutils::get<3>(labels()));}

      std::string LabelStr() const;

      // private:
    };


 class OperatorS0Subspace
  : public basis::BaseSpace<
    OperatorS0Subspace,
    OperatorT0Subspace,
    std::tuple<uint8_t,uint8_t>
  >
  {
   public:

    OperatorS0Subspace() = default;

    OperatorS0Subspace(
      const uint8_t S0,
      const uint8_t exchange_symm_bar,
      const u3shell::relative::OperatorParameters& operator_parameters,
      const std::array<std::vector<mcutils::bit_tuple<uint8_t,1,1,1,1>>,3>& operator_map
    ) :BaseSpace{{exchange_symm_bar,S0}}
    {
      // Restricted set of T0 values either given by Allowed_T0_values of
      // Include all possible values of T0.
      auto Allowed_T0_values = operator_parameters.Allowed_T0_values.size()?
        operator_parameters.Allowed_T0_values : std::set<uint8_t>{0,1,2};

      for(const auto& T0 : Allowed_T0_values)
        {
          const auto& state_label_vector = operator_map[T0];

          OperatorT0Subspace subspace(T0,state_label_vector);
          if(subspace.size()==0)
            continue;

          PushSubspace(std::move(subspace));
        }
    }

    // accessors
    unsigned int S0() const {return std::get<1>(labels());}
    unsigned int exchange_symm_bar() const {return std::get<0>(labels());}

    // // diagnostic output
    std::string LabelStr() const {return fmt::format("[{} {}]",S0(),exchange_symm_bar());}
    std::string DebugStr() const;
  // private:
  };

 class OperatorSpace
  : public basis::BaseSpace<OperatorSpace, OperatorS0Subspace>
  {
   public:
    OperatorSpace() = default;

    OperatorSpace(const u3shell::relative::OperatorParameters& operator_parameters);
    // Construct operator space from list of operator labels


    // // diagnostic output
    std::string DebugStr() const;

  // private:
  };

}// spin::twobody namespace



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
