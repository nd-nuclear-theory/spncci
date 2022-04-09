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
