/****************************************************************
  operator_indexing_spatial.h

  Indexing for operators defined by matrix elements in a harmonic
  oscillator basis.

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/23/22 (aem): Created.
  5/9/22 (aem): Templatized operator space on coordinate type.
****************************************************************/

#ifndef OPERATOR_INDEXING_SPATIAL_H_
#define OPERATOR_INDEXING_SPATIAL_H_

#include <unordered_set>

#include "basis/basis.h"
#include "basis/degenerate.h"
#include "sp3rlib/u3.h"
#include "u3shell/operator_parameters.h"

namespace u3shell::spatial
{

using OneCoordType = std::array<unsigned int, 1>;

void GenerateSpatialOperators(
    const int N0,
    const uint8_t parity_bar,
    const u3shell::relative::OperatorParameters& operator_parameters,
    std::map<u3::SU3, std::vector<OneCoordType>>& x0_Nbar_vector
  );


inline void GenerateSpatialOperators(
    const int N0,
    const uint8_t parity_bar,
    const unsigned int Nbar_max,
    std::map<u3::SU3, std::vector<OneCoordType>>& x0_Nbar_vector
  )
{
  return GenerateSpatialOperators(
      N0, parity_bar, {Nbar_max, u3shell::relative::kNone, {}, {}, {}, {}}, x0_Nbar_vector
    );
}

//////////////////////////////////////////////////////////////////////////////////////
// Space
//   OperatorSpace
//     -> OperatorParitySpace [parity_bar]
//       -> OperatorN0Space [N0]
//         -> OperatorL0Space [L0]
//           -> OperatorSU3Subspace (kappa0)[omega0]
//              ->OperatorStates [Nbar]  Nbarp fixed by N0.  (templatize based
//              on Nbar or Nbar1,Nbar2,Nbar3,Nbar4)
//
/////////////////////////////////////////////////////////////////////////////////////

// tOperatorStateLabelType is expected to be a tuple of integral type
template<typename tOperatorStateLabelType> class OperatorNbarStates;
template<typename tOperatorStateLabelType> class OperatorSU3Subspace;
template<typename tOperatorStateLabelType> class OperatorL0Space;
template<typename tOperatorStateLabelType> class OperatorN0Space;
template<typename tOperatorStateLabelType> class OperatorParitySpace;
template<typename tOperatorStateLabelType> class OperatorSpace;

////////////////////////////////////////////////////////////////
/// Subspace, labeled by x0, which contains states with label Nbar.
/// As a convenience N0 is also stored.
/// For each state, Nbarp obtained as Nbar+N0.
////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class OperatorU3Subspace
    : public basis::BaseSubspace<
          OperatorU3Subspace<tOperatorStateLabelType>,
          std::tuple<u3::SU3>,
          OperatorNbarStates<tOperatorStateLabelType>,
          std::tuple<tOperatorStateLabelType>
        >
{
 private:
  using BaseSubspaceType = basis::BaseSubspace<
      OperatorU3Subspace<tOperatorStateLabelType>,
      std::tuple<u3::SU3>,
      OperatorNbarStates<tOperatorStateLabelType>,
      std::tuple<tOperatorStateLabelType>
    >;

  using OperatorStateLabelType = tOperatorStateLabelType;

 public:
  /// Default OperatorU3Subspace constructor.
  OperatorU3Subspace() = default;


  /// OperatorU3Subspace constructor.
  /// Input:
  ///   N0 : N0 of operator
  ///   x0 : SU(3) tensor character of operator
  ///   Nbar_vector : Vector containing allowed Nbar values
  ///     for the operator subject to the contraint that
  ///     (Nbarp,0)x(0,Nbar)->x0, where Nbarp=Nbar+N0.
  ///     Can be obtained using function GenerateSpatialOperators.
  OperatorU3Subspace(
      const int N0,
      const u3::SU3& x0,
      const std::vector<OperatorStateLabelType>& Nbar_vector
    )
      : BaseSubspaceType{x0}, N0_{N0}
  {
    for (const auto& Nbar_labels : Nbar_vector)
    {
      BaseSubspaceType::PushStateLabels(Nbar_labels);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////////////////////
  /// Returns SU(3) label x0=$(\lambda_0,\mu_0)$ of subspace
  u3::SU3 x0() const { return std::get<0>(BaseSubspaceType::labels()); }

  /// Returns N0 of operator
  int N0() const { return N0_; }

  std::string LabelStr() const { return x0().Str(); }

  std::string DebugStr(const std::string& indent = "") const
  {
    std::string debug_str = fmt::format("{}x0: {}\n", indent, x0());

    std::string state_indent = indent + "  ";

    for (std::size_t i_state = 0; i_state < BaseSubspaceType::size(); ++i_state)
    {
      const auto& NbarState = BaseSubspaceType::GetState(i_state);

      debug_str += fmt::format("{}{}", i_state, NbarState.LabelStr());
    }

    return debug_str;
  }


 private:
  int N0_;
};


////////////////////////////////////////////////////////////////
// States given by Nbar,Nbar' pairs
////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class OperatorNbarStates
    : public basis::BaseState<OperatorU3Subspace<tOperatorStateLabelType>>
{
 private:
  using BaseStateType =
      basis::BaseState<OperatorU3Subspace<tOperatorStateLabelType>>;

  using OperatorStateLabelType = tOperatorStateLabelType;

 public:
  using SubspaceType = typename BaseStateType::SubspaceType;

  OperatorNbarStates() = default;
  // pass-through constructors

  OperatorNbarStates(const SubspaceType& subspace, std::size_t index)
  // Construct state by index.
      : BaseStateType(subspace, index)
  {}

  OperatorNbarStates(
      const SubspaceType& subspace,
      const typename BaseStateType::StateLabelsType& state_labels
    )
  // Construct state by reverse lookup on labels.
      : BaseStateType(subspace, state_labels)
  {}

  ////////////////////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////////////////////
  inline const OperatorStateLabelType& Nbar_labels() const
  {
    return std::get<0>(BaseStateType::labels());
  }

  inline unsigned int N1bar() const { return Nbar_labels()[0]; }

  template<
    typename U = tOperatorStateLabelType,
    std::enable_if_t<(std::tuple_size<U>::value >= 2)>* = nullptr
  >
  inline unsigned int N2bar() const { return Nbar_labels()[1]; }

  template<
    typename U = tOperatorStateLabelType,
    std::enable_if_t<(std::tuple_size<U>::value >= 3)>* = nullptr
  >
  inline unsigned int N1barp() const { return Nbar_labels()[2]; }

  template<
    typename U = tOperatorStateLabelType,
    std::enable_if_t<(std::tuple_size<U>::value >= 4)>* = nullptr
  >
  inline unsigned int N2barp() const { return Nbar_labels()[3]; }


  // Wrapper for N1bar() one coordinate
  template<
    typename U = tOperatorStateLabelType,
    std::enable_if_t<(std::tuple_size<U>::value <= 2)>* = nullptr
  >
  inline unsigned int Nbar() const { return N1bar(); }

  // Wrapper for N2bar() one coordinate
  template<
    typename U = tOperatorStateLabelType,
    std::enable_if_t<(std::tuple_size<U>::value == 2)>* = nullptr
  >
  inline unsigned int Nbarp() const { return N2bar(); }


  // Reconstruct Nbarp from Nbar and N0;
  template<
    typename U = tOperatorStateLabelType,
    std::enable_if_t<(std::tuple_size<U>::value == 1)>* = nullptr
  >
  inline unsigned int Nbarp() const
  {
    return static_cast<unsigned int>(Nbar() + BaseStateType::subspace().N0());
  }

  std::string LabelStr() const
  {
    std::string label_str = fmt::format("{:2d}", Nbar_labels()[0]);
    for (std::size_t i = 1; i < Nbar_labels().size(); ++i)
      label_str += fmt::format(" {:2d}", Nbar_labels()[i]);

    return fmt::format("[{}]", label_str);
  }

};

////////////////////////////////////////////////////////////////
// L0 subspace: constains x0 subspaces with degeneracy kappa0_max
////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class OperatorL0Space
    : public basis::BaseDegenerateSpace<
          OperatorL0Space<tOperatorStateLabelType>,
          OperatorU3Subspace<tOperatorStateLabelType>,
          std::tuple<unsigned int>
        >
{
 private:
  using BaseDegenerateSpaceType = basis::BaseDegenerateSpace<
      OperatorL0Space<tOperatorStateLabelType>,
      OperatorU3Subspace<tOperatorStateLabelType>,
      std::tuple<unsigned int>
    >;

 public:
  // constructors
  OperatorL0Space() = default;

  OperatorL0Space(
      const int N0,
      const unsigned int L0,
      const std::map<u3::SU3, std::vector<tOperatorStateLabelType>>& x0_Nbar_vector
    )
      : BaseDegenerateSpaceType{L0}
  {
    for (const auto& [x0, Nbar_vector] : x0_Nbar_vector)
    {
      unsigned int kappa0_max =
          L0 == u3shell::relative::kNone ? 1
                                         : u3::BranchingMultiplicitySO3(x0, L0);

      auto subspace =
          OperatorU3Subspace<tOperatorStateLabelType>(N0, x0, Nbar_vector);
      if (subspace.size() > 0 && kappa0_max > 0)
        BaseDegenerateSpaceType::PushSubspace(std::move(subspace), kappa0_max);
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////////////////////
  unsigned int L0() const
  {
    return std::get<0>(BaseDegenerateSpaceType::labels());
  }
  unsigned int kappa0_max(std::size_t i) const
  {
    return BaseDegenerateSpaceType::GetSubspaceDegeneracy(i);
  }

  std::string LabelStr() const
  {
    if (L0() == u3shell::relative::kNone)
      return "None";
    else
      return fmt::format("{}", L0());
  }

  std::string DebugStr(const std::string& indent) const
  {
    std::string debug_str = fmt::format("{}L0: {}\n", indent, LabelStr());

    std::string subspace_indent = indent + "  ";
    for (std::size_t i_subspace = 0; i_subspace < BaseDegenerateSpaceType::size();
         ++i_subspace)
    {
      const auto& subspace = BaseDegenerateSpaceType::GetSubspace(i_subspace);

      debug_str += fmt::format(
          "{}index: {} degeneracy: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          BaseDegenerateSpaceType::GetSubspaceDegeneracy(i_subspace),
          subspace.size(),
          subspace.dimension()
        );

      debug_str += subspace.DebugStr(subspace_indent);
    }

    return debug_str;
  }
  // private:
};

//////////////////////////////////////////////////////////////////////////////////////
/// OperatorN0Space: operators grouped by N0. Contains OperatorL0Spaces
/////////////////////////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class OperatorN0Space
    : public basis::BaseSpace<
          OperatorN0Space<tOperatorStateLabelType>,
          OperatorL0Space<tOperatorStateLabelType>,
          std::tuple<int>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<
      OperatorN0Space<tOperatorStateLabelType>,
      OperatorL0Space<tOperatorStateLabelType>,
      std::tuple<int>
    >;

 public:
  OperatorN0Space() = default;

  OperatorN0Space(
      const uint8_t parity_bar,
      const int N0,
      const u3shell::relative::OperatorParameters& operator_parameters
    )
      : BaseSpaceType{N0}
  {
    // Extract S0max from operator_parameters.
    // if no S0 values given, then S0max is default value
    unsigned int S0max = 2;
    const auto& allowed_S0_values = operator_parameters.Allowed_S0_values;
    if (allowed_S0_values.size() > 0)
      // allowed_S0_values are in ordered set with largest value being the right-most
      S0max = *(allowed_S0_values.rbegin());

    // Extract allowed values for L0 from operator_parameters.
    std::set<unsigned int> allowed_L0_values;
    if (operator_parameters.Allowed_L0_values.size() > 0)
      allowed_L0_values = operator_parameters.Allowed_L0_values;
    else
    {
      unsigned int L0min = static_cast<unsigned int>(
          std::max(0, int(operator_parameters.J0 - S0max))
        );

      for (unsigned int L0 = L0min; L0 <= operator_parameters.J0 + S0max; ++L0)
      {
        allowed_L0_values.insert(L0);
      }
    }

    std::map<u3::SU3, std::vector<tOperatorStateLabelType>> x0_Nbar_vector;
    u3shell::spatial::GenerateSpatialOperators(
        N0, parity_bar, operator_parameters, x0_Nbar_vector
      );

    // Iterate over L0 values and push L0Subspaces
    for (const unsigned int L0 : allowed_L0_values)
    {
      auto subspace =
          OperatorL0Space<tOperatorStateLabelType>(N0, L0, x0_Nbar_vector);

      if (subspace.size() > 0)
        BaseSpaceType::PushSubspace(std::move(subspace));
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////////////////////
  int N0() const { return std::get<0>(BaseSpaceType::labels()); }

  std::string LabelStr() const { return fmt::format("{}", N0()); }

  std::string DebugStr(const std::string& indent = "") const
  {
    std::string debug_str = fmt::format("{}N0: {}\n", indent, N0());

    std::string subspace_indent = indent + "  ";
    for (std::size_t i_subspace = 0; i_subspace < BaseSpaceType::size();
         ++i_subspace)
    {
      const auto& subspace = BaseSpaceType::GetSubspace(i_subspace);

      debug_str += fmt::format(
          "{}index: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          subspace.size(),
          subspace.dimension()
        );

      debug_str += subspace.DebugStr(subspace_indent);
    }

    return debug_str;
  }
};


template<typename tOperatorStateLabelType>
class OperatorParitySpace
    : public basis::BaseSpace<
          OperatorParitySpace<tOperatorStateLabelType>,
          OperatorN0Space<tOperatorStateLabelType>,
          std::tuple<uint8_t>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<
      OperatorParitySpace<tOperatorStateLabelType>,
      OperatorN0Space<tOperatorStateLabelType>,
      std::tuple<uint8_t>
    >;

 public:
  OperatorParitySpace() = default;

  OperatorParitySpace(
      const uint8_t parity_bar,
      const u3shell::relative::OperatorParameters& operator_parameters
    )
      : BaseSpaceType{parity_bar}
  {
    // For parity conserving Nbar and Nbarp
    // If Nmax even,
    //    If parity_bar=1, then max(Nbarp-Nbar)=Nbar_max-2
    //    If parity_bar=0, then max(Nbarp-Nbar)=Nbar_max
    // If Nmax odd
    //    then max(Nbarp-Nbar) = Nbar_max-1
    int N0_max = static_cast<int>(
        operator_parameters.Nbar_max - 2 * parity_bar
        + (operator_parameters.Nbar_max % 2)
      );
    for (int N0 = -N0_max; N0 <= N0_max; N0 += 2)
    {
      auto subspace = OperatorN0Space<tOperatorStateLabelType>(
          parity_bar, N0, operator_parameters
        );
      if (subspace.size() > 0)
        BaseSpaceType::PushSubspace(std::move(subspace));
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////////////////////
  uint8_t parity_bar() const { return std::get<0>(BaseSpaceType::labels()); }

  std::string LabelStr() const { return fmt::format("{}", parity_bar()); }

  std::string DebugStr(const std::string& indent) const
  {
    std::string debug_str = fmt::format("parity_bar: {}\n", parity_bar());

    for (std::size_t i_subspace = 0; i_subspace < BaseSpaceType::size();
         ++i_subspace)
    {
      const auto& subspace = BaseSpaceType::GetSubspace(i_subspace);

      std::string subspace_indent = indent + "  ";
      debug_str += fmt::format(
          "{}index: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          subspace.size(),
          subspace.dimension()
        );

      debug_str += subspace.DebugStr(subspace_indent);
    }

    return debug_str;
  }
};


//////////////////////////////////////////////////////////////////////////////////////
/// Space
///   OperatorSpace
///     -> OperatorParitySpace [parity_bar]
///       -> OperatorN0Space [N0]
///         -> OperatorL0Space [L0]
///           -> OperatorSU3Subspace (kappa0)[omega0]
///              ->OperatorStates [Nbar]  Nbarp fixed by N0.  (templatize based
///              on Nbar or Nbar1,Nbar2,Nbar3,Nbar4)
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
/// If the Allowed_{qn}_values.size()==0, then all possible values of qn are
/// considered. For e.g., ab initio Hamiltonin, operator_parameters would have
///
///    Nbar_max=Nmax+2*N1v
///    J0=0
///    Allowed_L0_values={} -> L0min=max(0,J0-2), L0=L0min,L0min+1,...,L0min+4;
///    Allowed_S0_values={} -> T0={0,1,2}
///    Allowed_T0_values={} -> T0={0,1,2}
///    Allowed_w0_values={} -> Given by GenerateSpatialOperators
///
/////////////////////////////////////////////////////////////////////////////////////
template<typename tOperatorStateLabelType>
class OperatorSpace
    : public basis::BaseSpace<
          OperatorSpace<tOperatorStateLabelType>,
          OperatorParitySpace<tOperatorStateLabelType>
        >
{
 private:
  using BaseSpaceType = basis::BaseSpace<
      OperatorSpace<tOperatorStateLabelType>,
      OperatorParitySpace<tOperatorStateLabelType>
    >;

 public:
  OperatorSpace() = default;

  OperatorSpace(const u3shell::relative::OperatorParameters& operator_parameters)
      : BaseSpaceType{}
  {
    for (uint8_t parity_bar : {0, 1})
    {
      auto subspace = OperatorParitySpace<tOperatorStateLabelType>(
          parity_bar, operator_parameters
        );
      if (subspace.size() > 0)
        BaseSpaceType::PushSubspace(std::move(subspace));
    }
  }

  std::string DebugStr() const
  {
    std::string debug_str;
    for (std::size_t i_subspace = 0; i_subspace < BaseSpaceType::size();
         ++i_subspace)
    {
      const auto& subspace = BaseSpaceType::GetSubspace(i_subspace);

      std::string indent = "";
      debug_str += fmt::format(
          "{}index: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          subspace.size(),
          subspace.dimension()
        );

      debug_str += subspace.DebugStr(indent);
    }
    return debug_str;
  }
};

}  // namespace u3shell::spatial

namespace u3shell::spatial::onecoord
{

// using OperatorSpace = u3shell::spatial::OperatorSpace<u3shell::spatial::OneCoordType>;
// using OperatorL0Space = u3shell::spatial::OperatorL0Space<OneCoordType>;
}  // namespace u3shell::spatial::onecoord


#endif
