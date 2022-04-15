/****************************************************************
  operator_indexing_spatial.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3shell/operator_indexing_spatial.h"

#include <cppitertools/imap.hpp>
#include <cppitertools/unique_everseen.hpp>

#include "fmt/format.h"

namespace u3shell
{
namespace spatial::onecoord
{
std::map<u3::SU3, std::vector<unsigned int>>
GenerateSpatialOperators(
    const int N0,
    const uint8_t parity_bar,
    const u3shell::relative::OperatorParameters& operator_parameters
  )
{
  unsigned int Nbar_min = parity_bar;
  unsigned int Nbar_max = operator_parameters.Nbar_max;
  const auto& Allowed_w0_values = operator_parameters.Allowed_w0_values;
  bool restrict_x0 = Allowed_w0_values.size() > 0;
  //////////////////////////////////////////////////////////////////////////////

  std::map<u3::SU3, std::vector<unsigned int>> x0_Nbar_vector;

  for (unsigned int Nbar = Nbar_min; Nbar <= Nbar_max; Nbar += 2) {
    // Nbarp must be greater than zero
    if ((N0 + static_cast<int>(Nbar)) < 0) continue;

    unsigned int Nbarp = N0 + Nbar;
    if (Nbarp <= Nbar_max) {
      MultiplicityTagged<u3::SU3>::vector possible_x0 =
          u3::KroneckerProduct({Nbarp, 0u}, {0u, Nbar});

      for (const auto& [x0, rho] : possible_x0) {
        if (!restrict_x0 || Allowed_w0_values.count({N0, x0}))
          x0_Nbar_vector[x0].push_back(Nbar);
      }
    }
  }

  return x0_Nbar_vector;
}


    std::string OperatorU3Subspace::DebugStr(const std::string& indent) const
    {
      std::string debug_str=fmt::format("{}x0: {}\n",indent,x0());

      std::string subspace_indent=indent+"  ";
      for (int i_subspace=0; i_subspace<size(); ++i_subspace)
      {
        const auto& NbarState = GetState(i_subspace);

        debug_str += fmt::format("{}[{} {}]\n", subspace_indent,
                                 NbarState.Nbar(), NbarState.Nbarp());
      }

      return debug_str;
    }


    std::string OperatorL0Space::DebugStr(const std::string& indent) const
    {
      std::string debug_str=fmt::format("{}L0: {}\n",indent,L0());

      std::string subspace_indent=indent+"  ";
      for (int i_subspace=0; i_subspace<size(); ++i_subspace)
      {
        const auto& subspace = GetSubspace(i_subspace);

        debug_str+=fmt::format("{}index: {} degeneracy: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          GetSubspaceDegeneracy(i_subspace),
          subspace.size(),
          subspace.dimension()
          );

        debug_str += subspace.DebugStr(subspace_indent);
      }

      return debug_str;
    }

    std::string OperatorN0Space::DebugStr(const std::string& indent) const
    {
      std::string debug_str=fmt::format("{}N0: {}\n",indent,N0());

      std::string subspace_indent=indent+"  ";
      for (int i_subspace=0; i_subspace<size(); ++i_subspace)
      {
        const auto& subspace = GetSubspace(i_subspace);

        debug_str+=fmt::format("{}index: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          subspace.size(),
          subspace.dimension()
          );

        debug_str += subspace.DebugStr(subspace_indent);
      }

      return debug_str;
    }

    std::string OperatorParitySpace::DebugStr(const std::string& indent) const
    {
      std::string debug_str=fmt::format("parity_bar: {}\n",parity_bar());

      for (int i_subspace=0; i_subspace<size(); ++i_subspace)
      {
        const auto& subspace = GetSubspace(i_subspace);

        std::string subspace_indent=indent+"  ";
        debug_str+=fmt::format("{}index: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          subspace.size(),
          subspace.dimension()
          );

        debug_str+=subspace.DebugStr(subspace_indent);

      }

      return debug_str;
    }

    std::string OperatorSpace::DebugStr() const
    {

      std::string debug_str;
      for (int i_subspace=0; i_subspace<size(); ++i_subspace)
      {
        const auto& subspace = GetSubspace(i_subspace);

        std::string indent="";
        debug_str+=fmt::format("{}index: {} size: {}  dimensions: {}\n",
          indent,
          i_subspace,
          subspace.size(),
          subspace.dimension()
          );

        debug_str+=subspace.DebugStr(indent);

      }

      return debug_str;
    }

}  // namespace spatial::onecoord
}//end u3shell namespace




