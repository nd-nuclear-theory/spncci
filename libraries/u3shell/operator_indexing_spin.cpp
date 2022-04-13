/****************************************************************
  operator_indexing_spin.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3shell/operator_indexing_spin.h"

#include "fmt/format.h"
#include <cppitertools/imap.hpp>
#include <cppitertools/unique_everseen.hpp>

namespace u3shell::spin::twobody{

std::string OperatorSubspace::DebugStr(const std::string& indent) const
{
  std::string debug_str = fmt::format(
      "{}S0Subspace: S0: {} parity_bar: {} size: {}  dimension: {}\n",
      indent,
      S0(),
      exchange_symm_bar(),
      size(),
      dimension()
    );

  std::string state_indent = indent + "  ";
  for (int i_subspace = 0; i_subspace < size(); ++i_subspace)
  {
    const auto& index = std::get<0>(GetStateLabels(i_subspace));

    const auto& labels = spin::twobody::OperatorLabelsST::ALLOWED_LABELS[index];
    debug_str += fmt::format(
        "{}[{} {} {} {} {} {}]\n",
        state_indent,
        labels.S0(),
        labels.T0(),
        labels.Sbar(),
        labels.Sbarp(),
        labels.Tbar(),
        labels.Tbarp()
      );
  }

  return debug_str;
}


std::string OperatorSpace::DebugStr() const
{
  std::string debug_str="";
  for (int i_subspace=0; i_subspace<size(); ++i_subspace)
    {
      const auto& subspace = GetSubspace(i_subspace);
      debug_str+=subspace.DebugStr();
    }
  return debug_str;
}

}





