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

namespace u3shell::spatial
{

void GenerateSpatialOperators(
    const int N0,
    const uint8_t parity_bar,
    const u3shell::relative::OperatorParameters& operator_parameters,
    std::map<u3::SU3, std::vector<OneCoordType>>& x0_Nbar_vector
  )
{
  unsigned int Nbar_min = parity_bar;
  unsigned int Nbar_max = operator_parameters.Nbar_max;
  const auto& Allowed_w0_values = operator_parameters.Allowed_w0_values;
  const auto& Allowed_L0_values = operator_parameters.Allowed_L0_values;
  bool restrict_x0 = Allowed_w0_values.size() > 0;
  bool restrict_by_L0 = (Allowed_L0_values.size() != 0);
  // Check that the last value is not kNone.  kNone is defined as max value of
  // unsigned int, so it will always be the last value of the set if in the set.
  if (restrict_by_L0)
    restrict_by_L0 &= (*Allowed_L0_values.rbegin()) != u3shell::relative::kNone;
  //////////////////////////////////////////////////////////////////////////////
  for (unsigned int Nbar = Nbar_min; Nbar <= Nbar_max; Nbar += 2)
  {
    // Nbarp must be greater than zero
    if ((N0 + static_cast<int>(Nbar)) < 0)
      continue;

    unsigned int Nbarp = static_cast<unsigned int>(N0 + Nbar);
    if (Nbarp <= Nbar_max)
    {
      MultiplicityTagged<u3::SU3>::vector possible_x0 =
          u3::KroneckerProduct({Nbarp, 0u}, {0u, Nbar});

      for (const auto& [x0, rho] : possible_x0)
      {
        bool keep_x0 = (!restrict_x0 || Allowed_w0_values.count({N0, x0}));

        // Check restrictions on x0 branching to L0 if any
        if (restrict_by_L0 && keep_x0)
        {
          keep_x0 = false;
          for (const auto& L0 : Allowed_L0_values)
          {
            if (u3::BranchingMultiplicitySO3(x0, L0) > 0)
            {
              keep_x0 = true;
              break;
            }
          }
        }
        // If x0 passed all checks, add to map
        if (keep_x0)
          x0_Nbar_vector[x0].push_back({Nbar});
      }
    }
  }
}

}  // namespace u3shell::spatial
