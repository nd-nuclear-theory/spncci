/****************************************************************
  operator_parameters.h

  Struct defining allowed quantum numbers for an operator
                                  
  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  4/6/22 (aem): Extracted from recurrence_indexing_operator.
****************************************************************/

#ifndef OPERATOR_PARAMETERS_H_
#define OPERATOR_PARAMETERS_H_

#include <unordered_set>
#include "sp3rlib/u3.h"


namespace u3shell::relative
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
      const std::set<uint8_t>& Allowed_S0_values_,
      const std::set<uint8_t>& Allowed_T0_values_
    )
        :
        Nbar_max{Nmax_+2*N1v_},
        J0{J0_},
        Allowed_w0_values{Allowed_w0_values_},
        Allowed_L0_values{Allowed_L0_values_},
        Allowed_S0_values{Allowed_S0_values_},
        Allowed_T0_values{Allowed_T0_values_}

    {
      for(const auto& S0 : Allowed_S0_values_)
        assert(S0==0 || S0==1 || S0==2);

      for(const auto& T0 : Allowed_T0_values_)
        assert(T0==0 || T0==1 || T0==2);
    }

    const int Nbar_max;
    const unsigned int J0;
    const std::set<unsigned int> Allowed_L0_values;
    const std::set<uint8_t> Allowed_S0_values;
    const std::set<uint8_t> Allowed_T0_values;
    const std::unordered_set<u3::U3> Allowed_w0_values;
  };

}//relative namespace



#endif
