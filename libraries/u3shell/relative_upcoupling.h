/****************************************************************
  relative_upcoupling.h

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/31/22 (aem): Created based on upcoupling
****************************************************************/

#ifndef U3SHELL_RELATIVE_UPCOUPLING_H_
#define U3SHELL_RELATIVE_UPCOUPLING_H_

#include "basis/basis.h"
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/operator_indexing_spin.h"
#include "u3shell/operator_indexing_spatial.h"
#include "u3shell/operator_indexing_sectors.h"


namespace u3shell::relative
{

  std::vector<double> UpcoupleU3ST(
      const unsigned int Nbar_max,
      const uint8_t g0,
      const u3shell::relative::OperatorSectors& sectors_u3st,
      const std::array<basis::RelativeSectorsLSJT,3>& component_sectors,
      std::array<basis::OperatorBlocks<double>,3>& component_blocks
    );

  std::pair<u3shell::relative::OperatorParameters,std::vector<double>> UpcoupleU3ST(
      const unsigned int Nbar_max,
      const std::string& input_filename,
      const unsigned int input_Nbar_max,
      const unsigned int input_Jmax
    );


  std::vector<double> UpcoupleU3ST(
      const unsigned int Nbar_max,
      const u3shell::relative::OperatorSectors& sectors_u3st,
      const std::string& input_filename,
      const unsigned int input_Nbar_max,
      const unsigned int input_Jmax
    );

}



#endif
