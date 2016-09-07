/****************************************************************
  lsu3shell_rme.h

  Input of RMEs in lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.
****************************************************************/

#ifndef LSU3SHELL_RME_H_
#define LSU3SHELL_RME_H_

#include "boost/functional/hash_fwd.hpp"
#include "eigen3/Eigen/Eigen"

#include "am/am.h"  
#include "basis/operator.h"
#include "u3shell/u3spn_scheme.h"

#include "lsu3shell/lsu3shell_basis.h"

namespace lsu3shell
{

  void 
  ReadLSU3ShellRMEs(
      std::ifstream& is,
      const u3shell::OperatorLabelsU3S& operator_labels,
      const LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& sectors,
      basis::MatrixVector& matrix_vector // in operator.h and initial to zero
    );


}
#endif
