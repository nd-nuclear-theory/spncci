/****************************************************************
  lsu3shell_operator.h

  Generation of operator files for input to lsu3shell recoupler.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.
  12/2/16 (aem): Added bool for U(N)->U(3) restriction on operators
****************************************************************/

#ifndef LSU3SHELL_OPERATOR_H_
#define LSU3SHELL_OPERATOR_H_

#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/relative_cm_xform.h" 

namespace lsu3shell
{
  void
  GenerateLSU3ShellOperator(
      int Nmax, 
      const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_tensor_expansion,
      std::string filename,
      bool un_u3_restrict=false
    );
  // Generates input files for LSU3Shell RecoupleSU3Interaction
  //
  // Arguments:

  void
  GenerateLSU3ShellOperator(
      int Nmax, 
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_tensor_labels,
      bool un_u3_restrict=false  
    );
  // Generate input files for LSUshell recoupler for all relative unit tensors 
  // tensor's which may have non-zero matrix elements between  LGI's. 

  void 
    GenerateLSU3ShellOperator(
        int Nmax, 
        const u3shell::TwoBodyUnitTensorCoefficientsU3ST& twobody_tensor_expansion,
        int operator_index,
        bool un_u3_restrict=false
      );

}
#endif
