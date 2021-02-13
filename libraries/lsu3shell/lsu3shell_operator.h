/****************************************************************
  lsu3shell_operator.h

  Generation of operator files for input to lsu3shell recoupler.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.
  12/2/16 (aem): Added bool for U(N)->U(3) restriction on operators
  2/23/17 (aem): Added openMP to GenerateLSU3ShellOperator
****************************************************************/

#ifndef LSU3SHELL_OPERATOR_H_
#define LSU3SHELL_OPERATOR_H_

#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/relative_cm_xform.h" 

namespace lsu3shell
{
  void GenerateModelSpaceFile(
      const std::string& model_space_filename, int Z, int N, int Nmax, int parity
    );
  // Generates model space for lsu3shell SU3RME with no restrictions on 
  // J taking full parity model space up to Nmax
  //
  // Arguments:
  //  model_space_filename (str) : model space filename
  //  Z,N (input) : number of protons/neutrons
  //  Nmax (input) : Truncation of oscillator quanta
  //  parity (input) : determines which parity space.
  //    parity=0 is natural parity, basis starts with Nex=0 
  //    parity=1 is unnatural parity, basis starts with Nex=1
  //    parity=-1 basis contains both parities and starts with Nex=0 

  void
  GenerateLSU3ShellOperator(
      int Nmax, 
      const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_tensor_expansion,
      std::string filename,
      bool un_u3_restrict=false
    );
  // Generates input files for LSU3Shell RecoupleSU3Interaction for an intrinsic or
  // relative operator.  Input files are save in "filename"
  //
  // Arguments:
  //    Nmax (input): basis truncation parameter
  //    relative_tensor_expansion (input) : relative unit tensor expansion of 
  //      relative or intrinsic operator
  //    un_u3_restrict (input) : if dealing with nucleus which will have single
  //      species only in a given shell, apply un->u3 branching restriction


  void
  GenerateLSU3ShellOperator(
      int Nmax, 
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_tensor_labels,
      bool un_u3_restrict=false  
    );
  // Overload of previous funciton
  //
  // Generate input files for LSU3Shell RecoupleSU3Interaction for all 
  // relative unit tensors tensor's which may have non-zero matrix elements
  // between  LGI's. 
  //
  // The "ith" unit tensor input file is written to "relative_unit_00000i.recoupler"
  //
  // Arguments:
  //    relative_tensor_labels (input) : list of relative unit tensor labels whose
  //      two-body biquad representation is to be generated

  void 
    GenerateLSU3ShellOperator(
        int Nmax, 
        const u3shell::TwoBodyUnitTensorCoefficientsU3ST& twobody_tensor_expansion,
        int operator_index,
        bool un_u3_restrict=false
      );
  // TODO: come up with better output files name 
    
  // Overload of previous funciton
  //
  // Generate input files for LSU3Shell RecoupleSU3Interaction for a two-body operator
  //
  // Currently writes input files to "two_body_unit_i.recoupler", where i is the operator index
  //
  // Arguments:
  //    relative_tensor_labels (input) : list of relative unit tensor labels whose
  //      two-body biquad representation is to be generated


}
#endif
