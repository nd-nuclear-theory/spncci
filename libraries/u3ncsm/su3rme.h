/****************************************************************
  su3rme.h

  Calculation of SU(3)xSU(2) reduced RMEs using lsu3shell

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  - 10/15/21 (aem): Created based on lsu3shell/tools/SU3RME_MPI.cpp
  
****************************************************************/
#ifndef SU3RME_H_
#define SU3RME_H_

#include "basis/operator.h"

// #include "SU3ME/proton_neutron_ncsmSU3Basis.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "UNU3SU3/UNU3SU3Basics.h"
namespace lsu3shell
{


  unsigned int get_num_U3PNSPN_irreps(const lsu3::CncsmSU3xSU2Basis& basis);

  basis::OperatorBlocks<double> CalculateRME( 
    const CInteractionPPNN& interactionPPNN,
    const CInteractionPN& interactionPN,
    const lsu3::CncsmSU3xSU2Basis& bra,
    const lsu3::CncsmSU3xSU2Basis& ket,
    int dN0,
    const SU3xSU2::LABELS& w0,
    int rhot_max
  );

}// end namespace 

#endif
