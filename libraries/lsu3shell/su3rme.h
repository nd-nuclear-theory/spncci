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

#include "lgi/lgi.h"
#include "SU3ME/CInteractionPN.h"
#include "SU3ME/InteractionPPNN.h"
#include "utilities/nuclide.h"


namespace lsu3shell
{


  using InteractionTerms = std::tuple<CInteractionPPNN,CInteractionPN>;

  InteractionTerms GetOperatorFromFile(
    const nuclide::NuclideType& nuclide,
    const int Nmax,
    const std::string& operator_filename
  );


}// end namespace 

#endif
