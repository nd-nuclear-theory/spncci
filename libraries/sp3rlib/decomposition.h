/****************************************************************
  decomposition.h

  Generates decompositions for sp3rlib.

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/

#ifndef SP3RLIB_DECOMPOSITION_H_
#define SP3RLIB_DECOMPOSITION_H_

#include <fstream>

#include "am/halfint.h"
#include "am/am.h"
#include "am/wigner_gsl.h"

#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "mcutils/io.h"

#include "sp3rlib/sp3r.h"
#include "sp3rlib/sp3rsj.h"

#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"

#include "basis/basis.h"
#include "basis/operator.h"
#include "basis/hypersector.h"

#include "sp3rlib/config_parameters.h"

namespace sp3r
{
  //
  void CalculateNexDecompositions(
    const std::vector<sp3r::Sp3RSJSpace>& spaces_sp3rsj,
    const std::vector<sp3r::Matrix>& eigenvectors,
    std::vector<sp3r::Matrix>& Nex_decompositions,
    HalfInt Nsigma0,
    int Nmax
  );
}

#endif