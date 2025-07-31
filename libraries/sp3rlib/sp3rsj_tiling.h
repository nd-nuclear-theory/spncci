/****************************************************************
  sp3rsj_tiling.h

  Construction and tiling of operator matrix in Sp3RSJ basis.

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  02/05/24 (cvc): Created.
  05/23/24 (cvc): Updated.
  06/03/24 (cvc): Added constructor for Sp3RSJ hyperblocks.
  06/03/24 (cvc): Refactored to .h/.cpp format.
  06/11/24 (cvc): Fixed looping over hypersectors.
  07/05/24 (cvc): Observables changed to linear combinations of operators.
  07/08/24 (cvc): ConstructSp3RHyperBlocks() computes Sp(3,R)xSU(2)>SU(2) RMEs.

****************************************************************/
#ifndef SP3RSJ_TILING_H_
#define SP3RSJ_TILING_H_

#include "sp3rlib/sp3r_operator.h"

#include <functional>
#include <fstream>
#include <tuple>
#include <vector>

#include "am/halfint.h"
#include "am/am.h"
#include "am/wigner_gsl.h"

#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "mcutils/eigen.h"

#include "sp3rlib/sp3r.h"
#include "sp3rlib/sp3rsj.h"
#include "basis/basis.h"
#include "basis/operator.h"
#include "basis/hypersector.h"

namespace sp3r
{
  
// Construct hyperblocks
void ConstructSp3RHyperBlocks(
  const sp3r::Sp3RSpace& sp3r_space,
  const sp3r::Sp3RSJSpace& bra_space,
  const sp3r::Sp3RSJSpace& ket_space,
  const sp3r::OperatorSpaceSp3RSJ& op_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  const int bra_index,
  const int ket_index,
  HalfInt& J_bra, HalfInt& J_ket, HalfInt& J_op,
  HalfInt& S_bra, HalfInt& S_ket, HalfInt& S_op,
  u3::UCoefCache& u_coef_cache,
  u3::WCoefCache& w_coef_cache,
  const sp3r::OperatorInfoSp3RSJ& operator_info,
  sp3r::Sp3RSJHypersectors& hypersectors,
  basis::OperatorHyperblocks<double>& operator_hyperblocks
);

// Get operator tile
void GetSp3ROperatorTile(
  const sp3r::Sp3RSpace& sp3r_space,
  const sp3r::Sp3RSJSpace& bra_space,
  const sp3r::Sp3RSJSpace& ket_space,
  const sp3r::OperatorSpaceSp3RSJ& op_space,
  const sp3r::U3Subspace& bra_subspace,
  const sp3r::U3Subspace& ket_subspace,
  const int bra_index,
  const int ket_index,
  HalfInt& J_bra, HalfInt& J_ket, HalfInt& J_op,
  HalfInt& S_bra, HalfInt& S_ket, HalfInt& S_op,
  u3::UCoefCache& u_coef_cache,
  u3::WCoefCache& w_coef_cache,
  const sp3r::OperatorInfoSp3RSJ& operator_info,
  sp3r::Sp3RSJHypersectors& hypersectors,
  basis::OperatorHyperblocks<double>& operator_hyperblocks,
  basis::OperatorBlock<double>& tile
);

// Construct Operator Matrix
void ConstructSp3ROperatorMatrix(
  const sp3r::Sp3RSpace& sp3r_space,
  const sp3r::Sp3RSJSpace& bra_space,
  const sp3r::Sp3RSJSpace& ket_space,
  const sp3r::ObservableSp3RSJ& observable_list,
  basis::OperatorBlock<double>& operator_matrix
);


}

#endif