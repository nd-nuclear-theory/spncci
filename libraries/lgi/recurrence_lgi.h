/****************************************************************
  recurrence_lgi.h

  Functions necessary for constructing spncci recurrence indexing
  for specific choice of lgi (U3SPN)

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  02/23/22 (aem): Created.
****************************************************************/
#ifndef RECURRENCE_LGI_H_
#define RECURRENCE_LGI_H_

#include <tuple>

#include "am/am.h"
#include "am/halfint.h"
#include "lgi/lgi.h"
#include "spncci_basis/recurrence_indexing.h"
#include "u3shell/operator_indexing_spin.h"

namespace lgi
{
inline bool UnitTensorAllowed(
    const lgi::LGI::UpstreamLabelsType& ket_labels,
    const u3shell::spin::twobody::OperatorLabelsST& tensor_labels,
    const lgi::LGI::UpstreamLabelsType& bra_labels
  )
// Check if bra and ket labels can be connected by a given ST tensor.
//
// Note: these checks were derived by brute force enumeration of the possible
// $\bar{S}_p' \bar{S}_n' \bar{S}_p \bar{S}_n$ for states with given S, T and
// Tz, followed by checking all possible products and selecting those which
// satisfy all the relevant triangle constraints for coupling to good $S_{p,0}
// \times S_{n,0} \rightarrow S_0$.
{
  const HalfInt &bra_Sp = bra_labels.Sp, bra_Sn = bra_labels.Sn;
  const int S0 = tensor_labels.S0();
  const HalfInt &ket_Sp = ket_labels.Sp, ket_Sn = ket_labels.Sn;

  if (am::AllowedTriangle(bra_Sp, 1, ket_Sp) && am::AllowedTriangle(bra_Sn, 1, ket_Sn))
    return true;

  if (S0 == 0)
  {
    if (bra_Sp == ket_Sp && bra_Sn == ket_Sn)
      return true;
  }
  else if (S0 == 1)
  {
    if (bra_Sp == ket_Sp && am::AllowedTriangle(bra_Sn, 1, ket_Sn))
      return true;

    if (am::AllowedTriangle(bra_Sp, 1, ket_Sp) && bra_Sn == ket_Sn)
      return true;
  }
  else if (S0 == 2 && tensor_labels.Tbarp() == 1 && tensor_labels.Tbar() == 1)
  {
    if (bra_Sp == ket_Sp && am::AllowedTriangle(bra_Sn, 2, ket_Sn))
      return true;

    if (am::AllowedTriangle(bra_Sp, 2, ket_Sp) && bra_Sn == ket_Sn)
      return true;
  }

  return false;
}

}  // namespace lgi


#endif
