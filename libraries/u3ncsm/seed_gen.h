/****************************************************************
  seed_gen.h

  Computes seeds for spncci recurrence using lsu3shell

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  2/24/22 (aem): Created.
****************************************************************/

#ifndef SEED_GEN_H_
#define SEED_GEN_H_

#include <unordered_map>

#include "LSU3/ncsmSU3xSU2Basis.h"
#include "am/halfint.h"
#include "fmt/format.h"
#include "lgi/lgi.h"
#include "lgi/recurrence_lgi.h"
#include "spncci_basis/recurrence_indexing.h"
#include "u3ncsm/u3ncsm_interface.h"
#include "u3shell/relative_operator.h"
#include "utilities/nuclide.h"

namespace spncci::seeds
{

basis::OperatorBlock<double> GenerateRecurrenceSeedBlock(
    const nuclide::NuclideType& nuclide,
    const HalfInt& Nsigma0,
    const int N1v,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
    const spncci::spin::RecurrenceSpace<lgi::LGI, u3shell::spin::twobody::OperatorLabelsST>&
        spin_recurrence_space,
    const spncci::spatial::RecurrenceSpace<u3shell::spatial::OneCoordType>& spatial_recurrence_space,
    const CBaseSU3Irreps& baseSU3Irreps,
    const std::string& operator_dir,
    const int recurrence_sp3r_space_index
  );


}  // namespace spncci::seeds
#endif
