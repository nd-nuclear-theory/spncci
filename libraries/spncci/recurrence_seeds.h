/****************************************************************
recurrence_seeds.h

    Populating recurrence seed matricex 


  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT

  09/24/21 (aem) : Created.
****************************************************************/

#ifndef RECURRENCE_SEEDS_H_
#define RECURRENCE_SEEDS_H_

#include <array>
#include <functional>  // for std::hash
#include <map>
#include <tuple>
#include <utility>

#include "am/am.h"
#include "am/halfint.h"
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "lgi/lgi.h"
#include "sp3rlib/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "recurrence_indexing_spatial.h"
#include "recurrence_indexing_spin.h"

namespace spncci
{
namespace recurrence
{

basis::OperatorBlocks<double>
GetRecurrenceSeedsFromFile(
  const lgi::MultiplicityTaggedLGIVector& lgi_vector,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::spatial::RecurrenceSpace& spatial_recurrence_space,
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>& spin_recurrence_space
);

} // namespace recurrence 
}  // namespace spncci


#endif  // RECURRENCE_SEEDS_H_