/****************************************************************
  recurrence_spatial.h

  Spatial recurrence matrix evaluation

  Patrick J. Fasano
  University of Notre Dame and LBNL

  SPDX-License-Identifier: MIT

  09/30/21 (pjf): Created.
****************************************************************/

#ifndef SPNCCI_RECURRENCE_SPATIAL_
#define SPNCCI_RECURRENCE_SPATIAL_

#include <memory>

#include "basis/basis.h"
#include "basis/operator.h"
#include "recurrence_indexing_spatial.h"
#include "vcs_cache.h"


namespace spncci::recurrence
{

class SpatialRecurrenceMatrix
{
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////
 public:
  SpatialRecurrenceMatrix()
      : lgi_recurrence_dimension_{0}
  {}

  SpatialRecurrenceMatrix(
      std::shared_ptr<const spncci::spatial::RecurrenceSp3RSpace> space_ptr
    );

  ////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////
 public:
  inline auto recurrence_space_ptr() const { return recurrence_space_ptr_; }
  inline const auto& recurrence_space() const { return *recurrence_space_ptr_; }

  inline const basis::OperatorBlocks<double>& GetRecurrenceBlock(unsigned int Nnsum) const
  {
#ifndef NDEBUG
    assert(recurrence_done_.at(Nnsum / 2));
#endif
    return recurrence_blocks_[Nnsum / 2];
  }

  ////////////////////////////////////////////////////////////////
  // computation
  ////////////////////////////////////////////////////////////////
 public:
  void GenerateRecurrenceBlock(unsigned int Nnsum);

  ////////////////////////////////////////////////////////////////
  // private accessors
  ////////////////////////////////////////////////////////////////
 private:
  inline basis::OperatorBlocks<double>& GetRecurrenceBlock(unsigned int Nnsum)
  {
#ifndef NDEBUG
    assert(recurrence_done_.size() > (Nnsum / 2));
#endif
    return recurrence_blocks_[Nnsum / 2];
  }

  ////////////////////////////////////////////////////////////////
  // local storage
  ////////////////////////////////////////////////////////////////

 private:
  std::shared_ptr<const spncci::spatial::RecurrenceSp3RSpace> recurrence_space_ptr_;
  std::vector<basis::OperatorBlocks<double>> recurrence_blocks_;
  const std::size_t lgi_recurrence_dimension_;
#ifndef NDEBUG
  std::vector<bool> recurrence_done_;
#endif
};

}  // namespace spncci::recurrence
#endif  // SPNCCI_RECURRENCE_SPATIAL_
