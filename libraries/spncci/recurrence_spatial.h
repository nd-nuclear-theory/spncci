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
#include "recurrence_indexing.h"


namespace spncci::recurrence
{

class SpatialRecurrenceMatrix
{
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////
 public:
  SpatialRecurrenceMatrix() = default;

  SpatialRecurrenceMatrix(
      const spncci::spatial::RecurrenceSp3RSpace& recurrence_space
    );

  ////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////
 public:
  auto recurrence_space_ptr() const { return recurrence_space_; }
  const auto recurrence_space() const { return *recurrence_space_; }

  const std::vector<basis::OperatorBlock<double>>& GetRecurrenceBlock(unsigned int Nnsum) const
  {
    assert(recurrence_done_.at(Nnsum / 2));
    return recurrence_blocks_.at(Nnsum / 2);
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
  std::vector<basis::OperatorBlock<double>>& GetRecurrenceBlock(unsigned int Nnsum)
  {
    assert(recurrence_done_.at(Nnsum / 2));
    return recurrence_blocks_.at(Nnsum / 2);
  }

  ////////////////////////////////////////////////////////////////
  // local storage
  ////////////////////////////////////////////////////////////////

 private:
  std::shared_ptr<const spncci::spatial::RecurrenceSp3RSpace> recurrence_space_;
  std::vector<basis::OperatorBlocks<double>> recurrence_blocks_;
  const std::size_t lgi_dimension_;
#ifndef NDEBUG
  std::vector<bool> recurrence_done_;
#endif
};

}  // namespace spncci::recurrence
#endif  // SPNCCI_RECURRENCE_SPATIAL_
