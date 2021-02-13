/****************************************************************
  spncci_common.h

  Some globally configurable numerical typedefs and configuration parameters for SpNCCI.

  Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  2/27/17 (mac): Created.
  7/5/17 (mac): Overhaul switchable matrix precision typedefs.
****************************************************************/

#ifndef SPNCCI_COMMON_H_
#define SPNCCI_COMMON_H_

#include <string>

#include "eigen3/Eigen/Core"

#include "basis/hypersector.h"
#include "basis/operator.h"

namespace spncci
{
  // matrix precision
  typedef double MatrixFloatType;
  typedef Eigen::Matrix<MatrixFloatType,Eigen::Dynamic,1> Vector;   // e.g., for vector of eigenvalues
  typedef Eigen::Matrix<MatrixFloatType,Eigen::Dynamic,Eigen::Dynamic> Matrix;  // for matrix which is not semantically an operator block, e.g., matrix of eigenvectors
  typedef basis::OperatorBlock<MatrixFloatType> OperatorBlock;
  typedef basis::OperatorBlocks<MatrixFloatType> OperatorBlocks;
  typedef basis::OperatorHyperblocks<MatrixFloatType> OperatorHyperblocks;

  // tolerance for zero-testing of SpNCCI matrix elements
  extern MatrixFloatType g_zero_tolerance;
  extern bool g_suppress_zero_sectors;

}  // namespace

#endif
