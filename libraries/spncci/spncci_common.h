/****************************************************************
  spncci_common.h

  Some globally configurable numerical typedefs and configuration parameters for SpNCCI.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/27/17 (mac): Created.
****************************************************************/

#ifndef SPNCCI_COMMON_H_
#define SPNCCI_COMMON_H_

#include <string>

#include "eigen3/Eigen/Core"

namespace spncci
{
  // numerics
  typedef double MatrixFloatType;
  typedef Eigen::VectorXd VectorType;
  typedef Eigen::MatrixXd MatrixType;

  // tolerance for zero-testing of SpNCCI matrix elements
  extern MatrixFloatType g_zero_tolerance;
  extern bool g_suppress_zero_sectors;

  // logging
  void WriteLog(const std::string& message);

}  // namespace

#endif
