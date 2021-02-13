/****************************************************************
  spncci_common.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/spncci_common.h"

#include <fstream>
#include <iostream>

namespace spncci
{
  spncci::MatrixFloatType g_zero_tolerance = 1e-6;  

  bool g_suppress_zero_sectors = false;

}  // namespace
