/****************************************************************
  spncci_common.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/spncci_common.h"

#include <fstream>
#include <iostream>

namespace spncci
{
  spncci::MatrixFloatType g_zero_tolerance = 1e-6;  

  bool g_suppress_zero_sectors = false;

  const std::string log_filename("spncci.log");
  std::ofstream log_stream(log_filename);
 
  void WriteLog(const std::string& message)
  {
    std::cout << message << std::endl;
    spncci::log_stream << message << std::endl;
  }

}  // namespace
