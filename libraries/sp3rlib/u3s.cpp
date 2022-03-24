/****************************************************************
  u3s.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/

#include "sp3rlib/u3s.h"

#include <sstream>



namespace u3 
{
 
  std::string U3S::Str() const
  {
    return fmt::format("{}x{}",U3(),S());
  }

}  // namespace
