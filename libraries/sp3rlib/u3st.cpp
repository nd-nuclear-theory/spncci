/****************************************************************
  u3st.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/

#include "sp3rlib/u3st.h"

#include <sstream>



namespace u3 
{
  std::string U3ST::Str() const
  {
    return fmt::format("{}x{}x{}",U3(),S(),T());
  }

}  // u3 namespace
