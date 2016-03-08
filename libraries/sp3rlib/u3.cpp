/****************************************************************
  u3.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "sp3rlib/u3.h"

#include <sstream>


namespace u3 
{
  
  std::string SU3::Str() const
  {
    std::ostringstream ss;

    ss << "(" << lambda << "," << mu << ")";
    return ss.str();
  }

  std::string U3::Str() const
  {
    std::ostringstream ss;

    ss << N() << SU3().Str();

    // ss << "[" << f1 << "," << f2 << "," << f3 << "]";

    return ss.str();
  }

  std::string U3S::Str() const
  {
    std::ostringstream ss;

    ss << w.Str() << "x" << S;

    return ss.str();
  }

  std::string U3ST::Str() const
  {
    std::ostringstream ss;

    ss << w.Str() << "x" << S<< "x" << T;

    return ss.str();
  }



}  // namespace
