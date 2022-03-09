/****************************************************************
  u3boson_test.cpp

  Anna E. McCoy
  Institute for Nuclear Theory
  
  SPDX-License-Identifier: MIT 

  3/7/22 (aem): Created.

****************************************************************/
// #include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"
#include "fmt/format.h"
#include <iostream>
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3boson.h"


int main(int argc, char **argv)
{
  u3::U3CoefInit(39);

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////
  u3::U3 sigma = u3::U3(16,u3::SU3(2,1));
  int Nn_max = 10;

  std::vector<u3::U3> polynomial_labels = vcs::RaisingPolynomialLabels(Nn_max);
  for (const auto&n : polynomial_labels)
    fmt::print("{}\n",n);

  vcs::U3BosonSpace u3boson_space(sigma,20);
  std::cout<<u3boson_space.DebugStr()<<std::endl;
} //main


