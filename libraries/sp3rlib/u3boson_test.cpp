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
  int Nn_max = 6;

  // examine raising polynomials
  //
  // We generate and print them here for testing purposes only.
  // Normally sp3r::RaisingPolynomialLabels is called directly by the
  // Sp3RSpace constructor.
  std::vector<u3::U3> polynomial_labels = vcs::RaisingPolynomialLabels(Nn_max);
  for (const auto&n : polynomial_labels)
    fmt::print("{}\n",n);

  // Check that U(3) is unitary LGI
  std::vector<u3::U3>test_irreps ={
    {16,{2,1}},
    {3,{0,0}},
    {5,{2,0}},
    {6,5,1},
    {9,{6,0}}
  };
  
  
  for(const auto& sigma : test_irreps)
    {
      fmt::print("\n{}\n",sigma);
      vcs::U3BosonSpace irrep(sigma,Nn_max);
      std::cout << irrep.DebugStr();
    }

} //main


