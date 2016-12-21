/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.

****************************************************************/

#include "sp3rlib/sp3r.h"

#include <iostream>

int main(int argc, char **argv)
{


  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////

  // 16(2,1)  Nn_max=6

  u3::U3 sigma = u3::U3(16,u3::SU3(2,1));
  int Nn_max = 6;

  // examine raising polynomials
  //
  // We generate and print them here for testing purposes only.
  // Normally sp3r::RaisingPolynomialLabels is called directly by the
  // Sp3RSpace constructor.

  std::vector<u3::U3> polynomial_labels = sp3r::RaisingPolynomialLabels(Nn_max);
  for (auto it = polynomial_labels.begin(); it != polynomial_labels.end(); ++it)
    std::cout << " label " << (*it).Str() << std::endl;
    
  // examine Sp3RSpace object
  sp3r::Sp3RSpace irrep(sigma,Nn_max);
  std::cout << irrep.DebugString();
  std::cout<<irrep.size()<<std::endl;
} //main
