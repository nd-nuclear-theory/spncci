/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame
  
  SPDX-License-Identifier: MIT 

  3/10/16 (aem,mac): Created.

****************************************************************/
#include "sp3rlib/sp3r_operator.h"

#include <iostream>
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "fmt/format.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit(39);

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////
  // u3::U3 sigma = u3::U3(16,u3::SU3(2,1));
  int Nn_max = 6;

  std::vector<u3::U3>
  sigma_list = {
      u3::U3(16,{2,1}),
      u3::U3({29,2},{2,1})
    };

  for(const auto& sigma : sigma_list)
    {
      //TODO: Eventually replace by single function call
      u3boson::U3BosonSpace u3boson_space(sigma,Nn_max);
      sp3r::Sp3RSpace sp3r_space(sigma,Nn_max,u3boson_space);
    }

} //main


