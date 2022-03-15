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

  std::vector<u3::U3> polynomial_labels = u3boson::RaisingPolynomialLabels(Nn_max);
  for (const auto&n : polynomial_labels)
    fmt::print("{}\n",n);


  fmt::print("sigma: {}\n",sigma);
  fmt::print("U3BosonSpace\n");
  u3boson::U3BosonSpace u3boson_space(sigma,20);
  std::cout<<u3boson_space.DebugStr()<<std::endl;

  u3::U3 n(2,0,0), np(4,0,0);
  fmt::print("BosonCreationRME ({}||a^dagger||{}\n",np,n);
  fmt::print("{}\n",u3boson::BosonCreationRME(np,n));

  u3::U3 omega(18,{4,1}), omegap(20,{4,2});
  fmt::print(
    "U3BosonCreationRME ({}{}{}{}||a^dagger||{}{}{}{})\n",
    sigma,np,1,omegap,sigma,n,1,omega
  );

 // Calculate with two different functions
  fmt::print(
    "{}\n{}\n",
    u3boson::U3BosonCreationRME(sigma,{np,1},omegap,sigma,{n,1},omega),
    u3boson::U3BosonCreationRME(sigma,np,1,omegap,sigma,n,1,omega)
  );

  fmt::print(
    "U3BosonAnnihilationRME ({}{}{}{}||a^dagger||{}{}{}{})\n",
    sigma,np,1,omegap,sigma,n,1,omega
  );

  // Calculate with two different functions
  fmt::print(
    "{}\n",
    u3boson::U3BosonAnnihilationRME(sigma,n,1,omega,sigma,np,1,omegap)
  );

  // Testing default constructor
  std::map<u3::U3,u3boson::U3BosonSpace> test_map;
  test_map[sigma]=u3boson_space;
  // =u3boson::U3BosonSpace(sigma,Nn_max);


} //main


