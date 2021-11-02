/****************************************************************
  u3_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT 

  3/7/16 (aem,mac): Created.
  3/8/16 (aem,mac): Add tests for U3ST.
  3/9/16 (aem,mac): Update includes and namespaces.

****************************************************************/

#include "sp3rlib/u3.h"

#include <vector>
#include <algorithm>

#include "am/am.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"

int main(int argc, char **argv)
{
  ////////////////////////////////////////////////////////////////
  // construction, string conversion and comparisons
  ////////////////////////////////////////////////////////////////

  u3::SU3 x1(1,3);
  u3::SU3 x2(4,2);
  u3::SU3 x3(1,3);
  std::cout<<fmt::format("x1: {},  x2: {}, x3: {}",x1,x2,x3)<<std::endl; 
  std::cout<<fmt::format("x1 == x2: {}",x1==x2)<<std::endl;
  std::cout<<fmt::format("x1 == x3: {}",x1==x3)<<std::endl;


  u3::U3 w1(7,{2,1});//[4,2,1]
  u3::U3 w2(7,x1);   //[4,3,0]
  u3::U3 w3(10,x1);  //[5,4,1]
  u3::U3 w4(4,2,1);
   
  std::cout<<fmt::format("w1: {} <-> [{},{},{}]",w1,w1.f1(),w1.f2(),w1.f3())<<std::endl;
  std::cout<<fmt::format("w2: {} <-> [{},{},{}]",w2,w2.f1(),w2.f2(),w2.f3())<<std::endl;
  std::cout<<fmt::format("w3: {} <-> [{},{},{}]",w3,w3.f1(),w3.f2(),w3.f3())<<std::endl;
  std::cout<<fmt::format("w3: {} <-> [{},{},{}]",w4,w4.f1(),w4.f2(),w4.f3())<<std::endl;

  std::cout<<fmt::format("w1 == w2: {}",w1==w2)<<std::endl;
  std::cout<<fmt::format("w1 == w3: {}",w1==w3)<<std::endl;
  std::cout<<fmt::format("w2 == w3: {}",w2==w3)<<std::endl;
  std::cout<<fmt::format("w1 == w4: {}",w1==w4)<<std::endl;

  std::cout<<fmt::format("ConguationGrade(x1) = {}",u3::ConjugationGrade(x1))<<std::endl;
  std::cout<<fmt::format("ConguationGrade(x2) = {}",u3::ConjugationGrade(x2))<<std::endl;

  std::cout<<fmt::format("ConguationGrade(w1) = {}",u3::ConjugationGrade(w1))<<std::endl;
  std::cout<<fmt::format("ConguationGrade(w2) = {}",u3::ConjugationGrade(w2))<<std::endl;

  std::cout<<std::endl<<"Functions"<<std::endl<<"----------------------"<<std::endl;
  std::cout<<fmt::format("dim(x1) = {}",dim(x1))<<std::endl;
  std::cout<<fmt::format("dim(w1) = {}",dim(w1))<<std::endl;
  std::cout<<fmt::format("Casimir2(x1) = {:4.2f}",u3::Casimir2(x1))<<std::endl;

  MultiplicityTagged<u3::SU3>::vector su3_product = u3::KroneckerProduct(x1,x2);
  std::cout<<"Kronecker product of x1 x x2: "<<std::endl;
  for(const auto&[x,rho_max] : su3_product)
    std::cout<<fmt::format("x = {}, rho_max = {}",x,rho_max)<<std::endl;

  MultiplicityTagged<u3::U3>::vector u3_product = u3::KroneckerProduct(w1,w2);
  std::cout<<"Kronecker product of w1 x w2: "<<std::endl;
  for(const auto& [w,rho_max] : u3_product)
    std::cout<<fmt::format("w = {}, rho_max = {}",w,rho_max)<<std::endl;

  std::cout<<"Branching x2"<<std::endl;
  MultiplicityTagged<unsigned int>::vector branched_su3 = BranchingSO3(x2);
  for(const auto& [L,kappa_max] : branched_su3)
    std::cout<<fmt::format("L = {}, kappa_max = {}",L,kappa_max)<<std::endl;

} //main
