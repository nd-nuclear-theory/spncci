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
  fmt::print("Checking SU(3) construction, equality and print statements\n");
  u3::SU3 x1(1,3);
  u3::SU3 x2(2,2);
  u3::SU3 x3(1,3);
  std::vector<u3::SU3> su3_vector{x1,x2,x3};

  fmt::print("Checking accessors for {}\n",x1);
  fmt::print("lambda   = {}  mu   = {}\n",x1.lambda(),x1.mu());
  fmt::print("lambda_u = {}  mu_u = {}\n",x1.lambda_u(),x1.mu_u());
  assert(x1.lambda()==x1.lambda_u());
  assert(x1.mu()==x1.mu_u());

  fmt::print("{} {}\n",x1.Str(),x1);
  assert(x1.Str()==fmt::format("{}",x1));

  fmt::print("Checking equality\n");
  for(const auto& x : su3_vector)
    for(const auto& xp : su3_vector)
      {
        std::cout<<fmt::format("{:6} == {:6}: {}",x.Str(),xp.Str(),x==xp)<<std::endl;
        bool equals_querry = x.lambda()==xp.lambda();
        equals_querry &= x.mu()==xp.mu();
        assert(equals_querry == (x==xp));
      }
  
  fmt::print("u3::dim{}:  {}\n",x1.Str(),u3::dim(x1));
  assert(u3::dim(x1) == 24);

  fmt::print("u3::ConjugationGrade{} = {}\n",x1.Str(),u3::ConjugationGrade(x1));
  assert(u3::ConjugationGrade(x1)==4);
  fmt::print("u3::ConjugationGrade{} = {}\n",x2.Str(),u3::ConjugationGrade(x2));
  assert(u3::ConjugationGrade(x2)==4);

  fmt::print("u3::Conjugate{} = {}\n",x1.Str(),u3::Conjugate(x1));
  assert(u3::Conjugate(x1)==u3::SU3(3,1));  
  fmt::print("u3::Conjugate{} = {}\n",x2.Str(),u3::Conjugate(x2));
  assert(u3::Conjugate(x2)==u3::SU3(2,2));

  fmt::print("\nChecking U(3) construction\n");
  std::vector<u3::U3> u3_vector;
  u3::U3 w1(7,{2,1});//[4,2,1]
  u3_vector.push_back(w1);
  u3::U3 w2(7,x1);   //[4,3,0]
  u3_vector.push_back(w2);
  u3::U3 w3(10,x1);  //[5,4,1]
  u3_vector.push_back({10,x1});
  u3::U3 w4(4,2,1);
  u3_vector.push_back(w4);
  u3::U3 w5(HalfInt(5,2),HalfInt(3,2),HalfInt(1,2));
  u3_vector.push_back(w5);
  u3::U3 w6(HalfInt(11,2),{2,1});
  u3_vector.push_back(w6);

  


  for(const auto& w : u3_vector)
  {
    HalfInt f1,f2,f3;
    std::tie(f1,f2,f3) = w.f(); 
    fmt::print("w = {:10}  <-> f = [{},{},{}]\n",w.Str(),f1,f2,f3);
    assert((f1==w.f1()) && (f2==w.f2()) && (f3==w.f3()));
  }


  u3::U3 test(2,{2,0});
  assert(test == u3::U3(2,0,0));
  for(const auto& w : u3_vector)
    {
      fmt::print("{} == {}: {}\n",w1,w,w1==w);
    }


  std::cout<<fmt::format("ConguationGrade({}) = {}",w1.Str(),u3::ConjugationGrade(w1))<<std::endl;
  std::cout<<fmt::format("ConguationGrade({}) = {}",w2.Str(),u3::ConjugationGrade(w2))<<std::endl;

  std::cout<<std::endl<<"Functions"<<std::endl<<"----------------------"<<std::endl;
  std::cout<<fmt::format("dim(x1) = {}",dim(x1))<<std::endl;
  std::cout<<fmt::format("dim(w1) = {}",dim(w1))<<std::endl;
  std::cout<<fmt::format("Casimir2(x1) = {:4.2f}",u3::Casimir2(x1))<<std::endl;

  
  u3::SU3 x(3,3), xp(2,1);
  std::map<u3::SU3,unsigned int> products_from_su3lib{
      {u3::SU3(6,2),1},
      {u3::SU3(5,4),1},
      {u3::SU3(4,3),2},
      {u3::SU3(3,5),1},
      {u3::SU3(2,4),2},
      {u3::SU3(1,6),1},
      {u3::SU3(0,5),1},
      {u3::SU3(5,1),1},
      {u3::SU3(3,2),2},
      {u3::SU3(1,3),1},
      {u3::SU3(4,0),1},
      {u3::SU3(2,1),1}
    };
  fmt::print("Kronecker product of {} x {}\n",x,xp);
  auto test_products = u3::KroneckerProduct(x,xp);
  assert(test_products.size()==products_from_su3lib.size());
  for(const auto& [xxp,rho] : test_products)
    {
      fmt::print("x = {} rho = {}\n",xxp,rho);
      assert(products_from_su3lib[xxp]==rho);
    }

  {
  MultiplicityTagged<u3::U3>::vector u3_product = u3::KroneckerProduct(w1,w2);
  std::cout<<"Kronecker product of w1 x w2: "<<std::endl;
  fmt::print("Possible SU3 products: \n");
  MultiplicityTagged<u3::SU3>::vector su3_product = u3::KroneckerProduct(w1.SU3(),w2.SU3());
  for(const auto& [x,rho_max] : su3_product)
    fmt::print("x = {}  rho_max = {}\n",x,rho_max);
  std::cout<<"----------------"<<std::endl;
  for(const auto& [w,rho_max] : u3_product)
    std::cout<<fmt::format("  w = {}, rho_max = {}",w,rho_max)<<std::endl;
  }
  {
  MultiplicityTagged<u3::U3>::vector u3_product = u3::KroneckerProduct(w4,w5);
  fmt::print("Kronecker product of {} x {}: \n", w4,w5);
  for(const auto& [w,rho_max] : u3_product)
    std::cout<<fmt::format("  w = {} <-> [{},{},{}], rho_max = {}",w,w.f1(),w.f2(),w.f3(),rho_max)<<std::endl;
  }

  {
  MultiplicityTagged<u3::U3>::vector u3_product = u3::KroneckerProduct(w5,w6);
  fmt::print("Kronecker product of {} x {}: \n", w5,w6);
  for(const auto& [w,rho_max] : u3_product)
    std::cout<<fmt::format("  w = {} <-> [{},{},{}], rho_max = {}",w,w.f1(),w.f2(),w.f3(),rho_max)<<std::endl;
  }


  {
  fmt::print("Branching {}",x2);
  MultiplicityTagged<unsigned int>::vector branched_su3 = u3::BranchingSO3(x2);
  for(const auto& [L,kappa_max] : branched_su3)
    std::cout<<fmt::format("L = {}, kappa_max = {}",L,kappa_max)<<std::endl;
  }
} //main
