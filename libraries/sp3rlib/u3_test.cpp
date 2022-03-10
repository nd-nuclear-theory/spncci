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

  // std::array<HalfInt,3> test = {HalfInt(5,2),HalfInt(3,2),HalfInt(1,2)};
  std::map<u3::U3,std::array<HalfInt,3>>
  u3_cartesian_check = {
    {{7,{2,1}},{4,2,1}},
    {{7,{1,3}},{4,3,0}},
    {{10,{1,3}},{5,4,1}},
    {{4,2,1},{4,2,1}},
    {{{5,2},{3,2},{1,2}},{HalfInt(5,2),HalfInt(3,2),HalfInt(1,2)}},
    {{0,{1,1}},{1,0,-1}}
  };
  // u3::U3 w1(7,{2,1});//[4,2,1]
  // u3_vector.push_back(w1);
  // u3::U3 w2(7,x1);   //[4,3,0]
  // u3_vector.push_back(w2);
  // u3::U3 w3(10,x1);  //[5,4,1]
  // u3_vector.push_back({10,x1});
  // u3::U3 w4(4,2,1);
  // u3_vector.push_back(w4);
  // u3::U3 w5(HalfInt(5,2),HalfInt(3,2),HalfInt(1,2));
  // u3_vector.push_back(w5);
  // u3::U3 w6(HalfInt(11,2),{2,1});
  // u3_vector.push_back(w6);




  for(const auto& [w,f] : u3_cartesian_check)
    {
      const auto&[f1,f2,f3] = f;
      assert((f1==w.f1()) && (f2==w.f2()) && (f3==w.f3()));
      fmt::print("w = {:8}  <-> f = [{},{},{}]\n",w.Str(),f1,f2,f3);
    }


  std::map<u3::U3,std::tuple<int,int,double>>
  u3_function_check = {
    {{7,{2,1}},  {3,15,10.666667}},
    {{10,{1,3}}, {4,24,16.666667}},
    {{0,{1,1}},  {2, 8, 6.000000}},
    {{{5,2},{3,2},{1,2}}, {2,8,6.000000}}
  };

  for(const auto& [w,test_values] : u3_function_check)
    {
      std::cout<<fmt::format("ConguationGrade[{}] = {}",w,u3::ConjugationGrade(w))<<std::endl;
      std::cout<<fmt::format("dim[{}] = {}",w,dim(w))<<std::endl;
      std::cout<<fmt::format("Casimir2[{}] = {:4.6f}",w,u3::Casimir2(w.SU3()))<<std::endl;

      const auto&[g,dimension,casimir] = test_values;
      assert(u3::ConjugationGrade(w)==g);
      assert(u3::dim(w)==dimension);
      assert(fabs(u3::Casimir2(w.SU3())-casimir)<1e-5);
    }



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

  u3::U3 w0(0,{1,1});
  u3::U3 w(5,4,2);
  std::cout<<u3::OuterMultiplicity({1,1},{1,2},{1,2})<<std::endl;
  for(const auto& [wp,rho] : u3::KroneckerProduct(w,w0))
    {
      std::cout<<fmt::format("{} : {}",wp,rho)<<std::endl;
    }

  // {
  // fmt::print("Branching {}",x2);
  // MultiplicityTagged<unsigned int>::vector branched_su3 = u3::BranchingSO3(x2);
  // for(const auto& [L,kappa_max] : branched_su3)
  //   std::cout<<fmt::format("L = {}, kappa_max = {}",L,kappa_max)<<std::endl;
  // }
} //main
