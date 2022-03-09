/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame
  
  SPDX-License-Identifier: MIT 

  3/10/16 (aem,mac): Created.

****************************************************************/
#include "sp3rlib/sp3r.h"

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

  // Comparison of new with old.
  for(const auto& sigma : sigma_list)
  {
    u3boson::U3BosonSpace u3boson_space(sigma,Nn_max);

    sp3r::Sp3RSpace sp3r_space1(sigma, Nn_max);
    sp3r::Sp3RSpace sp3r_space2(sigma, Nn_max,u3boson_space);
    assert(sp3r_space2.size()==sp3r_space1.size());

    std::cout<<sp3r_space1.DebugStr()<<std::endl;
    std::cout<<"--------------------------------"<<std::endl;
    std::cout<<sp3r_space2.DebugStr()<<std::endl;

    // Checking accessors
    assert(sp3r_space1.sigma()==sp3r_space2.sigma());
    assert(sp3r_space2.sigma()==sigma);
    assert(sp3r_space1.Nn_max()==sp3r_space2.Nn_max());
    assert(sp3r_space2.Nn_max()==Nn_max);

    for(int i=0; i<sp3r_space1.size(); ++i)
      {
        const auto& subspace1 = sp3r_space1.GetSubspace(i);
        const auto& subspace2 = sp3r_space2.GetSubspace(i);
        assert(subspace2.size()==subspace1.size());
        assert(subspace2.U3()==subspace1.U3());
        assert(subspace2.labels()==subspace2.U3());
        assert(subspace1.upsilon_max() == subspace2.upsilon_max());
        assert(subspace2.upsilon_max()<=subspace2.size());

        assert(subspace2.upsilon_max()==subspace2.K_matrix().cols());
        assert(subspace2.upsilon_max()==subspace2.Kinv_matrix().rows());
        assert(subspace2.size()==subspace2.K_matrix().rows());
        assert(subspace2.size()==subspace2.Kinv_matrix().cols());
        assert(subspace2.K_matrix().rows()==subspace2.Kinv_matrix().cols());
      }

  }

  // Check that U(3) is unitary LGI
  std::map<u3::U3,bool> unitary_test_map=
  {
    {{16,{2,1}},true},
    {{3,{0,0}},true},
    {{6,5,1},false}
  };
  
  fmt::print("Is Sp(3,R) irrep unitary?\n");
  for(const auto& [s,check] : unitary_test_map)
    {

      fmt::print("{}: {}\n",s,sp3r::IsUnitary(s));
      assert(sp3r::IsUnitary(s)==check);
    }

  fmt::print("Does branching need to be restricted?\n");
    // Check that U(3) is unitary LGI
  std::map<u3::U3,bool> modify_basis_test_map =
    {
      {{{3,2},{0,0}},true},
      {{3,{0,0}},true},
      {{5,{2,0}},true},
      {{16,{2,1}},false}
    };

  for(const auto& [s,check] : modify_basis_test_map)
    {
      fmt::print("{}\n",s);
      fmt::print("Is unitary? {}\n",sp3r::IsUnitary(s));
      fmt::print("Modify branching? {}\n",sp3r::ModifySp3RBranching(s));
      assert(sp3r::ModifySp3RBranching(s)==check);
      u3boson::U3BosonSpace u3boson_space(s,Nn_max);
      sp3r::Sp3RSpace sp3r_space(s, Nn_max,u3boson_space);
      std::cout<<sp3r_space.DebugStr()<<std::endl;
    }

} //main


