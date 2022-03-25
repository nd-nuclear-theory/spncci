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
    bool modify_branching = sp3r::ModifySp3RBranching(sigma);
    sp3r::Sp3RSpace sp3r_space(sigma, Nn_max);
    std::cout<<sp3r_space.DebugStr()<<std::endl;

    // Checking accessors
    assert(sp3r_space.sigma()==sigma);
    assert(sp3r_space.Nn_max()==Nn_max);

    // Checking subspaces
    for(int i=0; i<sp3r_space.size(); ++i)
      {
        const auto& subspace = sp3r_space.GetSubspace(i);
        assert(subspace.u3boson_subspace().dimension()==subspace.K_matrix().cols());
        assert(subspace.omega()==subspace.U3());
        assert(subspace.upsilon_max()<=subspace.dimension());
        if(!modify_branching)
          assert(subspace.upsilon_max()<=subspace.dimension());

        assert(subspace.upsilon_max()==subspace.K_matrix().rows());
        assert(subspace.upsilon_max()==subspace.Kinv_matrix().cols());
        assert(subspace.dimension()==subspace.K_matrix().cols());
        assert(subspace.dimension()==subspace.Kinv_matrix().rows());
        assert(subspace.K_matrix().rows()==subspace.Kinv_matrix().cols());
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

  std::map<u3::U3,sp3r::Sp3RSpace> test_map;
  std::vector<sp3r::Sp3RSpace> test_vector;
  for(const auto& [sigma,check] : modify_basis_test_map)
    {
      fmt::print("{}\n",sigma);
      fmt::print("Is unitary? {}\n",sp3r::IsUnitary(sigma));
      fmt::print("Modify branching? {}\n",sp3r::ModifySp3RBranching(sigma));
      assert(sp3r::ModifySp3RBranching(sigma)==check);
      

      sp3r::Sp3RSpace sp3r_space(sigma, Nn_max);
      std::cout<<sp3r_space.DebugStr()<<std::endl;
      std::cout<<"labels only"<<std::endl;
      bool subspace_labels_only = true;
      sp3r::Sp3RSpace sp3r_space_light(sigma, Nn_max,subspace_labels_only);
      std::cout<<sp3r_space_light.DebugStr()<<std::endl;
      
      const std::tuple<u3::U3> s{};
      test_vector.push_back(sp3r_space);

      // ATTENTION PATRICK:
      // Attempting to add sp3r_space to a map.  This fails with normal basis constructor with
      // sigma stored as std::tuple<u3::U3>.
      test_map[sigma]=sp3r::Sp3RSpace(sigma, Nn_max);

    }

  u3::U3 sigma(16,{2,1});
  std::map<u3::U3,bool> omega0_list = {
    {{0,{1,1}},true},
    {{2,{2,0}},false},
    {{-2,{0,2}},false}
  };

  sp3r::Sp3RSpace sp3r_space(sigma, Nn_max);
  for(const auto& [omega0,su3_generator] : omega0_list)
    {
      sp3r::Sp3RSectors sectors(sp3r_space,omega0,su3_generator);
      std::cout<<sectors.DebugStr()<<std::endl;
    }

} //main


