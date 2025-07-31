/****************************************************************
  sp3rl_test.cpp
  
  SPDX-License-Identifier: MIT 

  10/16/23 (cvc): Created.

****************************************************************/
#include "sp3rlib/sp3r.h"
#include "sp3rlib/sp3rl.h"

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
  unsigned int Nn_max = 6;

  std::vector<u3::U3>
  sigma_list = {
      u3::U3(16,{2u,1u}),
      u3::U3({29,2},{2u,1u})
    };
  std::vector<unsigned int> l_list = {0u,1u,2u,3u,4u,5u,6u};

  /**/
  // Testing
  for(const auto& sigma : sigma_list)
    for(const auto& l_val : l_list)
    {
      //bool modify_branching = sp3r::ModifySp3RBranching(sigma);
      sp3r::Sp3RLSpace sp3rl_space(sigma,Nn_max,l_val);
      
      std::cout<<sp3rl_space.DebugStr()<<std::endl;

      // Checking accessors
      assert(sp3rl_space.sigma()==sigma);
      assert(sp3rl_space.Nn_max()==Nn_max);
      assert(sp3rl_space.L()==l_val);
      
      /*

      // Checking subspaces
      for(std::size_t i=0; i<sp3rl_space.size(); ++i)
        {
          const auto& subspace = sp3rl_space.GetSubspace(i);
          assert(subspace.nonorthogonal_basis().dimension()==subspace.K_matrix().cols());
          assert(subspace.omega()==subspace.U3());
          assert(subspace.upsilon_max()<=subspace.dimension());
          if(!modify_branching)
            assert(subspace.upsilon_max()<=subspace.dimension());

          assert(subspace.upsilon_max()==subspace.K_matrix().rows());
          assert(subspace.upsilon_max()==subspace.Kinv_matrix().cols());
          assert(subspace.dimension()==subspace.K_matrix().cols());
          assert(subspace.dimension()==subspace.Kinv_matrix().rows());
          assert(subspace.K_matrix().rows()==subspace.Kinv_matrix().cols());

          // Test restriction of U3Subspace
          sp3r::U3Subspace(
            subspace.omega(),
            subspace.upsilon_max(),
            subspace.nonorthogonal_basis_ptr(),
            subspace.K_matrix(),
            subspace.Kinv_matrix(),
            true,
            {l_val,l_val}
          );


          sp3r::U3Subspace(
            subspace.omega(),
            subspace.upsilon_max(),
            subspace.nonorthogonal_basis_ptr(),
            subspace.K_matrix(),
            subspace.Kinv_matrix(),
            false
          );

        }*/

    }

  /*
  u3::U3 sigma(16,{2u,1u});
  std::map<u3::U3,bool> omega0_list = {
    {{0,{1u,1u}},true},
    {{2,{2u,0u}},false},
    {{-2,{0u,2u}},false}
  };

  sp3r::Sp3RLSpace sp3rl_space(sigma, Nn_max);
  for(const auto& [omega0,su3_generator] : omega0_list)
    {
      sp3r::Sp3RLSectors sectors(sp3rl_space,omega0,su3_generator);
      std::cout<<sectors.DebugStr()<<std::endl;
    }
  */

} //main


