/****************************************************************
  sp3rsj_test.cpp

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  10/17/23 (cvc): Created.
  06/04/24 (cvc): Modified tests.

****************************************************************/
#include "sp3rlib/sp3r.h"
#include "sp3rlib/sp3rsj.h"

#include <iostream>
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "fmt/format.h"
#include "am/halfint.h"
#include "am/am.h"



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

  HalfInt st1 = HalfInt(3,2);
  HalfInt st2 = HalfInt(4,2);
  HalfInt::vector s_list = am::ProductAngularMomenta(st1,st2);
  HalfInt::vector j_list = am::ProductAngularMomenta(st1,st2);
  //HalfInt::vector s_list = am::ProductAngularMomenta(HalfInt(0,2),HalfInt(0,2));
  //HalfInt::vector j_list = am::ProductAngularMomenta(HalfInt(12,2),HalfInt(0,2));

  /**/
  // Testing
  for(const auto& sigma : sigma_list)
    for(const auto& s_val : s_list)
      for(const auto& j_val : j_list)
      {
        //bool modify_branching = sp3r::ModifySp3RBranching(sigma);
        sp3r::Sp3RSJSpace sp3rsj_space(sigma,Nn_max,s_val,j_val);

        // Test Sp3RSJSpace
        std::cout<<"Testing Sp3RSJSpace"<<std::endl;
        std::cout<<sp3rsj_space.DebugStr()<<std::endl;

        // Checking accessors
        assert(sp3rsj_space.sigma()==sigma);
        assert(sp3rsj_space.Nn_max()==Nn_max);
        assert(sp3rsj_space.S()==s_val);
        assert(sp3rsj_space.J()==j_val);

        // Test Sp3RSpace
        std::cout<<"Testing Sp3RSJSpace.Sp3RSpace"<<std::endl;
        sp3r::Sp3RSpace sp3r_space_test = sp3rsj_space.Sp3RSpace();
        std::cout<<sp3r_space_test.DebugStr()<<std::endl;
        
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

  sp3r::Sp3RSJSpace sp3rsj_space(sigma, Nn_max);
  for(const auto& [omega0,su3_generator] : omega0_list)
    {
      sp3r::Sp3RSJSectors sectors(sp3rsj_space,omega0,su3_generator);
      std::cout<<sectors.DebugStr()<<std::endl;
    }
  */


  ////////////////////////////////////////////////////////////////
  // Sp(3,R) Operator Space construction test
  ////////////////////////////////////////////////////////////////
  // u3::U3 sigma = u3::U3(16,u3::SU3(2,1));
  std::vector<u3::U3>
  op_list = {
      u3::U3(2,{2u,0u}),
      u3::U3(-2,{0u,2u}),
      u3::U3(0,{0u,0u}),
      //u3::U3(0,{1u,1u}),
      u3::U3(0,{2u,2u}),
      u3::U3(4,{4u,0u}),
      u3::U3(-4,{0u,4u})
    };

  /**/
  // Testing
  for(const auto& s_val : s_list)
    for(const auto& j_val : j_list)
    {
      //bool modify_branching = sp3r::ModifySp3RBranching(sigma);
      sp3r::OperatorSpaceSp3RSJ op_space(s_val,j_val,op_list);
      
      std::cout<<op_space.DebugStr()<<std::endl;

      // Checking accessors
      assert(op_space.S0()==s_val);
      assert(op_space.J0()==j_val);

    }
  
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) Observable construction test
  ////////////////////////////////////////////////////////////////
  /*
  sp3r::OperatorSpaceSp3RSJ op_space(st1,st2,op_list);

  // Construct identity
  sp3r::OperatorInfoSp3RSJ id_info(op_space,&IdentitySp3R);
  // Construct raising operator
  sp3r::OperatorInfoSp3RSJ A_info(op_space,&RaisingOperatorSp3R);
  // Construct lowering operator
  sp3r::OperatorInfoSp3RSJ B_info(op_space,&LoweringOperatorSp3R);
  // Construct harmonic oscillator
  sp3r::OperatorInfoSp3RSJ H_info(op_space,&HarmonicOscillatorSp3R);

  u3::U3 sigma(16,{2,1});
  //u3::U3 sigma(5,{2,0});
  //int Nn_max=4;
  u3::UCoefCache u_coef_cache;
  sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
  int outer_multiplicity = 1;
  u3::U3 omega_Id(0,{0,0});
  u3::U3 omega0(0,{0,0});
  u3::U3 omega1(0,{1,1});
  u3::U3 omega2(0,{2,2});
  // check these labels
  u3::U3 omega_A(2,{2,0});
  u3::U3 omega_B(-2,{0,2});
  u3::U3 omega_C(0,{1,1});

  for (const auto& bra_subspace : sp3r_space)
  {
    for (const auto& ket_subspace : sp3r_space)
    {
      const auto& bra_labels = bra_subspace.omega();
      const auto& ket_labels = ket_subspace.omega();
      // Skip highest N values in tests
      if (((bra_labels.N() - sigma.N()) == Nn_max) || ((ket_labels.N() - sigma.N()) == Nn_max)) continue;
      std::cout<<fmt::format("Testing I")<<std::endl;
      basis::OperatorBlock<double> Imat = id_info.ComputeRMESector(omega_Id,outer_multiplicity,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      std::cout<<Imat<<std::endl;
      std::cout<<fmt::format("Testing A")<<std::endl;
      basis::OperatorBlock<double> Amat = A_info.ComputeRMESector(omega_A,outer_multiplicity,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      std::cout<<Amat<<std::endl;
      std::cout<<fmt::format("Testing B")<<std::endl;
      basis::OperatorBlock<double> Bmat = B_info.ComputeRMESector(omega_B,outer_multiplicity,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      std::cout<<Bmat<<std::endl;
      std::cout<<fmt::format("Testing H")<<std::endl;
      basis::OperatorBlock<double> Hmat = H_info.ComputeRMESector(omega_Id,outer_multiplicity,sp3r_space,bra_subspace,ket_subspace,u_coef_cache);
      std::cout<<Hmat<<std::endl;
    }
  }
  */



} //main


