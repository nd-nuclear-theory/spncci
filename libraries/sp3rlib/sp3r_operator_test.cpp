/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame
  
  SPDX-License-Identifier: MIT 

  3/10/16 (aem,mac): Created.

****************************************************************/
#include "sp3rlib/sp3r_operator.h"

#include <iostream>
#include <functional>

#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3boson.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3.h"


namespace sp3r{

  double Chi(const u3::U3& omega)
  {
    double N(omega.N());
    return u3::Casimir2(omega.SU3())+N*N/3.-4*N;
  }

  // Eigenvalue of [AxB]^{(0,0)}
  double AB00(const u3::U3& sigma, const u3::U3&omega)
    {
      return (Chi(omega)-Chi(sigma))/(2*std::sqrt(6));
    }

  using Matrix = basis::OperatorBlock<double>;


  sp3r::Matrix CommutatorBA00(
    const u3::U3& sigma, // Not needed, just here for consistency with other operators
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace
  )
  {
    const u3::U3& omega_bra=bra_subspace.omega();
    const u3::U3& omega_ket=ket_subspace.omega();
    if(omega_bra!=omega_ket)
      return sp3r::Matrix::Zero(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());

    sp3r::Matrix matrix
      = std::sqrt(2./3)*TwiceValue(omega_bra.N())
          *sp3r::Matrix::Identity(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());

    return matrix;
  }

  sp3r::Matrix CommutatorBA11(
    const u3::U3& sigma,
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace
    )
  {
    return std::sqrt(10./4)*sp3r::SU3Generator(sigma,bra_subspace,ket_subspace);
  }

  sp3r::Matrix CommutatorBA22(
    const u3::U3& sigma,
    const sp3r::U3Subspace& bra_subspace,
    const sp3r::U3Subspace& ket_subspace
    )
  {
    return sp3r::Matrix::Zero(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());
  }

}

int main(int argc, char **argv)
{
  u3::U3CoefInit(100);

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////


  // Reproducing values from jmp-25-1984-2662-Rowe
  if(false)
  {
    u3::U3 sigma(20,13,10);
    int Nn_max=4;
    sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
    u3::UCoefCache u_coef_cache;
    std::map<std::tuple<u3::U3,int>,std::map<std::tuple<u3::U3,u3::U3>,sp3r::Matrix>> Amatrix_test;

    // Reproducing values from jmp-25-1984-2662-Rowe
    Amatrix_test[{{20,13,10},4}]={
          {{{22,13,10},{20,13,10}},sp3r::Matrix{{6.3245}}},
          {{{21,14,10},{20,13,10}},sp3r::Matrix{{-5.5678}}},
          {{{22,15,10},{22,13,10}},sp3r::Matrix{{4.5634},{ 2.07044}}},
          {{{22,15,10},{22,13,10}},sp3r::Matrix{{4.7631},{-6.27253}}}
        };

    //Reproducing values from ap-126-1980-343-Rosensteel
    Amatrix_test[{{48.5_hi,{8,0}},10}]={
          {{{50.5_hi,{10,0}},{48.5_hi,{8,0}}},sp3r::Matrix{{6.55744}}},
          {{{54.5_hi,{12,1}},{52.5_hi,{10,1}}},sp3r::Matrix{{10.24695}}},
          {{{54.5_hi,{12,1}},{52.5_hi,{12,0}}},sp3r::Matrix{{-4.69042}}},
          {{{54.5_hi,{12,1}},{52.5_hi,{12,0}}},sp3r::Matrix{{-4.69042}}}
          // RMEs do not match in cases where upsilon_max!=1.  This may be due to
          // different approaches to calculating sqrt of K^2 matrix.
          // Checked SVD decompositions and singular values match:
          // {{{52.5_hi,{8,2}},{50.5_hi,{6,2}}},sp3r::Matrix{{6.04725},{2.66931}}},
          // Singular values: {6.610177}
          // {{{54.5_hi,{10,2}},{52.5_hi,{8,2}}},sp3r::Matrix{{9.01365,-0.99088},{2.45265,7.94941}}},
          // Singular values:{7.75045, 9.55859}
        };
      Amatrix_test[{{1.5_hi,{0,0}},4}] = {
        {{{3.5_hi,{2,0}},{1.5_hi,{0,0}}},sp3r::Matrix{{1.}}},
        {{{5.5_hi,{4,0}},{3.5_hi,{2,0}}},sp3r::Matrix{{std::sqrt(6)}}}
      };

    for(const auto&[irrep_labels, validation_map] : Amatrix_test)
      {
        const auto&[sigma,Nn_max] = irrep_labels;
        sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
        for(const auto&[subspace_labels, validation_matrix] : validation_map)
          {
            const auto&[omega_bra,omega_ket] = subspace_labels;
            const auto& subspace_bra=sp3r_space.LookUpSubspace({omega_bra});
            const auto& subspace_ket=sp3r_space.LookUpSubspace({omega_ket});

            sp3r::Matrix
              A_matrix = sp3r::Sp3rRaisingOperator(sigma,subspace_bra,subspace_ket,u_coef_cache);

            assert(mcutils::IsZero(A_matrix-validation_matrix,1e-4));

          }
      }
  }




  // Commutation relations given by ap-126-1980-343-Rosensteel
  if(true)
  {
    std::map<u3::U3,std::function<sp3r::Matrix(u3::U3,sp3r::U3Subspace,sp3r::U3Subspace)>>
      commutator_test_map = {
        {{0,{0,0}},sp3r::CommutatorBA00},
        {{0,{1,1}},sp3r::CommutatorBA11},
        {{0,{2,2}},sp3r::CommutatorBA22}
      };

    std::vector<u3::U3> sigma_list = {{16,{2,1}},{{3,2},{0,0}},{5,{2,0}}};

    for(const auto& sigma : sigma_list)
    {
      int Nn_max = 6;
      // u3::U3 sigma(16,{2,1});
      std::cout<<fmt::format("sigma: {}",sigma)<<std::endl;
      u3::UCoefCache u_coef_cache;
      sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
      for(const auto&[w0,validation_function] : commutator_test_map)
        {
          for(const auto& bra_subspace : sp3r_space)
          for(const auto& ket_subspace : sp3r_space)
            {
              const auto& omega_bra = bra_subspace.omega();
              const auto& omega_ket = ket_subspace.omega();

              if(u3::OuterMultiplicity(omega_ket,w0,omega_bra)==0)continue;
              if((omega_ket.N()-sigma.N())==Nn_max) continue;

              sp3r::Matrix commutator_matrix = sp3r::Matrix::Zero(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());
              sp3r::Matrix validation_matrix = validation_function(sigma,bra_subspace,ket_subspace);

              for(const auto& bar_subspace : sp3r_space)
                {
                  const auto& omega_bar=bar_subspace.omega();

                  if(u3::OuterMultiplicity(omega_ket,{2,{2,0}},omega_bar) && u3::OuterMultiplicity(omega_bar,{-2,{0,2}},omega_bra))
                    {
                      commutator_matrix +=
                      u3::UCached(u_coef_cache,omega_ket.SU3(),{2,0},omega_bra.SU3(),{0,2},omega_bar.SU3(),1,1,w0.SU3(),1,1)
                      *sp3r::Sp3rLoweringOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                      *sp3r::Sp3rRaisingOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                    }

                  if(u3::OuterMultiplicity(omega_ket,{-2,{0,2}},omega_bar) && u3::OuterMultiplicity(omega_bar,{2,{2,0}},omega_bra))
                    {
                      commutator_matrix-=
                        u3::UCached(u_coef_cache,omega_ket.SU3(),{0,2},omega_bra.SU3(),{2,0},omega_bar.SU3(),1,1,w0.SU3(),1,1)
                        *sp3r::Sp3rRaisingOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                        *sp3r::Sp3rLoweringOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                    }

                }
              assert(mcutils::IsZero(commutator_matrix-validation_matrix));
            }
        std::cout<<fmt::format("Validation of [B,A] coupled to {} complete",w0)<<std::endl;
        }
    }
  }


  // Checking [AxB]^(0,0)
  if(true)
  {
    int Nn_max = 6;
    // u3::U3 sigma(16,{2,1});
    u3::U3 sigma(5,{2,0});
    u3::UCoefCache u_coef_cache;
    sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
    for(const auto& bra_subspace : sp3r_space)
      for(const auto& ket_subspace : sp3r_space)
        {
          const auto& omega_bra = bra_subspace.omega();
          const auto& omega_ket = ket_subspace.omega();
          if(omega_bra!=omega_ket) continue;

          sp3r::Matrix validation_matrix
            = sp3r::AB00(sigma, omega_bra)*sp3r::Matrix::Identity(bra_subspace.upsilon_max(), ket_subspace.upsilon_max());
          sp3r::Matrix AB_matrix=sp3r::Matrix::Zero(bra_subspace.upsilon_max(), ket_subspace.upsilon_max());
          for(const auto& bar_subspace : sp3r_space)
            {
              const auto& omega_bar = bar_subspace.omega();
              if(u3::OuterMultiplicity(omega_ket,{-2,{0,2}},omega_bar)==0) continue;
              if(u3::OuterMultiplicity(omega_bar,{2,{2,0}},omega_bra)==0) continue;

              AB_matrix
                += u3::UCached(u_coef_cache,omega_ket.SU3(),{0,2},omega_bra.SU3(),{2,0},omega_bar.SU3(),1,1,{0,0},1,1)
                    *sp3r::Sp3rRaisingOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                    *sp3r::Sp3rLoweringOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
            }
          assert(mcutils::IsZero(AB_matrix-validation_matrix,1e-12));
        }
  }

  if(true)
  {
    // Commutation relations given by ap-126-1980-343-Rosensteel
    // [A,A]^(2,1)
    u3::U3 sigma(16,{2,1});
    int Nn_max=8;
    sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);

    u3::UCoefCache u_coef_cache;
    for(const auto& bra_subspace : sp3r_space)
      for(const auto& ket_subspace : sp3r_space)
        {
          u3::SU3 x0(2,1);
          const auto& omega_bra = bra_subspace.omega();
          const auto& omega_ket = ket_subspace.omega();
          if(u3::OuterMultiplicity(omega_ket,{4,x0},omega_bra))
            {

              sp3r::Matrix commutator_matrix
                =sp3r::Matrix::Zero(bra_subspace.upsilon_max(),ket_subspace.upsilon_max());
              for(const auto& bar_subspace : sp3r_space)
                {
                  const auto& omega_bar = bar_subspace.omega();
                  if(!u3::OuterMultiplicity(omega_ket,{2,{2,0}},omega_bar)) continue;
                  if(!u3::OuterMultiplicity(omega_bar,{2,{2,0}},omega_bra)) continue;

                  commutator_matrix
                    -=2*u3::UCached(u_coef_cache,omega_ket.SU3(),{2,0},omega_bra.SU3(),{2,0},omega_bar.SU3(),1,1,x0,1,1)
                      *sp3r::Sp3rRaisingOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                      *sp3r::Sp3rRaisingOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                }
              assert(mcutils::IsZero(commutator_matrix,1e-12));
            }
        }

  }

} //main


