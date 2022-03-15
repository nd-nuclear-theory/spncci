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
#include "mcutils/eigen.h"
#include "sp3rlib/u3.h"


namespace sp3r{
double Chi(const u3::U3& omega)
{
  double N(omega.N());
  return u3::Casimir2(omega.SU3())+N*N/3.-4*N;
}
double AB00(const u3::U3& sigma, const u3::U3&omega)
  {
    return (Chi(omega)-Chi(sigma))/(2*std::sqrt(6));
  }

using Matrix = basis::OperatorBlock<double>;
}

int main(int argc, char **argv)
{
  u3::U3CoefInit(100);

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////

// Reproducing values from jmp-25-1984-2662-Rowe
  if(true)
  {
    u3::U3 sigma(20,13,10);
    int Nn_max=4;
    sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
    u3::UCoefCache u_coef_cache;

    std::map<std::tuple<u3::U3,u3::U3>,sp3r::Matrix> Amatrix_test ={
      {{{22,13,10},{20,13,10}},sp3r::Matrix{{6.3245}}},
      {{{21,14,10},{20,13,10}},sp3r::Matrix{{-5.5678}}},
      {{{22,15,10},{22,13,10}},sp3r::Matrix{{4.5634},{ 2.07044}}},
      {{{22,15,10},{22,13,10}},sp3r::Matrix{{4.7631},{-6.27253}}}
    };

    for(const auto&[irreps, validation_matrix] : Amatrix_test)
      {
        const auto&[omega_bra,omega_ket] = irreps;
        const auto& subspace_bra=sp3r_space.LookUpSubspace({omega_bra});
        const auto& subspace_ket=sp3r_space.LookUpSubspace({omega_ket});

        sp3r::Matrix
          A_matrix = sp3r::Sp3rRaisingOperator(sigma,subspace_bra,subspace_ket,u_coef_cache);

        assert(mcutils::IsZero(A_matrix-validation_matrix,1e-4));

      }
  }


  // Checking [AxB]^(0,0)
  if(true)
  {
    int Nn_max = 6;
    u3::U3 sigma(16,{2,1});
    u3::UCoefCache u_coef_cache;
    sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
    for(const auto& bra_subspace : sp3r_space)
      for(const auto& ket_subspace : sp3r_space)
        {
          const auto& omega_bra = bra_subspace.omega();
          const auto& omega_ket = ket_subspace.omega();
          if(omega_bra!=omega_ket) continue;

          sp3r::Matrix validation_matrix
            = sp3r::AB00(sigma, omega_bra)*sp3r::Matrix::Identity(bra_subspace.dimension(), ket_subspace.dimension());
          sp3r::Matrix AB_matrix=sp3r::Matrix::Zero(bra_subspace.dimension(), ket_subspace.dimension());
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
    // [B,A]^{(2,2)}=0
    int Nn_max = 8;
    u3::U3 sigma(16,{2,1});
    u3::UCoefCache u_coef_cache;
    sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
    for(const auto& bra_subspace : sp3r_space)
      for(const auto& ket_subspace : sp3r_space)
        {
          const auto& omega_bra = bra_subspace.omega();
          const auto& omega_ket = ket_subspace.omega();
          int rho0_max = u3::OuterMultiplicity(omega_ket,{0,{2,2}},omega_bra);

          if((omega_ket.N()-sigma.N())==Nn_max) continue;
          // if(omega_ket!=omega_bra) continue;

          if(rho0_max==0) continue;


          // fmt::print("{} {}  {}\n",omega_bra,omega_ket,rho0_max);
          sp3r::Matrix commutator_matrix
            =sp3r::Matrix::Zero(bra_subspace.dimension(),ket_subspace.dimension());


          for(const auto& bar_subspace : sp3r_space)
            {
              const auto& omega_bar=bar_subspace.omega();

              if(u3::OuterMultiplicity(omega_ket,{2,{2,0}},omega_bar) && u3::OuterMultiplicity(omega_bar,{-2,{0,2}},omega_bra))
                {
                  commutator_matrix +=
                  u3::UCached(u_coef_cache,omega_ket.SU3(),{2,0},omega_bra.SU3(),{0,2},omega_bar.SU3(),1,1,{2,2},1,1)
                  *sp3r::Sp3rLoweringOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                  *sp3r::Sp3rRaisingOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                }

              if(u3::OuterMultiplicity(omega_ket,{-2,{0,2}},omega_bar) && u3::OuterMultiplicity(omega_bar,{2,{2,0}},omega_bra))
                {
                  commutator_matrix-=
                    u3::UCached(u_coef_cache,omega_ket.SU3(),{0,2},omega_bra.SU3(),{2,0},omega_bar.SU3(),1,1,{2,2},1,1)
                    *sp3r::Sp3rRaisingOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                    *sp3r::Sp3rLoweringOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                }
            }
          assert(mcutils::IsZero(commutator_matrix,1e-6));
          // std::cout<<mcutils::FormatMatrix(commutator_matrix,"+.3f")<<std::endl;
          // std::cout<<mcutils::FormatMatrix(std::sqrt(3./8)*commutator_matrix,"+.3f")<<std::endl;
        }

  }

if(true)
  {
  // Commutation relations given by ap-126-1980-343-Rosensteel
  // [A,B]^{(0,0)}
    u3::U3 sigma(16,{2,1});
    int Nn_max=8;
    sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
    u3::UCoefCache u_coef_cache;
    for(const auto& bra_subspace : sp3r_space)
      for(const auto& ket_subspace : sp3r_space)
        {
          const auto& omega_bra = bra_subspace.omega();
          const auto& omega_ket = ket_subspace.omega();
          if(omega_ket!=omega_bra) continue;
          if((omega_ket.N()-sigma.N())==Nn_max) continue;

          // fmt::print("{} {}\n",omega_bra,omega_ket);
          sp3r::Matrix commutator_matrix
            =sp3r::Matrix::Zero(bra_subspace.dimension(),ket_subspace.dimension());

          sp3r::Matrix validation_matrix
            =std::sqrt(2./3)*TwiceValue(omega_bra.N())*sp3r::Matrix::Identity(bra_subspace.dimension(),ket_subspace.dimension());

          for(const auto& bar_subspace : sp3r_space)
            {
              const auto& omega_bar=bar_subspace.omega();

              if(u3::OuterMultiplicity(omega_ket,{2,{2,0}},omega_bar) && u3::OuterMultiplicity(omega_bar,{-2,{0,2}},omega_bra))
                {
                  commutator_matrix +=
                  u3::UCached(u_coef_cache,omega_ket.SU3(),{2,0},omega_bra.SU3(),{0,2},omega_bar.SU3(),1,1,{0,0},1,1)
                  *sp3r::Sp3rLoweringOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                  *sp3r::Sp3rRaisingOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                }

              if(u3::OuterMultiplicity(omega_ket,{-2,{0,2}},omega_bar) && u3::OuterMultiplicity(omega_bar,{2,{2,0}},omega_bra))
                {
                  commutator_matrix-=
                    u3::UCached(u_coef_cache,omega_ket.SU3(),{0,2},omega_bra.SU3(),{2,0},omega_bar.SU3(),1,1,{0,0},1,1)
                    *sp3r::Sp3rRaisingOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                    *sp3r::Sp3rLoweringOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                }
            }

          // std::cout<<mcutils::FormatMatrix(commutator_matrix,"+.3f")<<std::endl;
          // std::cout<<mcutils::FormatMatrix(validation_matrix,"+.3f")<<std::endl;

          assert(mcutils::IsZero(commutator_matrix-validation_matrix));

        }

  }


  // TODO: Finish setting up test

  std::vector<u3::U3>
  sigma_list = {
      {16,{2,1}}//,
      // {{29,2},{2,1}},
      // {{3,2},{0,0}},
      // {3,{0,0}}
    };
  int Nn_max=4;


  if(false)
  {
    // Commutation relations given by ap-126-1980-343-Rosensteel
    // [A,A]^(2,1)
    for(const auto& sigma : sigma_list)
      {
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
                    =sp3r::Matrix::Zero(bra_subspace.dimension(),ket_subspace.dimension());
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
  }













} //main


