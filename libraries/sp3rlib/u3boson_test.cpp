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
#include "mcutils/eigen.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3boson.h"

namespace u3boson
{
  using Matrix = basis::OperatorBlock<double>;
}

int main(int argc, char **argv)
{
  u3::U3CoefInit(39);

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////
  u3::U3 sigma = u3::U3(16,u3::SU3(2u,1u));
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

  u3::U3 omega(18,{4u,1u}), omegap(20,{4u,2u});
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

  double comparison_rme
    =ParitySign(u3::ConjugationGrade(omegap)+u3::ConjugationGrade(omega))
      *std::sqrt(double(u3::dim(omegap.SU3()))/double(u3::dim(omega.SU3())))
      *u3boson::U3BosonCreationRME(sigma,np,1,omegap,sigma,n,1,omega);


  fmt::print(
    "{}  {}\n",
    u3boson::U3BosonAnnihilationRME(sigma,n,1,omega,sigma,np,1,omegap),
    comparison_rme
  );




  // Testing default constructor
  std::map<u3::U3,u3boson::U3BosonSpace> test_map;
  test_map[sigma]=u3boson_space;
  // =u3boson::U3BosonSpace(sigma,Nn_max);




// Checking commutation relations [a,a^\dagger]^(0,0)
  if(true)
  {
    int Nn_max = 6;
    u3::U3 sigma(16,{2u,1u});
    // u3::U3 sigma(5,{2u,0u});
    u3::UCoefCache u_coef_cache;
    u3boson::U3BosonSpace boson_space(sigma,Nn_max);
    for(const auto& bra_subspace : boson_space)
      for(const auto& ket_subspace : boson_space)
        {
          const auto& omega_bra = bra_subspace.omega();
          const auto& omega_ket = ket_subspace.omega();
          if(omega_bra!=omega_ket) continue;
          if((omega_bra.N()-sigma.N())==Nn_max) continue;

          u3boson::Matrix validation_matrix
            = std::sqrt(6)
              *u3boson::Matrix::Identity(bra_subspace.dimension(), ket_subspace.dimension());

          u3boson::Matrix operator_block
            = u3boson::Matrix::Zero(bra_subspace.dimension(),ket_subspace.dimension());

          for(const auto& bar_subspace : boson_space)
            {
              const auto& omega_bar = bar_subspace.omega();
              // a x a^\dagger term
              if(u3::OuterMultiplicity(omega_bar,{-2,{0u,2u}},omega_bra)>0)
                if(u3::OuterMultiplicity(omega_ket,{2,{2u,0u}},omega_bar)>0)
                  {
                    operator_block
                      += //std::sqrt(u3::dim(omega_bar.SU3())/(6.0*u3::dim(omega_bra.SU3())))
                      u3::UCached(u_coef_cache,omega_ket.SU3(),{2u,0u},omega_bra.SU3(),{0u,2u},omega_bar.SU3(),1,1,{0u,0u},1,1)
                          *u3boson::U3BosonAnnihilationOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                          *u3boson::U3BosonCreationOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                  }

              // a^\dagger x a term
              if(u3::OuterMultiplicity(omega_ket,{-2,{0u,2u}},omega_bar)>0)
                if(u3::OuterMultiplicity(omega_bar,{2,{2u,0u}},omega_bra)>0)
                  {
                    operator_block
                      -= //std::sqrt(u3::dim(omega_bar.SU3())/(6.0*u3::dim(omega_bra.SU3())))
                          u3::UCached(u_coef_cache,omega_ket.SU3(),{0u,2u},omega_bra.SU3(),{2u,0u},omega_bar.SU3(),1,1,{0u,0u},1,1)
                          *u3boson::U3BosonCreationOperator(sigma,bra_subspace,bar_subspace,u_coef_cache)
                          *u3boson::U3BosonAnnihilationOperator(sigma,bar_subspace,ket_subspace,u_coef_cache);
                  }

            }
          assert(mcutils::IsZero(operator_block-validation_matrix,1e-8));
          // std::cout<<"---------"<<std::endl;
          // std::cout<<fmt::format("{}  {} ",omega_bra,omega_ket)<<std::endl;
          // std::cout<<mcutils::FormatMatrix(validation_matrix,".2f")<<std::endl<<std::endl;
          // std::cout<<mcutils::FormatMatrix(operator_block,".2f")<<std::endl;

        }
  }

} //main


