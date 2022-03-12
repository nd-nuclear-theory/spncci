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
      {16,{2,1}}//,
      // {{29,2},{2,1}},
      // {{3,2},{0,0}},
      // {3,{0,0}}
    };

  for(const auto& sigma : sigma_list)
    {
      sp3r::Sp3RSpace sp3r_space(sigma,Nn_max);
      u3::UCoefCache u_coef_cache;
      for(const auto& bra_subspace : sp3r_space)
        for(const auto& ket_subspace : sp3r_space)
          {
            if(bra_subspace.omega()<=ket_subspace.omega())
              continue;

            std::cout<<fmt::format("({}|| ||{})",
              bra_subspace.omega(), ket_subspace.omega()
            )<<std::endl;

            basis::OperatorBlock<double> A_matrix 
              = Sp3rRaisingOperator(
                  sigma,
                  bra_subspace,
                  ket_subspace,
                  u_coef_cache
                );
            std::cout<<A_matrix<<std::endl;
          }

    }

} //main


