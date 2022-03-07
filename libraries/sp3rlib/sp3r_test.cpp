/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame
  
  SPDX-License-Identifier: MIT 

  3/10/16 (aem,mac): Created.

****************************************************************/

#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"
#include "fmt/format.h"
#include <iostream>
#include "sp3rlib/u3coef.h"
int main(int argc, char **argv)
{
  u3::U3CoefInit(39);

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////
  u3::U3 sigma = u3::U3(16,u3::SU3(2,1));
  int Nn_max = 6;

  // examine raising polynomials
  //
  // We generate and print them here for testing purposes only.
  // Normally sp3r::RaisingPolynomialLabels is called directly by the
  // Sp3RSpace constructor.
  std::vector<u3::U3> polynomial_labels = sp3r::RaisingPolynomialLabels(Nn_max);
  for (const auto&n : polynomial_labels)
    fmt::print("{}\n",n);

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
    }

  for(const auto& [s,dummy] : modify_basis_test_map)
    {

      sp3r::Sp3RSpace irrep(s,Nn_max);
      std::cout << irrep.DebugStr();
      std::cout<<std::endl;
    }

  std::map<u3::U3,bool> helium3_set = {
    {{3,{0,0}},true},
    {{5,{2,0}},true},
    {{5,{0,1}},true},
    {{7,{4,0}},true},
    {{7,{2,1}},false},
    {{7,{1,0}},true},
    {{7,{0,2}},false},
    {{9,{6,0}},true},
    {{9,{3,0}},true},
    {{9,{4,1}},false},
    {{9,{2,2}},false},
    {{9,{1,1}},true}
  };

  for(const auto& [s,is_unitary] : helium3_set)
    {
      fmt::print("\n{}\n",s);
      fmt::print("Is unitary? {}\n",sp3r::IsUnitary(s));
      assert(sp3r::IsUnitary(s)==is_unitary);
      fmt::print("Modify branching? {}\n",sp3r::ModifySp3RBranching(s));
      if(sp3r::IsUnitary(s))
        {
          sp3r::Sp3RSpace irrep(s,Nn_max);
          std::cout << irrep.DebugStr();
        }

    }

} //main


