/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.

****************************************************************/

#include "sp3rlib/sp3r.h"

#include <iostream>

int main(int argc, char **argv)
{


  ////////////////////////////////////////////////////////////////
  // Sp(3,R) irrep construction test
  ////////////////////////////////////////////////////////////////

  // 16(2,1)  Nn_max=6

  u3::U3 sigma = u3::U3(16,u3::SU3(2,1));
  int Nn_max = 6;

  // raising polynomials
  std::vector<u3::U3> polynomial_labels = sp3r::RaisingPolynomialLabels(Nn_max);

  for (auto it = polynomial_labels.begin(); it != polynomial_labels.end(); ++it)
    std::cout << " label " << (*it).Str() << std::endl;
    
  sp3r::Sp3RSpace irrep(sigma,Nn_max);

  for (int i_subspace=0; i_subspace<irrep.size(); ++i_subspace)
    {
      const sp3r::U3Subspace& subspace = irrep.GetSubspace(i_subspace);
      u3::U3 w = subspace.GetSubspaceLabels();
      std::cout << "subspace " << w.Str() << std::endl;

      for (int i_state=0; i_state<subspace.Dimension(); ++i_state)
	{
	  MultiplicityTagged<u3::U3> n_tagged = subspace.GetStateLabels(i_state);
	  std::cout << "  " << i_state << " " << n_tagged.Str() << std::endl;
	}
    }

} //main
