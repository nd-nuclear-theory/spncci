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

  // examine raising polynomials
  //
  // We generate and print them here for testing purposes only.
  // Normally sp3r::RaisingPolynomialLabels is called directly by the
  // Sp3RSpace constructor.

  std::vector<u3::U3> polynomial_labels = sp3r::RaisingPolynomialLabels(Nn_max);
  for (auto it = polynomial_labels.begin(); it != polynomial_labels.end(); ++it)
    std::cout << " label " << (*it).Str() << std::endl;
    
  // examine Sp3RSpace object
  sp3r::Sp3RSpace irrep(sigma,Nn_max);
  std::cout << irrep.DebugString();
  //for (int i_subspace=0; i_subspace<irrep.size(); ++i_subspace)
  //  {
  //    // set up reference to subspace of interest
  //    //
  //    // Using a reference avoids copying the U3Subspace object (and
  //    // all its lookup tables).
  //    const sp3r::U3Subspace& subspace = irrep.GetSubspace(i_subspace);
  //
  //    // print subspace labels
  //    u3::U3 omega = subspace.GetSubspaceLabels();
  //    std::cout << "subspace " << omega.Str() << std::endl;
  //
  //    // enumerate state labels within subspace
  //    for (int i_state=0; i_state<subspace.Dimension(); ++i_state)
  //      {
  //        MultiplicityTagged<u3::U3> n_tagged = subspace.GetStateLabels(i_state);
  //        std::cout << "  " << i_state << " " << n_tagged.Str() << std::endl;
  //      }
  //  }

} //main
