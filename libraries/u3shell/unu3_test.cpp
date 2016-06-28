/****************************************************************
  sp3r_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.

****************************************************************/

#include "sp3rlib/un.h"

#include <iostream>
#include "cppformat/format.h"

int main(int argc, char **argv)
{
  int n=2;
  int A=2;
  MultiplicityTagged<u3::U3S>::vector allowed_irreps;
  un::GenerateAllowedSU3xSU2Irreps(n, A,allowed_irreps);

  for(int i=0; i<allowed_irreps.size(); ++i)
    {
      u3::U3S state=allowed_irreps[i].irrep;
      int multiplicity=allowed_irreps[i].tag;
      std::cout<<fmt::format("{} {}", state.Str(),multiplicity)<<std::endl;
    }

  /// TwoBodyTest
  MultiplicityTagged<u3::U3ST>::vector allowed_two_body_irreps;
  un::GenerateAllowedSU3xSU2xSU2TwoBodyIrreps(n,allowed_two_body_irreps);
  for(int i=0; i<allowed_two_body_irreps.size(); ++i)
  {
    u3::U3ST state=allowed_two_body_irreps[i].irrep;
    int multiplicity=allowed_two_body_irreps[i].tag;
    std::cout<<fmt::format("{} {}", state.Str(),multiplicity)<<std::endl;
  }
}