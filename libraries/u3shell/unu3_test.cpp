/****************************************************************
  unu3_test.cpp

  Anna E. McCoy
  University of Notre Dame

  3/10/16 (aem): Created.
****************************************************************/

#include "u3shell/unu3.h"

#include <iostream>
#include "fmt/format.h"

int main(int argc, char **argv)
{
  if(argc<2)
  {
    std::cout<<"Syntax: A shell_num"<<std::endl;
    std::exit(EXIT_FAILURE);
  }
  // single particle cutoff
  // int n=1;
  // int A=2;
  int A=std::stoi(argv[1]);
  int n=std::stoi(argv[2]);

  un::SingleShellAllowedU3SIrreps allowed_irreps;
  un::GenerateAllowedSU3xSU2Irreps(n, A,allowed_irreps);

  for(auto i=allowed_irreps.begin(); i!=allowed_irreps.end(); ++i)
    {
      u3::U3S state=i->first;
      int multiplicity=i->second;
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