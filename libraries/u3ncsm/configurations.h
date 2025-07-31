/****************************************************************
  configurations.h

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT

  + 01/17/25 (aem): Created.
****************************************************************/
#ifndef CONFIGURATIONS_H_
#define CONFIGURATIONS_H_

#include <vector>
#include <set>

namespace shell
{
  inline int Omega(int N) {return int(2*(N+1)*(N+2));}

  typedef std::vector<int> Configuration;
  typedef std::set<shell::Configuration> ShellConfigurations;


  std::vector<shell::ShellConfigurations>
  generate_configurations(const unsigned int A, const unsigned int Nmax,const unsigned int Nshell_max);


}
#endif