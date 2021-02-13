/****************************************************************
  sp3rcalc.cpp

  Calculate and output Sp(3,R) quantities, e.g., irrep contents.

  Syntax and examples:

    raising Nnmax

    >>> raising 4
    Nnmax = 4
    0(0,0)
    2(2,0)
    4(0,2)
    4(4,0)

    irrep N lambda mu Nnmax

    >>> irrep 16 2 1 4
    sigma = 16(2,1), Nnmax = 4
    # omega n rho_n
      16(2,1)      0(0,0)         1
      18(0,3)      2(2,0)         1
      18(1,1)      2(2,0)         1
      18(2,2)      2(2,0)         1
      18(3,0)      2(2,0)         1
      18(4,1)      2(2,0)         1
      20(0,1)      4(0,2)         1
      20(1,2)      4(0,2)         1
      20(1,2)      4(4,0)         1
      20(2,0)      4(0,2)         1
      20(2,3)      4(0,2)         1
      20(2,3)      4(4,0)         1
      20(3,1)      4(0,2)         1
      20(3,1)      4(4,0)         1
      20(4,2)      4(4,0)         1
      20(5,0)      4(4,0)         1
      20(6,1)      4(4,0)         1

    subspaces N lambda mu Nnmax

    >>> subspaces 16 2 1 4
    sigma = 16(2,1), Nnmax = 4
    # omega.N omega.lambda omega.mu   dim
       16.0   2   1      1
       18.0   0   3      1
       18.0   1   1      1
       18.0   2   2      1
       18.0   3   0      1
       18.0   4   1      1
       20.0   0   1      1
       20.0   1   2      2
       20.0   2   0      1
       20.0   2   3      2
       20.0   3   1      2
       20.0   4   2      1
       20.0   5   0      1
       20.0   6   1      1

    lattice N lambda mu Nnmax

    >>> lattice 16 2 1 4

    # omega1.N omega1.lambda omega1.mu   omega2.N omega2.lambda omega2.mu
       16.0   2   1    18.0   0   3
       16.0   2   1    18.0   1   1
       16.0   2   1    18.0   2   2
       16.0   2   1    18.0   3   0
       16.0   2   1    18.0   4   1
       18.0   0   3    20.0   0   1
       18.0   0   3    20.0   1   2
       18.0   0   3    20.0   2   3
       18.0   1   1    20.0   0   1
       18.0   1   1    20.0   1   2
       18.0   1   1    20.0   2   0
       18.0   1   1    20.0   3   1
       18.0   2   2    20.0   1   2
       18.0   2   2    20.0   2   0
       18.0   2   2    20.0   2   3
       18.0   2   2    20.0   3   1
       18.0   2   2    20.0   4   2
       18.0   3   0    20.0   1   2
       18.0   3   0    20.0   3   1
       18.0   3   0    20.0   5   0
       18.0   4   1    20.0   2   3
       18.0   4   1    20.0   3   1
       18.0   4   1    20.0   4   2
       18.0   4   1    20.0   5   0
       18.0   4   1    20.0   6   1

  Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  6/25/17 (mac): Created, with framework from su3calc.

****************************************************************/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "sp3rlib/sp3r.h"

////////////////////////////////////////////////////////////////
// control code for each interactive keyword option 
////////////////////////////////////////////////////////////////

void DoHelp()
{
  std::cout << "Syntax:" << std::endl;
  std::cout << "  raising Nnmax" << std::endl;
  std::cout << "  irrep N lambda mu Nnmax" << std::endl;
  std::cout << "  subspaces N lambda mu Nnmax" << std::endl;
  std::cout << "  lattice N lambda mu Nnmax" << std::endl;
  std::cout << std::endl;
  std::cout << "Enter \"help\" for this help message." << std::endl;
  std::cout << "Enter \"exit\" to exit." << std::endl;
}

void DoRaising(int Nnmax)
{
  std::cout << fmt::format("Nnmax = {}",Nnmax) << std::endl;
  std::vector<u3::U3> raising_polynomial_labels = sp3r::RaisingPolynomialLabels(Nnmax);
  for (const u3::U3& n : raising_polynomial_labels)
    {
      std::cout << n.Str() << std::endl;
    }
}

void DoIrrep(int N, int lambda, int mu, int Nnmax)
{

  // Note this essentially provides the same output as the irrep's
  // DebugStr, but in neater tabulation.

  u3::U3 sigma(N,u3::SU3(lambda,mu));
  std::cout << fmt::format("sigma = {}, Nnmax = {}",sigma.Str(),Nnmax) << std::endl;
  std::cout << "# omega n rho_n" << std::endl;

  sp3r::Sp3RSpace irrep(sigma,Nnmax);
  for (int subspace_index=0; subspace_index<irrep.size(); ++subspace_index)
    {
      const sp3r::U3Subspace& subspace = irrep.GetSubspace(subspace_index);
      const u3::U3 omega = subspace.U3();
      for (int state_index=0; state_index<subspace.size(); ++state_index)
        {
          MultiplicityTagged<u3::U3> n_rho = subspace.GetStateLabels(state_index);
          const u3::U3 n = n_rho.irrep;
          int rho = n_rho.tag;
          std::cout << fmt::format("  {:12s} {:12s} {:3d}",omega.Str(),n.Str(),rho) << std::endl;
        }
    }
}

void DoSubspaces(int N, int lambda, int mu, int Nnmax)
{
  u3::U3 sigma(N,u3::SU3(lambda,mu));
  std::cout << fmt::format("sigma = {}, Nnmax = {}",sigma.Str(),Nnmax) << std::endl;

  sp3r::Sp3RSpace irrep(sigma,Nnmax);

  std::cout << "# omega.N omega.lambda omega.mu   dim" << std::endl;
  for (int subspace_index=0; subspace_index<irrep.size(); ++subspace_index)
    {
      const sp3r::U3Subspace& subspace = irrep.GetSubspace(subspace_index);
      const u3::U3 omega = subspace.U3();
      int dimension = subspace.size();
      std::cout 
        << fmt::format("  {:5.1f} {:3d} {:3d}   {:4d}",float(omega.N()),omega.SU3().lambda(),omega.SU3().mu(),dimension)
        << std::endl;
    }

}

void DoLattice(int N, int lambda, int mu, int Nnmax)
{
  u3::U3 sigma(N,u3::SU3(lambda,mu));
  std::cout << fmt::format("sigma = {}, Nnmax = {}",sigma.Str(),Nnmax) << std::endl;

  sp3r::Sp3RSpace irrep(sigma,Nnmax);

  std::cout << "# omega1.N omega1.lambda omega1.mu   omega2.N omega2.lambda omega2.mu" << std::endl;
  for (int subspace1_index=0; subspace1_index<irrep.size(); ++subspace1_index)
    for (int subspace2_index=0; subspace2_index<irrep.size(); ++subspace2_index)
      {
        const sp3r::U3Subspace& subspace1 = irrep.GetSubspace(subspace1_index);
        const sp3r::U3Subspace& subspace2 = irrep.GetSubspace(subspace2_index);
        const u3::U3 omega1 = subspace1.U3();
        const u3::U3 omega2 = subspace2.U3();

        // check U(1)xSU(3) connection by +2(2,0)
        bool connected =  (
            (omega2.N() == omega1.N()+2)
            && u3::OuterMultiplicity(omega1.SU3(),u3::SU3(2,0),omega2.SU3())
          );

        if (connected)
          {
            std::cout 
              << fmt::format(
                  "  {:5.1f} {:3d} {:3d}   {:5.1f} {:3d} {:3d}",
                  float(omega1.N()),omega1.SU3().lambda(),omega1.SU3().mu(),
                  float(omega2.N()),omega2.SU3().lambda(),omega2.SU3().mu()
                )
              << std::endl;
          }
      }
}

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // initialize su3lib
  // u3::U3CoefInit();

  // initial output
  std::cout << "sp3rcalc" << std::endl;
  std::cout << std::endl;
  DoHelp();
  std::cout << std::endl;


  // process user input lines
  std::string line;
  int line_count = 0;
  while (std::getline(std::cin,line))
    // while (std::getline(std::cin,line) && (line.size()>0))
    {
      ++line_count;

      // std::cout << line << std::endl;

      // set up for line parsing
      std::istringstream line_stream(line);
      std::string keyword;
      line_stream >> keyword;

      // terminate on "exit"
      if ((keyword=="exit") || (keyword=="exit()") || (keyword=="bye")
          || (keyword=="q") || (keyword==":q!"))
        break;

      // skip blank line or hash comment line
      if ((keyword=="") || (keyword=="#"))
        continue;

      // select action based on keyword
      if (keyword=="help")
        {
          DoHelp();
        }
      else if (keyword=="raising")
        {
          int Nnmax;
          
          line_stream >> Nnmax; 
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoRaising(Nnmax);
        }
      else if (keyword=="irrep")
        {
          int N, lambda, mu, Nnmax;
          
          line_stream >> N >> lambda >> mu >> Nnmax; 
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoIrrep(N,lambda,mu,Nnmax);
        }
      else if (keyword=="subspaces")
        {
          int N, lambda, mu, Nnmax;
          
          line_stream >> N >> lambda >> mu >> Nnmax; 
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoSubspaces(N,lambda,mu,Nnmax);
        }
      else if (keyword=="lattice")
        {
          int N, lambda, mu, Nnmax;
          
          line_stream >> N >> lambda >> mu >> Nnmax; 
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoLattice(N,lambda,mu,Nnmax);
        }
      else
        {
          std::cout << "Unrecognized keyword" << std::endl;
        }

      std::cout << std::endl;
    }

} //main
