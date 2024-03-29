/****************************************************************
  su3calc.cpp

  Calculate and output SU(3) quantities, e.g., product of two SU(3)
  irreps.

  Syntax and examples:

    product lambda1 mu1 lambda2 mu2

      >>> product 1 1 1 1
      (1,1) x (1,1):
        ((0,0),1)
        ((0,3),1)
        ((1,1),2)
        ((2,2),1)
        ((3,0),1)

    branching lambda mu

      >>> branching 2 2
      (2,2):
        (0,1)
        (2,2)
        (3,1)
        (4,1)

    W lambda1 mu1 kappa1 L1 lambda2 mu2 kappa2 L2 lambda3 mu3 kappa3 L3 rho0

      >>> W 2 0 1 0   0 0 1 0   2 0 1 0   1
      (2,0) 1 0 x (0,0) 1 0 -> (0,0) 1 0 1
  
    U lambda1 mu1 lambda2 mu2 lambda mu lambda3 mu3 lambda12 mu12 rho12 rho12_3 lambda23 mu23 rho23 rho1_23

      >>> U 0 0 2 0 2 0 0 0   2 0 1 1   2 0 1 1
      U[(0,0),(2,0),(2,0),(0,0),(2,0),1,1,(2,0),1,1]
      +1.00000000

      >>> U 2 2 2 0 2 1 2 0   2 0 1 1   0 2 1 1
      U[(2,2),(2,0),(2,1),(2,0),(2,0),1,1,(0,2),1,1]
      +0.36514837

    Z lambda2 mu2 lambda1 mu1 lambda mu lambda3 mu3 lambda12 mu12 rho12 rho12_3 lambda13 mu13 rho13 rho13_2

      >>> Z 0 0 2 0 2 0 0 0   2 0 1 1   2 0 1 1
      Z[(0,0),(2,0),(2,0),(0,0),(2,0),1,1,(2,0),1,1]
      +1.00000000

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  10/17/16 (mac): Create, extending su3_coupler.
  06/25/17 (mac): Update keywords.  Add SO(3) branching.
  11/06/19 (aem): Fix bug.
  03/21/21 (mac): Add Z coefficient.

****************************************************************/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"

////////////////////////////////////////////////////////////////
// control code for each interactive keyword option 
////////////////////////////////////////////////////////////////

void DoHelp()
{
  std::cout << "Syntax:" << std::endl;
  std::cout << "  product lambda1 mu1 lambda2 mu2" << std::endl;
  std::cout << "  branching lambda mu" << std::endl;
  std::cout << "  W lambda1 mu1 kappa1 L1 lambda2 mu2 kappa2 L2 lambda3 mu3 kappa3 L3 rho0" << std::endl;
  std::cout << "  U lambda1 mu1 lambda2 mu2 lambda mu lambda3 mu3 lambda12 mu12 rho12 rho12_3 lambda23 mu23 rho23 rho1_23" << std::endl;
  std::cout << "  Z lambda2 mu2 lambda1 mu1 lambda mu lambda3 mu3 lambda12 mu12 rho12 rho12_3 lambda13 mu13 rho23 rho13_2" << std::endl;
  std::cout << std::endl;
  std::cout << "Enter \"help\" for this help message." << std::endl;
  std::cout << "Enter \"exit\" to exit." << std::endl;
}
void DoKroneckerProduct(int lambda1, int mu1, int lambda2, int mu2)
{
  // label
  u3::SU3 x1(lambda1,mu1); 
  u3::SU3 x2(lambda2,mu2); 
  std::cout << fmt::format("{} x {}:",x1.Str(), x2.Str()) << std::endl;

  // calculate
  MultiplicityTagged<u3::SU3>::vector product = u3::KroneckerProduct(x1,x2);
  for(int i=0; i<product.size(); ++i)  
    std::cout << "  " << product[i].Str() << std::endl;
}

void DoBranchingSO3(int lambda, int mu)
{
  // label
  u3::SU3 x(lambda,mu); 
  std::cout << fmt::format("{}:",x.Str()) << std::endl;

  // calculate
  MultiplicityTagged<int>::vector branching = u3::BranchingSO3(x);
  for(int i=0; i<branching.size(); ++i)  
    std::cout << "  " << branching[i].Str() << std::endl;
}

void DoW(int lambda1, int mu1, int kappa1, int L1, int lambda2, int mu2, int kappa2, int L2, int lambda3, int mu3, int kappa3, int L3, int rho0)
{
  // label
  u3::SU3 x1(lambda1,mu1); 
  u3::SU3 x2(lambda2,mu2); 
  u3::SU3 x3(lambda3,mu3); 
  std::cout
    << fmt::format(
      "{} {} {} x {} {} {} -> {} {} {} {}",
      x1.Str(),kappa1,L1,
      x2.Str(),kappa2,L2,
      x3.Str(),kappa3,L3,
      rho0
    )
    << std::endl;

  // calculate
  double value = u3::W(x1,kappa1,L1,x2,kappa2,L2,x3,kappa3,L3,rho0);
  std::cout << fmt::format("{:+.8f}",value) << std::endl;
}

void DoU(
    int lambda1, int mu1, int lambda2, int mu2, int lambda, int mu, int lambda3, int mu3,
    int lambda12, int mu12, int rho12, int rho12_3,
    int lambda23, int mu23, int rho23, int rho1_23
  )
{
  // label
  u3::SU3 x1(lambda1,mu1); 
  u3::SU3 x2(lambda2,mu2); 
  u3::SU3 x(lambda,mu); 
  u3::SU3 x3(lambda3,mu3); 
  u3::SU3 x12(lambda12,mu12); 
  u3::SU3 x23(lambda23,mu23); 

  std::cout
    << fmt::format(
      "U[{},{},{},{},{},{},{},{},{},{}]",
      x1.Str(),x2.Str(),x.Str(),x3.Str(),x12.Str(),rho12,rho12_3,x23.Str(),rho23,rho1_23
      )
    << std::endl;

  // calculate
  double value = u3::U(x1,x2,x,x3,x12,rho12,rho12_3,x23,rho23,rho1_23);
  std::cout << fmt::format("{:+.8f}",value) << std::endl;
}

void DoZ(
    int lambda2, int mu2, int lambda1, int mu1, int lambda, int mu, int lambda3, int mu3,
    int lambda12, int mu12, int rho12, int rho12_3,
    int lambda13, int mu13, int rho13, int rho13_2
  )
{
  // label
  u3::SU3 x2(lambda2,mu2); 
  u3::SU3 x1(lambda1,mu1); 
  u3::SU3 x(lambda,mu); 
  u3::SU3 x3(lambda3,mu3); 
  u3::SU3 x12(lambda12,mu12); 
  u3::SU3 x13(lambda13,mu13); 

  std::cout
    << fmt::format(
      "Z[{},{},{},{},{},{},{},{},{},{}]",
      x2.Str(),x1.Str(),x.Str(),x3.Str(),x12.Str(),rho12,rho12_3,x13.Str(),rho13,rho13_2
      )
    << std::endl;

  // calculate
  double value = u3::Z(x2,x1,x,x3,x12,rho12,rho12_3,x13,rho13,rho13_2);
  std::cout << fmt::format("{:+.8f}",value) << std::endl;
}

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // initialize su3lib
  u3::U3CoefInit();

  // initial output
  std::cout << "su3calc" << std::endl;
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
      else if (keyword=="branching")
        {
          int lambda, mu;
          
          line_stream >> lambda >> mu; 
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoBranchingSO3(lambda, mu);
        }
      else if (keyword=="product")
        {
          int lambda1, mu1, lambda2, mu2;
          
          line_stream >> lambda1 >> mu1 >> lambda2 >> mu2; 
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoKroneckerProduct(lambda1, mu1, lambda2, mu2);
        }
      else if (keyword=="W")
        {
          int lambda1, mu1, kappa1, L1, lambda2, mu2, kappa2, L2, lambda3, mu3, kappa3, L3, rho0;
          
          line_stream
            >> lambda1 >> mu1 >> kappa1 >> L1
            >> lambda2 >> mu2 >> kappa2 >> L2
            >> lambda3 >> mu3 >> kappa3 >> L3
            >> rho0; 
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoW(
              lambda1, mu1, kappa1, L1,
              lambda2, mu2, kappa2, L2,
              lambda3, mu3, kappa3, L3,
              rho0
            );
        }
      else if (keyword=="U")
        {

          int lambda1, mu1, lambda2, mu2, lambda, mu, lambda3, mu3, lambda12, mu12, rho12, rho12_3, lambda23, mu23, rho23, rho1_23;
          
          line_stream
            >> lambda1 >> mu1 >> lambda2 >> mu2 >> lambda >> mu >> lambda3 >> mu3
            >> lambda12 >> mu12 >> rho12 >> rho12_3
            >> lambda23 >> mu23 >> rho23 >> rho1_23;
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoU(
              lambda1, mu1, lambda2, mu2, lambda, mu, lambda3, mu3,
              lambda12, mu12, rho12, rho12_3,
              lambda23, mu23, rho23, rho1_23
            );
        }
      else if (keyword=="Z")
        {

          int lambda2, mu2, lambda1, mu1, lambda, mu, lambda3, mu3, lambda12, mu12, rho12, rho12_3, lambda13, mu13, rho13, rho13_2;
          
          line_stream
            >> lambda2 >> mu2 >> lambda1 >> mu1 >> lambda >> mu >> lambda3 >> mu3
            >> lambda12 >> mu12 >> rho12 >> rho12_3
            >> lambda13 >> mu13 >> rho13 >> rho13_2;
          mcutils::ParsingCheck(line_stream,line_count,line);

          DoZ(
              lambda2, mu2, lambda1, mu1, lambda, mu, lambda3, mu3,
              lambda12, mu12, rho12, rho12_3,
              lambda13, mu13, rho13, rho13_2
            );
        }
      else
        {
          std::cout << "Unrecognized keyword" << std::endl;
        }

      std::cout << std::endl;
    }

} //main
