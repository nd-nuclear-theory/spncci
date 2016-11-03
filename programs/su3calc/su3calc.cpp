/****************************************************************
  .cpp

  Calculate and output product of two SU(3) irreps.

  Test input:
    KroneckerProduct 1 1 1 1
    (1,1) x (1,1):
      ((0,0),1)
      ((0,3),1)
      ((1,1),2)
      ((2,2),1)
      ((3,0),1)

    W 2 0 1 0 0 0 1 0 2 0 1 0 1
    (2,0) 1 0 x (0,0) 1 0 -> (0,0) 1 0 1
  
    U 0 0 2 0 2 0 0 0 2 0 1 1 2 0 1 1
    U[(0,0),(2,0),(2,0),(0,0),(2,0),1,1,(2,0),1,1
    +1.00000000

    U 2 2 2 0 2 1 2 0 2 0 1 1 0 2 1 1
    U[(2,2),(2,0),(2,1),(2,0),(2,0),1,1,(0,2),1,1
    +0.36514837

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  10/17/16 (mac): Create, extending su3_coupler.

****************************************************************/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "cppformat/format.h"
#include "utilities/parsing.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"

void DoHelp()
{
  std::cout << "Syntax:" << std::endl;
  std::cout << "  KroneckerProduct lambda1 mu1 lambda2 mu2" << std::endl;
  std::cout << "  W lambda1 mu1 kappa1 L1 lambda2 mu2 kappa2 L2 lambda3 mu3 kappa3 L3 rho0" << std::endl;
  std::cout << "  U lambda1 mu1 lambda2 mu2 lambda mu lambda3 mu3 lambda12 mu12 rho12 rho12_3 lambda23 mu23 rho23 rho1_23" << std::endl;
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
  MultiplicityTagged<u3::SU3>::vector product=KroneckerProduct(x1,x2);
  for(int i=0; i<product.size(); ++i)  
    std::cout << "  " << product[i].Str() << std::endl;
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
  double value = u3::W(x1,kappa1,L2,x2,kappa2,L2,x3,kappa3,L3,rho0);
  std::cout << fmt::format("{:+.8f}",value) << std::endl;
}

void DoU(
    int lambda1, int mu1, int lambda2, int mu2, int lambda, int mu, int lambda3, int mu3,
    int lambda12, int mu12, int rho12, int rho12_3, int lambda23, int mu23, int rho23, int rho1_23
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
      "U[{},{},{},{},{},{},{},{},{},{}",
      x1.Str(),x2.Str(),x.Str(),x3.Str(),x12.Str(),rho12,rho12_3,x23.Str(),rho23,rho1_23
      )
    << std::endl;

  // calculate
  double value = u3::U(x1,x2,x,x3,x12,rho12,rho12_3,x23,rho23,rho1_23);
  std::cout << fmt::format("{:+.8f}",value) << std::endl;
}

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
      if ((keyword=="exit") || (keyword=="bye") || (keyword==":q!"))
        break;

      // skip blank line or hash comment line
      if ((keyword=="") || (keyword=="#"))
        continue;

      // select action based on keyword
      if (keyword=="help")
        {
          DoHelp();
        }
      else if (keyword=="KroneckerProduct")
        {
          int lambda1, mu1, lambda2, mu2;
          
          line_stream >> lambda1 >> mu1 >> lambda2 >> mu2; 
          ParsingCheck(line_stream,line_count,line);

          DoKroneckerProduct(lambda1, mu1, lambda2, mu2);
        }
      else if (keyword=="KroeneckerProduct")
        {
          std::cout << "No, that's not how it's spelled." << std::endl;
        }
      else if (keyword=="W")
        {
          int lambda1, mu1, kappa1, L1, lambda2, mu2, kappa2, L2, lambda3, mu3, kappa3, L3, rho0;
          
          line_stream >> lambda1 >> mu1 >> kappa1 >> L1 >> lambda2 >> mu2 >> kappa2 >> L2 >> lambda3 >> mu3 >> kappa3 >> L3 >> rho0; 
          ParsingCheck(line_stream,line_count,line);

          DoW(lambda1, mu1, kappa1, L1, lambda2, mu2, kappa2, L2, lambda3, mu3, kappa3, L3, rho0);
        }
      else if (keyword=="U")
        {

          int lambda1, mu1, lambda2, mu2, lambda, mu, lambda3, mu3, lambda12, mu12, rho12, rho12_3, lambda23, mu23, rho23, rho1_23;
          
          line_stream >> lambda1 >> mu1 >> lambda2 >> mu2 >> lambda >> mu >> lambda3 >> mu3
                      >> lambda12 >> mu12 >> rho12 >> rho12_3 >> lambda23 >> mu23 >> rho23 >> rho1_23;
          ParsingCheck(line_stream,line_count,line);

          DoU(lambda1, mu1, lambda2, mu2, lambda, mu, lambda3, mu3, lambda12, mu12, rho12, rho12_3, lambda23, mu23, rho23, rho1_23);
        }
      else
        {
          std::cout << "Unrecognized keyword" << std::endl;
        }

      std::cout << std::endl;
    }

} //main
