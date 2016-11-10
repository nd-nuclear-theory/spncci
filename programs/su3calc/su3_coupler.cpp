/****************************************************************
  su3_coupler.cpp

  Calculate and output product of two SU(3) irreps.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/9/16 (aem,mac): Created.
  10/17/16 (mac): Update to read input from stdin.
  Note: End active development. Functionality is absorbed into su3calc.

****************************************************************/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

#include "mcutils/parsing.h"
#include "sp3rlib/u3.h"

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // construction and string conversion
  ////////////////////////////////////////////////////////////////

  std::string line;
  int line_count = 0;
  while (std::getline(std::cin,line) && (line.size()>0))
    {
      ++line_count;

      // parse line
      int lambda1, mu1, lambda2, mu2;
      std::istringstream line_stream(line);
      line_stream >> lambda1 >> mu1 >> lambda2 >> mu2; 
      ParsingCheck(line_stream,line_count,line);

      // label
      u3::SU3 x1(lambda1,mu1); 
      u3::SU3 x2(lambda2,mu2); 
      std::cout << x1.Str() << " x " << x2.Str() << std::endl;

      // calculate
      MultiplicityTagged<u3::SU3>::vector product=KroneckerProduct(x1,x2);
      for(int i=0; i<product.size(); ++i)  
        std::cout << "  " << product[i].Str() << std::endl;
    }

} //main
