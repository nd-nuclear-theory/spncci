/****************************************************************
  su3_coupler.cpp

  Calculate and output product of two SU(3) irreps.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/9/16 (aem,mac): Created.

****************************************************************/

#include <cstdlib>
#include "sp3rlib/u3.h"

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // construction and string conversion
  ////////////////////////////////////////////////////////////////

  // quick and dirty...  to neaten!
  int lambda1 = atoi(argv[1]);
  int mu1 = atoi(argv[2]);
  int lambda2 = atoi(argv[3]);
  int mu2 = atoi(argv[4]);

  u3::SU3 x1(lambda1,mu1); 
  u3::SU3 x2(lambda2,mu2); 
  std::cout << x1.Str() << " x " << x2.Str() << std::endl;

  MultiplicityTagged<u3::SU3>::vector product=KroneckerProduct(x1,x2);
  for(int i=0; i<product.size(); ++i)  
    std::cout << "  " << product[i].Str() << std::endl;

} //main
