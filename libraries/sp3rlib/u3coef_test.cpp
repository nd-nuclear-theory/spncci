/****************************************************************
  u3coef_test.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.
****************************************************************/

#include "am/halfint.h"
#include "utilities/utilities.h"
#include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"


int main(int argc, char **argv)
{

  u3::SU3 x1(2,0);
  u3::SU3 x2(2,0);
  u3::SU3 x12(4,0);
  u3::SU3 x3(2,0);
  u3::SU3 x4(1,0);
  u3::SU3 x23(4,0);
  u3::SU3 x34(3,0);
  u3::SU3 x13(4,0);
  u3::SU3 x24(3,0);
  u3::SU3 x(6,0);
  u3::SU3 xx(7,0);


  u3::U3CoefInit();
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << U(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;
  std::cout << Z(x1, x2, x, x3, x12,1, 1, x23, 1, 1) << std::endl;

  std::cout << W(x1,1,2,x2,1,2,x12,1,4,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,2,1) << std::endl;
  std::cout << W(x1,1,2,x2,1,2,x12,1,0,1) << std::endl;

  std::cout << Unitary9LambdaMu(
    x1,  x2,  x12, 1,
    x3,  x4,  x34, 1,
    x13, x24, xx,  1,
    1,    1,   1
    )<<std::endl;







}
