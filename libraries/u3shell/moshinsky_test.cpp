/****************************************************************
  moshinsky_test.cpp
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created to test moshinky.h, moshinsky.cpp.
  5/11/16 (aem,mac): Update namespace.
****************************************************************/

#include "u3shell/moshinsky.h"

int main(int argc, char **argv)
{

  std::cout<< u3shell::MoshinskyCoefficient(2,2,4,0,u3::SU3(4,0))<<std::endl;
  //std::cout<< u3shell::MoshinskyCoefficient(u3::SU3(2,0),u3::SU3(2,0),u3::SU3(4,0),u3::SU3(0,0),u3::SU3(4,0))<<std::endl;

}
