/****************************************************************
  moshinsky.h
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created to test moshinky.h, moshinsky.cpp
****************************************************************/

#include "am/halfint.h"
#include "u3.h"
#include "moshinsky.h"

int main(int argc, char **argv)
{

	std::cout<< u3::MoshinskyCoefficient(2,2,4,0,u3::SU3(4,0))<<std::endl;
	std::cout<< u3::MoshinskyCoefficient(u3::SU3(2,0),u3::SU3(2,0),u3::SU3(4,0),u3::SU3(0,0),u3::SU3(4,0))<<std::endl;

}