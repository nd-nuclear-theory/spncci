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
  // initialize su3lib
  u3::U3CoefInit();


  std::cout<< u3shell::MoshinskyCoefficient(2,2,4,0,u3::SU3(4,0))<<std::endl;
	u3shell::RelativeStateLabelsU3ST bra(2,1,0);
	u3shell::RelativeStateLabelsU3ST ket(2,1,0);
  // testing moshinksy transformation 
  u3shell::RelativeUnitTensorLabelsU3ST tensor(u3::SU3(2,2),1,0,bra,ket);
   MoshinskyTransformation(tensor, 10);

}
