/****************************************************************
  tensor_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iostream>

#include "cppformat/format.h"
#include "u3shell/tensor.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  if (true)
    {

      u3shell::OperatorLabelsU3ST op(2,u3::SU3(2,0),0,0,0);
      std::cout << op.Str() << std::endl;

      u3shell::RelativeStateLabelsU3ST relative_state(1,1,0); // eta,S,T
      std::cout << relative_state.Str() << std::endl;

      u3shell::TwoBodyStateLabelsU3ST tb_state(1,1,u3::SU3(2,0),1,0);
      std::cout << tb_state.Str() << std::endl;

      // Anna -- Do these label values make sense?
      u3shell::RelativeUnitTensorLabelsU3ST
        unit_tensor_relative(
                             u3::SU3(1,1),1,1,
                             u3shell::RelativeStateLabelsU3ST(1,1,0),
                             u3shell::RelativeStateLabelsU3ST(1,1,0)
                             ); // x0,S0,T0,bra,ket 
      std::cout << unit_tensor_relative.Str() << " N0 " << unit_tensor_relative.N0() << std::endl;

      // Anna -- Do these label values make sense?
      u3shell::TwoBodyUnitTensorLabelsU3ST
        unit_tensor_two_body(
                             u3::SU3(4,0),1,1,
                             1,
                             u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0),
                             u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0)
                             ); // x0,S0,T0,rho0,bra,ket 
      std::cout << unit_tensor_two_body.Str() << " N0 " << unit_tensor_two_body.N0() << std::endl;


    }
  

  // termination
  return 0;
}
