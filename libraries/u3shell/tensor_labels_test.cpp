/****************************************************************
  tensor_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame
  
  SPDX-License-Identifier: MIT

  5/25/16 (aem): Created.
****************************************************************/

#include <iostream>

#include "fmt/format.h"
#include "u3shell/tensor_labels.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // operator/state construction
  if (true)
    {
      u3shell::OperatorLabelsU3ST op(2,u3::SU3(2,0),0,0,0);
      std::cout << op.Str() << std::endl;
      
      u3shell::RelativeStateLabelsU3ST relative_state(1,1,0); // eta,S,T
      std::cout << relative_state.Str() << std::endl;
      
      u3shell::TwoBodyStateLabelsU3ST tb_state(1,1,u3::SU3(2,0),1,0);
      std::cout << tb_state.Str() << std::endl;
    }


  // Coming up with some sensible couplings for the examples:
  //
  // % su3_coupler 2 0 0 2                                                                                                                              
  // (2,0) x (0,2)
  //   ((0,0),1)
  //   ((1,1),1)
  //   ((2,2),1)

  // relative U3ST unit tensor: construction from individual operator labels (DEPRECATED)
  // if (true)
  //   {
  //     std::cout << "relative U3ST: construction from individual operator labels (DEPRECATED)"
  //               << std::endl;
  //     u3shell::RelativeUnitTensorLabelsU3ST
  //       unit_tensor(
  //           u3::SU3(1,1),1,0,
  //           u3shell::RelativeStateLabelsU3ST(2,1,0),
  //           u3shell::RelativeStateLabelsU3ST(2,1,0)
  //         ); // x0,S0,T0,bra,ket 
  //     std::cout << unit_tensor.Str()
  //               << " N0 " << unit_tensor.N0() << std::endl;
  //   }

  // relative U3ST unit tensor
  if (true)
    {
      std::cout << "relative U3ST: construction using base class instance"
                << std::endl;
      u3shell::RelativeUnitTensorLabelsU3ST
        unit_tensor(
            u3shell::OperatorLabelsU3ST(0,u3::SU3(1,1),1,0,0),
            u3shell::RelativeStateLabelsU3ST(2,1,0),
            u3shell::RelativeStateLabelsU3ST(2,1,0)
          );
      std::cout << unit_tensor.Str()
                << " N0 " << unit_tensor.N0() << std::endl;
    }



  // two-body U3ST unit tensor: construction from individual operator labels (DEPRECATED)
  // if (true)
  //   {
  //     std::cout << "two-body U3ST: construction from individual operator labels (DEPRECATED)"
  //               << std::endl;
  //     u3shell::TwoBodyUnitTensorLabelsU3ST
  //       unit_tensor(
  //           u3::SU3(2,2),1,1,
  //           1,
  //           u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0),
  //           u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0)
  //         );
  //     std::cout << unit_tensor.Str()
  //               << " N0 " << unit_tensor.N0() << std::endl;
  //   }

  // two-body U3ST unit tensor
  if (true)
    {
      std::cout << "two-body U3ST: construction using base class instance"
                << std::endl;
      u3shell::TwoBodyUnitTensorLabelsU3ST
        unit_tensor(
            u3shell::OperatorLabelsU3ST(0,u3::SU3(2,2),1,1,0),
            1,
            u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0),
            u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0)
          );
      std::cout << unit_tensor.Str()
                << " N0 " << unit_tensor.N0() << std::endl;
    }

  // two-body U3ST unit tensor: test FAILED assertion
  if (false)
    {
      std::cout << "two-body U3ST: construction using base class instance"
                << std::endl;
      u3shell::TwoBodyUnitTensorLabelsU3ST
        unit_tensor(
            u3shell::OperatorLabelsU3ST(999,u3::SU3(2,2),1,1,0),
            1,
            u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0),
            u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(2,0),1,0)
          );
      std::cout << unit_tensor.Str()
                << " N0 " << unit_tensor.N0() << std::endl;
    }


  // two-body U3S unit tensor: construction from individual operator labels (DEPRECATED)
  // if (true)
  //   {
  //     std::cout << "two-body U3S: construction from individual operator labels (DEPRECATED)"
  //               << std::endl;
  //     u3shell::TwoBodyUnitTensorLabelsU3S
  //       unit_tensor_pn(
  //           u3::SU3(2,2),1,
  //           1,
  //           u3shell::TwoBodyStateLabelsU3S(1,1,u3::SU3(2,0),1),
  //           u3shell::TwoBodyStateLabelsU3S(1,1,u3::SU3(2,0),1)
  //         ); // x0,S0,rho0,bra,ket 
  //     std::cout << unit_tensor_pn.Str()
  //               << " N0 " << unit_tensor_pn.N0() << std::endl;

  //   }

  // two-body U3S unit tensor
  if (true)
    {
      std::cout << "two-body U3S: construction using base class instance"
                << std::endl;
      u3shell::TwoBodyUnitTensorLabelsU3S
        unit_tensor(
            u3shell::OperatorLabelsU3S(0,u3::SU3(2,2),1,0),
            1,
            u3shell::TwoBodyStateLabelsU3S(1,1,u3::SU3(2,0),1),
            u3shell::TwoBodyStateLabelsU3S(1,1,u3::SU3(2,0),1)
          );
      std::cout << unit_tensor.Str()
                << " N0 " << unit_tensor.N0() << std::endl;
    }

 

  // termination
  return 0;
}
