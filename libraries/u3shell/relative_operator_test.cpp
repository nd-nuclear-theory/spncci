/****************************************************************
  relative_operator_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iostream>

#include "cppformat/format.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/relative_operator.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

void test_relative()
{

      // declare coefficient containers
      u3shell::RelativeUnitTensorCoefficientsU3ST relative_unit_tensor_coefficients;

      // populate operator
      //
      // TO CHECK: if labels make sense
      u3shell::RelativeUnitTensorLabelsU3ST
        relative_unit_tensor_labels(
            u3::SU3(1,1),1,0,  // x0,S0,T0
            u3shell::RelativeStateLabelsU3ST(1,1,1), // bra
            u3shell::RelativeStateLabelsU3ST(1,1,1)  // ket
          );
      relative_unit_tensor_coefficients[relative_unit_tensor_labels] = 1;

      // dump operator
      std::cout << "relative_unit_tensor_coefficients" << std::endl;
      for (auto key_value : relative_unit_tensor_coefficients)
        {

          // extract unit tensor labels and coefficients
          auto labels= key_value.first;
          double coefficient = key_value.second;
        
          std::cout 
            << fmt::format(
                "  {} {:e}",
                labels.Str(),
                coefficient
              )
            << std::endl;
        }

}


int main(int argc, char **argv)
{

  test_relative();

  // termination
  return 0;
}