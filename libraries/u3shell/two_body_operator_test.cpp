/****************************************************************
  two_body_operator_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <iostream>

#include "cppformat/format.h"

#include "sp3rlib/u3coef.h"
#include "u3shell/two_body_operator.h"

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

void test_moshinsky()
{
  // declare coefficient containers
  u3shell::RelativeUnitTensorCoefficientsU3ST relative_unit_tensor_coefficients;
  u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;


  // populate operator
  u3shell::RelativeUnitTensorLabelsU3ST
    relative_unit_tensor_labels(
        u3::SU3(1,1),1,1,  // x0,S0,T0
        u3shell::RelativeStateLabelsU3ST(1,1,0),  // bra
        u3shell::RelativeStateLabelsU3ST(1,1,0)   //ket
      );
  relative_unit_tensor_coefficients[relative_unit_tensor_labels] = 1;

  // carry out transformations
  // u3shell::TransformRelativeUnitTensorToTwoBodyUnitTensor(
  //     relative_unit_tensor_coefficients,
  //     two_body_unit_tensor_coefficients
  //   );

}

void test_two_body()
{

      // declare coefficient containers
      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
      u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
      u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;

      // populate operator
      u3shell::TwoBodyUnitTensorLabelsU3ST
        two_body_unit_tensor_labels(
            u3::SU3(1,1),1,1,  // x0,S0,T0
            1,  //rho0
            u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(1,1),1,0), // bra
            u3shell::TwoBodyStateLabelsU3ST(1,1,u3::SU3(1,1),1,0)  // ket
          );
      two_body_unit_tensor_coefficients[two_body_unit_tensor_labels] = 1;

      // dump operator
      std::cout << "two_body_unit_tensor_coefficients" << std::endl;
      for (auto key_value : two_body_unit_tensor_coefficients)
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

      // convert to biquads
      u3shell::TransformTwoBodyUnitTensorToBiquad(
          two_body_unit_tensor_coefficients,
          biquad_coefficients
        );

      // dump operator
      std::cout << "biquad_coefficients" << std::endl;
      for (auto key_value : biquad_coefficients)
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

      // convert to pn scheme
      u3shell::TransformBiquadToPNScheme(
          biquad_coefficients,
          biquad_coefficients_pn
        );

      // dump operator
      std::cout << "biquad_coefficients_pn" << std::endl;
      for (auto key_value : biquad_coefficients_pn)
        {

          // extract unit tensor labels and coefficients
          auto labels= key_value.first;
          auto coefficient = key_value.second;
        
          std::cout 
            << fmt::format(
                "  {} {:e} {:e} {:e}",
                labels.Str(),
                coefficient.pppp,
                coefficient.nnnn,
                coefficient.pnnp
              )
            << std::endl;
        }

}


int main(int argc, char **argv)
{

  u3::U3CoefInit();

  test_moshinsky();
  test_two_body();

  // termination
  return 0;
}
