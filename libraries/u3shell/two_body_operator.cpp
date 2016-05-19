/****************************************************************
  two_body_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "u3shell/two_body_operator.h"

namespace u3shell {

  void TransformTwoBodyUnitTensorToBiquad
    (
     const u3shell::TwoBodyUnitTensorCoefficientsU3ST& unit_tensor_coefficients,
     u3shell::TwoBodyUnitTensorCoefficientsU3ST& biquad_coefficients
     )
  {
    for (auto key_value : unit_tensor_coefficients)
      {

        // extract unit tensor labels and coefficients
        u3shell::TwoBodyUnitTensorLabelsU3ST unit_tensor_labels = key_value.first;
        double coefficient = key_value.second;

        
      }
  }


 
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
