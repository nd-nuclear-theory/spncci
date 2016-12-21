/****************************************************************
  lsu3shell_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lsu3shell/lsu3shell_operator.h"

#include <fstream>
#include <iostream>
#include <algorithm>

#include "cppformat/format.h"
#include "moshinsky/moshinsky_xform.h"


namespace lsu3shell
{

  void 
  GenerateLSU3ShellOperator(
      int Nmax, 
      const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_tensor_expansion,
      std::string filename,
      bool un_u3_restrict
    )
  {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    // declare coefficient containers
    u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
    u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
    //moshinsky transform and accumulate coefficients
    u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
// CHECK THIS FIRST
    u3shell::TransformRelativeTensorToTwobodyTensor(
        relative_tensor_expansion,
        twobody_space,
        two_body_unit_tensor_coefficients,
        "NAS"   
      );
    std::cout<<"TransformRelativeTensorToTwobodyTensor"<<std::endl;
    for(auto it=two_body_unit_tensor_coefficients.begin(); it!=two_body_unit_tensor_coefficients.end(); ++it)
    {
      std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
    }
    // convert to biquads
    u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
    // convert to pn scheme
    u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn,un_u3_restrict);
    
    std::ofstream operator_stream(filename);
    WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
    operator_stream.close();
  }

  void
  GenerateLSU3ShellOperator(
      int Nmax, 
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_tensor_labels,
      bool un_u3_restrict
    )
  {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    int i=0;
    for(auto tensor : relative_tensor_labels)
      {
        // declare coefficient containers
        u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
        u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
        u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
        
        // Moshinsky transform unit tensors to two-body operators 
        MoshinskyTransformTensor(tensor, twobody_space, 
          two_body_unit_tensor_coefficients, "NAS");
        
        // convert to biquads
        u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
        // convert biquads to pn scheme
        u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn,un_u3_restrict);
        std::string operator_stream_filename = fmt::format("relative_unit_{:06d}.recoupler",i);
        std::ofstream operator_stream(operator_stream_filename);
        WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
        operator_stream.close();
        ++i;
      }
  }

  void GenerateLSU3ShellOperator(
      int Nmax, 
      const u3shell::TwoBodyUnitTensorCoefficientsU3ST& twobody_tensor_expansion,
      int operator_index,
      bool un_u3_restrict
      )
  {
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
    // declare coefficient containers
    u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
    u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
    u3shell::TransformTwoBodyUnitTensorToBiquad(twobody_tensor_expansion,biquad_coefficients);
    // convert to pn scheme
    u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn,un_u3_restrict);
    std::string operator_stream_filename =fmt::format("two_body_unit_{:06d}.recoupler",operator_index);
    std::ofstream operator_stream(operator_stream_filename);
    WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
    operator_stream.close();
  }

}// end namespace
