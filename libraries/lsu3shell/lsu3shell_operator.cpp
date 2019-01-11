/****************************************************************
  lsu3shell_operator.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lsu3shell/lsu3shell_operator.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <omp.h>

#include "fmt/format.h"
#include "moshinsky/moshinsky_xform.h"


namespace lsu3shell
{

  void GenerateModelSpaceFile(
      const std::string& model_space_filename, int Z, int N, int Nmax, int parity
    )
  {
    // If parity is -1, then full model space with no parity restriction,
    // otherwise construct model space by parity
    int Nmin=(parity==-1)?0:parity;
    int Nstep=(parity==-1)?1:2;

    std::ofstream model_stream(model_space_filename);
    // Model space with no restirction on J
    model_stream<<Z<<"  "<<N<<"  "<<-1<<std::endl;
    // Get full space up to Nmax in steps of Nex=2
    for(int Nex=Nmin; Nex<=Nmax; Nex+=Nstep)
      model_stream<<Nex<<"  "<<-1<<std::endl;
    
    model_stream.close();
  }

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

    u3shell::TransformRelativeTensorToTwobodyTensor(
        relative_tensor_expansion,
        twobody_space,
        two_body_unit_tensor_coefficients,
        "NAS"   
      );
    // std::cout<<"TransformRelativeTensorToTwobodyTensor"<<std::endl;
    // for(auto it=two_body_unit_tensor_coefficients.begin(); it!=two_body_unit_tensor_coefficients.end(); ++it)
    // {
    //   std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
    // }
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

    // container for computed biquad expansion of unit tensor 
    std::vector<u3shell::TwoBodyUnitTensorCoefficientsU3SPN> 
      biquad_coefficients_pn_list(relative_tensor_labels.size());

    #pragma omp parallel for schedule(runtime)  
    for(int i=0; i<relative_tensor_labels.size(); ++i)
      {
        const u3shell::RelativeUnitTensorLabelsU3ST& tensor=relative_tensor_labels[i];
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
        
        biquad_coefficients_pn_list[i]=biquad_coefficients_pn;
      }
      
    // writing out recoupler files
      for(int i=0; i<relative_tensor_labels.size(); ++i)
      {
        u3shell::TwoBodyUnitTensorCoefficientsU3SPN& biquad_coefficients_pn=biquad_coefficients_pn_list[i];
        std::string operator_stream_filename = fmt::format("relative_unit_{:06d}.recoupler",i);
        std::ofstream operator_stream(operator_stream_filename);
        WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
        operator_stream.close();
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
