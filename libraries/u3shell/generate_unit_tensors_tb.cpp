/****************************************************************
  two_body_generator.cpp

  Construct files in Tomas's recoupler format for complete set of
  two-body unit tensors, in two-body truncated space.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/25/16 (mac): Created.

****************************************************************/

#include <iostream>

#include "cppformat/format.h"

#include "u3shell/indexing_u3st.h"
#include "u3shell/two_body_operator.h"

////////////////////////////////////////////////////////////////
// two-body unit tensor enumeration
////////////////////////////////////////////////////////////////

std::vector<std::vector<TwoBodyUnitTensorLabelsU3ST>>

  // // Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 
  // // for rp>=r.  To get the other half, use conjugation 

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // initialize su3lib
  u3::U3CoefInit();

  // define relevant two-body space
  const int Nmax = 4;
  u3shell::TwoBodySpaceU3ST two_body_space(Nmax);

  // generate list of unit tensors
  u3shell::TwoBodySectorsU3ST two_body_sectors(two_body_space);
  std::vector<std::vector<u3shell::TwoBodyUnitTensorLabelsU3ST>> two_body_unit_tensor_labels;
  GenerateTwoBodyUnitTensorLabelsU3ST(
      two_body_sectors,
      two_body_unit_tensor_labels
    );
  //       int Nmax, 
  //       std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
  //       );


  // write list of unit tensors

  // generate unit tensors

  // termination
  return 0;
}
