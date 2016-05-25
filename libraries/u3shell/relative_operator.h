/****************************************************************
  relative_operator.h

  Relative operator representation and enumeration.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/25/16 (mac): Created with code from two_body_operator and
    tensor_labels.

****************************************************************/

#ifndef RELATIVE_OPERATOR_H_
#define RELATIVE_OPERATOR_H_

#include "u3shell/tensor_labels.h"

namespace u3shell
{

  ////////////////////////////////////////////////////////////////
  // coefficient storage -- relative
  ////////////////////////////////////////////////////////////////

  typedef
    std::map<u3shell::RelativeUnitTensorLabelsU3ST,double>
    RelativeUnitTensorCoefficientsU3ST;
  // Storage of coefficients of relative unit tensors in U(3)xSxT scheme.



  ////////////////////////////////////////////////////////////////
  // generation of unit tensor label lists
  ////////////////////////////////////////////////////////////////

  // Note: If we define a RelativeSpaceU3ST (or RelativeSectorsU3ST),
  // that would be a more natural argument than Nmax, and it would
  // take care of much of the iteration.  This is the model I am
  // following with the two-body unit tensors.

  // // Generates map of RelativeUnitTensorLabelsU3ST for a given Nmax truncation,   
  // // stored in relative_unit_tensor_labels.
  // // Map keys are N0: number of oscillator quanta carried by the operator
  // // values are vectors of RelativeUnitTensorLabelsU3ST with oscillator quanta N0
  // void GenerateRelativeUnitTensorLabelsU3ST(
  //       int Nmax, 
  //       std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
  //       );

}  // namespace

#endif
