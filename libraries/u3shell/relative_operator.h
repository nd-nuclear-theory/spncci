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

  std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>
    GenerateRelativeUnitTensorLabelsU3ST(int Nmax);
  // Generate labels for U3ST-scheme relative tensors acting within
  // the relative space of a given Nmax truncation.
  //
  // The resulting unit tensor labels are grouped by N0, i.e., the
  // number of oscillator quanta caried by the operator.  Although N0
  // may in general vary from -2*Nmax to +2*Nmax, only unit tensors
  // with nonnegative N0 (etap>eta) are considered here.
  //
  // Arguments:
  //   Nmax (int) : maximum oscillator truncation
  //
  // Returns:
  //   (std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>)
  //   : map from N0 -> vector of relative unit tensor labels


  // COMMENTS to Anna

  // Comment #1: If we define a RelativeSpaceU3ST (or RelativeSectorsU3ST),
  // that would be a more natural argument than Nmax, and it would
  // take care of much of the iteration.  This is the model I am
  // following with the two-body unit tensors.  Not sure how useful it
  // would be here...

  // Comment #2: Given that the map indices are restricted to
  // nonnegative values 0, 1, 2, ..., and I believe they run
  // sequentially, why not just a do a vector (of vectors)?
  //
  //  std::vector<std::vector<RelativeUnitTensorLabelsU3ST>>
  //
  // Note: Then you just size the vector to Nmax+1 at the beginning,
  // and all the entries are default constructed as empty vectors.  Or
  // you push_back vectors for each No one-by-one.
  //
  // On the other hand, I realize you might want to keep generality to
  // possibly store negative N0.  I am actually doing so for the
  // two-body case, to provide more complete testing.

  // OLD COMMENT from .h file (to delete):
  //
  // // Generates map of RelativeUnitTensorLabelsU3ST for a given Nmax truncation,   
  // // stored in relative_unit_tensor_labels.
  // // Map keys are N0: number of oscillator quanta carried by the operator
  // // values are vectors of RelativeUnitTensorLabelsU3ST with oscillator quanta N0
  //
  // OLD COMMENT from .cpp file (to delete):
  //
  // // Generates a map containing (key, value) pair (N0, operator_labels) of the unit tensors 
  // // for rp>=r.  To get the other half, use conjugation 


  // Generates map of RelativeUnitTensorLabelsU3ST for a given Nmax truncation,   
  // stored in relative_unit_tensor_labels.
  // Map keys are N0: number of oscillator quanta carried by the operator
  // values are vectors of RelativeUnitTensorLabelsU3ST with oscillator quanta N0
  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax, 
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
        );




}  // namespace

#endif
