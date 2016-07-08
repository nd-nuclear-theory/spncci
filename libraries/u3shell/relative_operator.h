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
#include "u3shell/u3st_scheme.h"
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


  void GenerateRelativeUnitTensorLabelsU3ST(
        const u3shell::RelativeSpaceU3ST& space,
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels
        );


  std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>
    GenerateRelativeUnitTensorLabelsU3ST(int Nmax);
    // Depreciated version

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
  //Depreciated version

  // Overload function if we don't need operators sorted by N0
  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax, 
        std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
        );


  double RelativeNumberOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

  double RelativeSp3rRaisingOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

  double RelativeSp3rLoweringOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

  double RelativeKineticEnergyOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);


  // define a pointer to a function 


}  // namespace

#endif
