/****************************************************************
  relative_operator.h

  Relative operator representation and enumeration.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/25/16 (mac): Created with code from two_body_operator and
                 tensor_labels.
  12/6/16 (aem): Added optional parameters to GenerateRelativeUnitTensors
****************************************************************/

#ifndef RELATIVE_OPERATOR_H_
#define RELATIVE_OPERATOR_H_

#include "sp3rlib/sp3r.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/u3st_scheme.h"
#include "u3shell/tensor_labels.h"
// #include "u3shell/moshinsky.h"

namespace u3shell
{

  ////////////////////////////////////////////////////////////////////////
  // coefficient container for expansion in terms of relative unit tensors
  ////////////////////////////////////////////////////////////////////////

  // typedef
  //   std::map<u3shell::RelativeUnitTensorLabelsU3ST,double>
  //   RelativeUnitTensorCoefficientsU3ST;

  typedef
    std::unordered_map<
      u3shell::RelativeUnitTensorLabelsU3ST,
      double,
      boost::hash<u3shell::RelativeUnitTensorLabelsU3ST>
    > RelativeUnitTensorCoefficientsU3ST;


  ////////////////////////////////////////////////////////////////
  // generation of unit tensor label lists
  ///////////////////////////////////////////////////////////////

  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax,
        int N1v, 
        std::vector<RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels,
        int J0=-1,
        int T00=-1,
        bool restrict_positive_N0=false
        );

  // Generate labels for U3ST-scheme relative tensors acting within
  // the relative space of a given Nmax truncation, J0 and T0 and stores
  // label in vector
  // 
  //  Tensors in vector are ordered by: 
  //    N0, Sp,Tp,S,T,S0,T0,etap (eta=etap-N0)
  //   
  //  Tensors are subject to trianglarity constraints on 
  //    (Sp,S0,S) and (Tp,T0,T)
  //
  //  and parity constraint on bra and ket
  //    (eta+S+T)~1 and (etap+Sp+Tp)~1
  //
  //  Arguments:
  //    space (input) : space on which the unit tensors are represented
  //    relative_unit_tensor_labels (output) : container for unit tensor labels grouped by N0
  //    J0 (optional) : Angular momentum of operator. If default value of -1, operators are not
  //                    restricted to those which can branch to J0. 
  //    T0 (optional) : Isospin componenet of operator.  If default value of -1, then T0 takes
  //                    all values that satisfy triangularity constraint. 
  //    bool (optional) : if true restricts, N0 or operator to positive values only


  void GenerateRelativeUnitTensorLabelsU3ST(
        int Nmax,
        int N1v, 
        std::map<int,std::vector<RelativeUnitTensorLabelsU3ST>>& relative_unit_tensor_labels,
        int J0=-1,
        int T00=-1,
        bool restrict_positive_N0=false
      );
  // Overload of above function where container is map with unit tensors sorted by N0. 


  double RelativeNumberOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

  double RelativeSp3rRaisingOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

  double RelativeSp3rLoweringOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

  double RelativeKineticEnergyOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

  double RelativeMassQuadrupoleOperator(const u3shell::RelativeStateLabelsU3ST& bra, const u3shell::RelativeStateLabelsU3ST& ket);

}  // namespace

#endif
