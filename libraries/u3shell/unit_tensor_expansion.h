/****************************************************************
  unit_tensor_expansion.h

  Two-body operator representation and second-quantization.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  4/9/16 (mac): Created.
  9/7/16 (aem): Added Trel operator and renamed Nrel to Nintr
****************************************************************/

#ifndef UNIT_TENSOR_EXPANSION_H_
#define UNIT_TENSOR_EXPANSION_H_

#include <map>
#include <unordered_map>  

#include "u3shell/tensor_labels.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"
#include "moshinsky/moshinsky_xform.h"

namespace u3shell
{

  void 
    BrelRelativeUnitTensorExpansion(
        int Nmin, int Nmax,
        u3shell::RelativeUnitTensorCoefficientsU3ST& Brel_operator,
        int A=2
      );
  // kappa0 and L0 index labels are just dummy labels since we aren't branching

  void ArelRelativeUnitTensorExpansion(
      int Nmin, int Nmax, 
      u3shell::RelativeUnitTensorCoefficientsU3ST& Arel_operator, int A=2
    );


  void RaisingPolynomialRelativeUnitTensorExpansion(
      const u3::U3& n0, int Nmin, int Nmax, 
      u3shell::RelativeUnitTensorCoefficientsU3ST& Prel_operator, int A=2
    );


  void 
    NrelRelativeUnitTensorExpansion(
  	int Nmin, int Nmax,
  	u3shell::RelativeUnitTensorCoefficientsU3ST& Nrel_operator,
  	int A=2
      );

  void 
    IdentityRelativeUnitTensorExpansion(
  	int Nmin, int Nmax,
  	u3shell::RelativeUnitTensorCoefficientsU3ST& Nrel_operator,
  	int A=2
      );

  void 
    TrelRelativeUnitTensorExpansion(
        int Nmin,int Nmax,
        u3shell::RelativeUnitTensorCoefficientsU3ST& Trel_operator,
        int A=2);
  // Generates the unit tensor expansion for the K^2 operator
  // or kinetic energy operator without factor of 1/(2AM),
  // where M is the nucleon mass. 


}
#endif
