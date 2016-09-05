/****************************************************************
  unit_tensor_expansion.h

  Two-body operator representation and second-quantization.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/9/16 (mac): Created.
****************************************************************/

#ifndef UNIT_TENSOR_EXPANSION_H_
#define UNIT_TENSOR_EXPANSION_H_

#include <map>
#include <unordered_map>  

#include "u3shell/tensor_labels.h"
#include "u3shell/relative_operator.h"
#include "u3shell/two_body_operator.h"

namespace u3shell
{
	// void 
	// NumberOperatorUnitTensorExpansion(
	//   int Nmin, int Nmax, 
	//   u3shell::TwoBodyUnitTensorCoefficientsU3ST& N_operator_expansion, 
	//   int A=2
	// );

	void 
	BrelRelativeUnitTensorExpansion(
		int Nmin, int Nmax,
		u3shell::RelativeUnitTensorCoefficientsU3ST& Brel_operator,
		int A=2
	);
  
  void 
  NrelRelativeUnitTensorExpansion(
  	int Nmin, int Nmax,
  	u3shell::RelativeUnitTensorCoefficientsU3ST& Nrel_operator,
  	int A=2
  );

  // void 
  // NcmUnitTensorExpansion(
  // 	int Nmin, int Nmax,
  // 	u3shell::TwoBodyUnitTensorCoefficientsU3ST& Ncm_operator,
  // 	int A=2
  // );

  void 
  IdentityRelativeUnitTensorExpansion(
  	int Nmin, int Nmax,
  	u3shell::RelativeUnitTensorCoefficientsU3ST& Nrel_operator,
  	int A=2
  );

}
#endif