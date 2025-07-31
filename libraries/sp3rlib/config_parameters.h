/****************************************************************
  config_parameters.h

  Configuration parameters for sp3rlib calculations

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/


#ifndef SP3R_CONFIG_PARAMS_H_
#define SP3R_CONFIG_PARAMS_H_

#include <Eigen/Core>

#include "sp3rlib/sp3r.h"


namespace sp3r
{

typedef double MatrixFloatType;
typedef Eigen::Matrix<MatrixFloatType,Eigen::Dynamic,1> Vector;   
typedef Eigen::Matrix<MatrixFloatType,Eigen::Dynamic,Eigen::Dynamic> Matrix;

}

#endif