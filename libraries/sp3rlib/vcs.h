  /****************************************************************
  vcs.h

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  Created by Anna E. McCoy on 3/9/16.

  3/9/16 (aem): Created.
  1/6/17 (aem): Changed phase convention for U3BosonCreation from
                R->L to L->R.
  2/20/17 (aem): Added is_intrinsic option to Kmatrices to allow
    for computation of normalization for intrinsic operators
  7/1/17 (aem): Removed intrinsic option.  Normalization of
    intrinsic operators is now enforced by choice of Nsigma0.
  10/4/17 (aem):
    + Factored Kmatrix and Smatrix calculations
    + Implemented construction of Kmatrix in basis with
      redundent subspaces
    + Removed broken OpenMP version of K matrix calculation
****************************************************************/

#ifndef VCS_H_
#define VCS_H_

#include <Eigen/Eigen>
#include <unordered_map>
#include "basis/operator.h"
#include "mcutils/arithmetic.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3boson.h"

namespace vcs
{
  typedef long double smatrix_float_type;
  typedef basis::OperatorBlock<long double> SMatrixType;
  typedef std::unordered_map<u3::U3,SMatrixType, boost::hash<u3::U3> > SMatrixCache;

  #ifdef HASH_UNIT_TENSOR
  typedef std::unordered_map<u3::U3,Eigen::MatrixXd, boost::hash<u3::U3> > MatrixCache;
  #else
  typedef std::map<u3::U3,Eigen::MatrixXd> MatrixCache;
  #endif



  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map);
  //Calculates the K matrix
  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map, vcs::MatrixCache& Kinv_matrix_map);
  // Generates K matrices and Kinv matrices, for A<6

  
  using KmatrixMap = std::unordered_map<u3::U3, std::array<basis::OperatorBlock<double>, 2>>;

  KmatrixMap GetKMatrices(
    const u3::U3& sigma,
    const vcs::U3BosonSpace& space,
    const double zero_threshold = 1e-12
  );

}  //  namespace

#endif
