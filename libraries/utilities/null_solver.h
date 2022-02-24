/****************************************************************
  null_solver.h

  Null solver for LGI extraction.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  1/19/17 (aem,mac): Created.
  10/10/19 (aem): Updated null solver to handel case when matrix
    is zero matrix
  2/23/22 (aem): Move from lgi library to utilities and changed namespace
****************************************************************/
#ifndef NULL_SOLVER_H_
#define NULL_SOLVER_H_

#include <Eigen/Dense>

#include "mcutils/eigen.h"

namespace utils
{

  void FindNullSpaceSVD(Eigen::MatrixXd& A, Eigen::MatrixXd& V, double threshold, bool verbose=false);
  // Extract null vectors by SVD.
  //
  // Assumes A is a "tall" matrix (rows>=cols).
  //
  // Note: Nullity is null_vectors.cols().  If nullity is zero,
  // returned matrix will be (dimensions)x0.
  //
  // Arguments:
  //   A (matrix, input): ratrix representation of operator
  //   null_vectors (matrix, output): matrix with null vectors as columns
  //   threshold (double): threshold for nonzero singular values


  Eigen::MatrixXd FindNullSpaceSVD(Eigen::MatrixXd& A,  int nullity, double threshold=1e-10);

}

#endif
