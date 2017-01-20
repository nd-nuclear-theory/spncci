/****************************************************************
  null_solver.h

  Null solver for LGI extraction.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/19/17 (aem,mac): Created.
****************************************************************/
#ifndef NULL_SOLVER_H_
#define NULL_SOLVER_H_

#include "eigen3/Eigen/Dense"

#include "mcutils/eigen.h"

namespace lgi
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

}

#endif
