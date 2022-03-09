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
  3/9/22 (aem): Extracted U(3)-boson functions to u3boson and
    added new Kmatrix generator which uses new U3BosonSpace.
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
  typedef std::unordered_map<u3::U3,SMatrixType > SMatrixCache;

  #ifdef HASH_UNIT_TENSOR
  typedef std::unordered_map<u3::U3,Eigen::MatrixXd > MatrixCache;
  #else
  typedef std::map<u3::U3,Eigen::MatrixXd> MatrixCache;
  #endif

  ////////////////////////////////////////////////////////////////////////////////////////////////

  inline double Omega(const u3::U3& n, const u3::U3& omega)

  // Calculate Omega factor used in Kmatrix calculations.
  //
  // Based on protopye vcs.py and equation given in
  //   D. J. Rowe, J. Math Phys. 25 (1984) 2662.
  //
  // Returns:
  //   (double) : Omega factor
  {
    const auto& [n1,n2,n3] = std::tuple<int,int,int>(n.f());
    const auto& [w1,w2,w3] = omega.f();

    double value=0;
    value += double(int(2*w1)*w1-n1*n1+8*(w1-n1)-2*(2*w1-n1));
    value += double(int(2*w2)*w2-n2*n2+8*(w2-n2)-4*(2*w2-n2));
    value += double(int(2*w3)*w3-n3*n3+8*(w3-n3)-6*(2*w3-n3));
    return value/4.;
  }

  // DEPRECATED
  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map);
  //Calculates the K matrix
  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map, vcs::MatrixCache& Kinv_matrix_map);
  // Generates K matrices and Kinv matrices

  
  using KmatrixMap = std::unordered_map<u3::U3, std::array<basis::OperatorBlock<double>, 2>>;

  KmatrixMap GenerateKmatrices(
    const u3::U3& sigma,
    const u3boson::U3BosonSpace& space,
    const double zero_threshold = 1e-12
  );
  // Calculate the VCS K matrix and its inverse.
  //
  // Input:
  //  sigma : LGI of the Sp(3,R) irrep
  //  space : Corresponding U(3)-boson space
  //  zero_threshold(optional) : Threshold for determining if a U(3) irrep should be included in the basis
  //    For more detailed discusion, see D. J. Rowe, A. E. McCoy and M. A. Caprio. Phys. Script. 91 (2016) 033003.
  //
  // Output
  //  Returns map keyed by omega.  Value is 2D array containing K and Kinv.
  //    K is a u3bosn_dimension x upsilon_max matrix
  //    Kinv is a upsilon_max x u3boson_dimension matrix

}  //  namespace

#endif
