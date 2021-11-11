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

  double BosonCreationRME(const u3::U3& np, const u3::U3& n);
  // SU(3) Reduced matrix element of a^\dagger boson creation operator
  //
  // Based on protoype u3boson.py  Formula is given by:
  //   G. Rosensteel and D. J. Rowe. J. Math Phys. 24 (1983) 2461.
  //
  // Returns:
  //    rme: (double) reduced matrix element of boson creation operator.

  // double SMatrix(const u3::U3& s, const u3::U3& omega, MultiplicityTagged<u3::U3>& n1_tagged, MultiplicityTagged<u3::U3>& n2_tagged);
  // Calculate the K^2 matrix elements

  double U3BosonCreationRME(
  const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop, const u3::U3& omegap,
  const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega);

  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map);
  //Calculates the K matrix
  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map, vcs::MatrixCache& Kinv_matrix_map);
  // Generates K matrices and Kinv matrices, for A<6

  
  using KmatrixMap = std::unordered_map<u3::U3, std::array<basis::OperatorBlock<double>, 2>>;
  
  KmatrixMap
  GenerateKMatrices(
    const u3::U3& sigma,
    const std::map<u3::U3, MultiplicityTagged<u3::U3>::vector>& u3_subspaces
  );



}  //  namespace

#endif
