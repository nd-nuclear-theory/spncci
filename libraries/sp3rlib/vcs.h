  /****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy on 3/9/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py.

****************************************************************/

#ifndef VCS_H_
#define VCS_H_


#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues>  


#include "sp3rlib/u3coef.h" 
#include "sp3rlib/sp3r.h"  
  


namespace vcs
{
  inline double Omega(const u3::U3& n, const u3::U3& omega)

  // Calculate Omega factor used in Kmatrix calculations.
  //
  // Based on protopye vcs.py and equation given in
  //   D. J. Rowe, J. Math Phys. 25 (1984) 2662. 
  //
  // Returns:
  //   (double) : Omega factor
  {
    double value=0;
    value += (2*sqr(double(omega.f1))-sqr(double(n.f1))+8*(double(omega.f1)-double(n.f1))-2*(2*double(omega.f1)-double(n.f1)));
    value += (2*sqr(double(omega.f2))-sqr(double(n.f2))+8*(double(omega.f2)-double(n.f2))-4*(2*double(omega.f2)-double(n.f2)));
    value += (2*sqr(double(omega.f3))-sqr(double(n.f3))+8*(double(omega.f3)-double(n.f3))-6*(2*double(omega.f3)-double(n.f3)));
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
  const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
  const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega);

  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep,std::map<u3::U3,Eigen::MatrixXd>& K_matrix_map);
  //Calculates the K matrix 	

}  //  namespace

#endif
