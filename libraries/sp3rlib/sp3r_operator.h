  /****************************************************************
  sp3r_operators.h                       

  Sp(3,R) generator and operator reduced matrix elements.

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy and M. A. Caprio on 6/6/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py and sp3r.py.

****************************************************************/
#ifndef SP3R_OPERATORS_H_
#define SP3R_OPERATORS_H_

#include "sp3rlib/u3.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/vcs.h"

namespace sp3r
{
  // double NumberOperatorU3(
  //         const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
  //         const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
  //         );

  // double Sp3rRaisingOperator(
  //         const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
  //         const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
  //         );

  // double Sp3rLoweringOperator(
  //         const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
  //         const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
  //         );

  Eigen::MatrixXd  Sp3rRaisingOperator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, 
      const u3::U3& omega, 
      const vcs::MatrixCache& K_matrices
    );
  // Reduced matrix elements of symplectic raising operator between states in
  // omegap and omega subspace of Sp(3,R) irrep defiend by sp3r_space 

  Eigen::MatrixXd Sp3rLoweringOperator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, 
      const u3::U3& omega, 
      const vcs::MatrixCache& K_matrices
    );
  // Reduced matrix elements of symplectic lowering operator between states in
  // omegap and omega subspace of Sp(3,R) irrep defiend by sp3r_space 


  Eigen::MatrixXd  U3Operator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, 
      const u3::U3& omega
    );
  // Reduced matrix elements of U(3) generators operator between states in
  // omegap and omega subspace of Sp(3,R) irrep defiend by sp3r_space 


}

#endif