/****************************************************************
 sp3r_operators.cpp

  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/
#include "sp3rlib/sp3r_operator.h"



namespace sp3r
{

  Eigen::MatrixXd  Sp3rRaisingOperator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, 
      const u3::U3& omega
    )
  {
    // Get irrep label
    const u3::U3& sigma=sp3r_space.sigma();

    // Look up  relevant subspaces in Sp(3,R) irrep
    int subspace_index_bra=sp3r_space.LookUpSubspaceIndex(omegap);
    int subspace_index_ket=sp3r_space.LookUpSubspaceIndex(omega);
    const auto& subspace_bra = sp3r_space.GetSubspace(subspace_index_bra);
    const auto& subspace_ket = sp3r_space.GetSubspace(subspace_index_ket);

    //Construct boson matrix 
    int vp_max=subspace_bra.size();
    int v_max=subspace_ket.size();
    basis::OperatorBlock<double> A_boson
      = basis::OperatorBlock<double>::Zero(subspace_bra.dimension(),subspace_ket.dimension());


    for(int ip=0; ip<subspace_bra.size(); ++ip)
      for(int i=0; i<subspace_bra.size(); ++i)
        {
          const auto& [np,rhop] = subspace_bra.GetStateLabels(ip);
          const auto& [n,rho] = subspace_ket.GetStateLabels(i);
          A_boson(ip,i)=u3boson::U3BosonCreationRME(sigma,np,rhop,omegap,sigma,n,rho,omega);
        }

    const basis::OperatorBlock<double>& Kp = subspace_bra.K_matrix();
    const basis::OperatorBlock<double>& Kinv = subspace_bra.Kinv_matrix();
    
    //Calculate matrix element of symplectic raising operator 
    return Kp*A_boson*Kinv;
  }

  Eigen::MatrixXd Sp3rLoweringOperator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, 
      const u3::U3& omega
    )
  {
    int parity_sign=ParitySign(u3::ConjugationGrade(omega.SU3())-u3::ConjugationGrade(omegap.SU3()));
    return parity_sign*sqrt(1.0*u3::dim(omega)/u3::dim(omegap))
            *sp3r::Sp3rRaisingOperator(sp3r_space, omega, omegap);
  }

  Eigen::MatrixXd  U3Operator(
      const sp3r::Sp3RSpace& sp3r_space, 
      const u3::U3& omegap, 
      const u3::U3& omega
    )
  {
    // Look up  relevant subspaces in Sp(3,R) irrep
    int subspace_index_bra=sp3r_space.LookUpSubspaceIndex(omegap);
    int subspace_index_ket=sp3r_space.LookUpSubspaceIndex(omega);

    const sp3r::U3Subspace& subspace_bra =  sp3r_space.GetSubspace(subspace_index_bra);
    const sp3r::U3Subspace& subspace_ket =  sp3r_space.GetSubspace(subspace_index_ket);

    //Construct boson matrix 
    int vp_max=subspace_bra.size();
    int v_max=subspace_ket.size();

    if(omegap==omega)
      {
        double Crme=sqrt(2.0*Casimir2(omega.SU3()));
        return Crme*Eigen::MatrixXd::Identity(v_max,v_max);
      }

    else
      return Eigen::MatrixXd::Zero(vp_max,v_max);
  }

  // double NumberOperatorU3(
  //         const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
  //         const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
  //         )
  // {
  //   double rme=0;
  //   if (
  //     (sigmap==sigma)
  //     && (np_rhop==n_rho)
  //     && (omegap==omega)
  //     )
  //     rme=double(omega.N());

  //   return rme;
  // }
  
  // // TODO:: ANNA Needs to be edited to use Kmatrix 
  // double Sp3rRaisingOperator(
  //         const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
  //         const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
  //         )
  // {
  //   double rme=0;
  //   if(
  //       (sigmap==sigma)
  //       &&(u3::OuterMultiplicity(omega.SU3(),u3::SU3(2,0),omegap.SU3())!=0)
  //       && ((omega.N()+2)==omegap.N())
  //     )
  //       rme=sqrt(vcs::Omega(np_rhop.irrep, omegap)-vcs::Omega(n_rho.irrep, omega))
  //           *vcs::U3BosonCreationRME(sigma,np_rhop,omegap,sigma,n_rho,omega);
  //   return rme;
  // }


  // // TODO:: ANNA Needs to be edited to use Kmatrix 
  // double Sp3rLoweringOperator(
  //         const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
  //         const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
  //         )
  // {
  //   double rme=0;
  //   if(
  //       (sigmap==sigma)
  //       &&(u3::OuterMultiplicity(omega.SU3(),u3::SU3(0,2),omegap.SU3())!=0)
  //       && ((omega.N()-2)==omegap.N())
  //     )
  //       rme=parity(u3::ConjugationGrade(omegap))//.lambda()+omegap.mu+omega.lambda+omega.mu)
  //           *sqrt(1.*u3::dim(omega)/u3::dim(omegap))
  //           *sqrt(vcs::Omega(n_rho.irrep, omega)-vcs::Omega(np_rhop.irrep, omegap))
  //           *vcs::U3BosonCreationRME(sigma,n_rho,omega,sigma,np_rhop,omegap);
  //   return rme;
  // }
}
