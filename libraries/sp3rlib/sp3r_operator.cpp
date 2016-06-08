/****************************************************************
 sp3r_operators.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "sp3rlib/sp3r_operator.h"

#include "sp3rlib/vcs.h"

  double NumberOperatorU3(
          const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
          const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
          )
  {
    double rme=0;
    if (
      (sigmap==sigma)
      && (np_rhop==n_rho)
      && (omegap==omega)
      )
      rme=double(omega.N());

    return rme;
  }
  
  // TODO:: ANNA Needs to be edited to use Kmatrix 
  double Sp3rRaisingOperator(
          const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
          const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
          )
  {
    double rme=0;
    if(
        (sigmap==sigma)
        &&(u3::OuterMultiplicity(omega.SU3(),u3::SU3(2,0),omegap.SU3())!=0)
        && ((omega.N()+2)==omegap.N())
      )
        double rme=sqrt(vcs::Omega(np_rhop.irrep, omegap)-vcs::Omega(n_rho.irrep, omega))
                    *vcs::U3BosonCreationRME(sigma,np_rhop,omegap,sigma,n_rho,omega);
    return rme;
  }


  // TODO:: ANNA Needs to be edited to use Kmatrix 
  double Sp3rLoweringOperator(
          const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
          const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
          )
  {
    double rme=0;
    if(
        (sigmap==sigma)
        &&(u3::OuterMultiplicity(omega.SU3(),u3::SU3(0,2),omegap.SU3())!=0)
        && ((omega.N()-2)==omegap.N())
      )
        double rme=parity(u3::ConjugationGrade(omegap))//.lambda()+omegap.mu+omega.lambda+omega.mu)
                    *sqrt(1.*u3::dim(omega)/u3::dim(omegap))
                    *sqrt(vcs::Omega(n_rho.irrep, omega)-vcs::Omega(np_rhop.irrep, omegap))
                    *vcs::U3BosonCreationRME(sigma,n_rho,omega,sigma,np_rhop,omegap);
    return rme;
  }