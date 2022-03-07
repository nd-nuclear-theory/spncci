/****************************************************************
  vcs.cpp

  Anna E. McCoy
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/
#include "sp3rlib/vcs.h"

#include <numeric>
#include "fmt/format.h"
#include <Eigen/Eigenvalues>
#include "sp3rlib/u3coef.h"
#include "mcutils/eigen.h"
#include "cppitertools/itertools.hpp"

namespace vcs{

  std::vector<u3::U3> RaisingPolynomialLabels(int Nn_max)
  {
    std::vector<u3::U3> poly_labels;

  // prime with Nn=0 entry
    if (Nn_max>=0)
      poly_labels.push_back(u3::U3(0,0,0));

  // Should't N start with 2?
  // append remaining entries
    for (int N=0; N<=Nn_max; N+=2)
      for (int a=N-2; a>=0; a-=2)
        for (int b=2*(a/4); b>=std::max((2*a-N),0); b-=2)
         poly_labels.push_back(u3::U3(N-a,a-b,b));

       return poly_labels;

     }


  U3BosonSpace::U3BosonSpace(const u3::U3& sigma, const int Nn_max)
      : BaseSpace{sigma}
  {
    Nn_max_=Nn_max;

    // find all raising polynomials
    std::vector<u3::U3> n_vector = vcs::RaisingPolynomialLabels(Nn_max);

    // temporary container
    std::map<u3::U3, MultiplicityTagged<u3::U3>::vector> states;

    // For each raising polynomial n obtain all allowed couplings
    // omega (sigma x n -> omega) with outer multiplicities rho_max.
    for (const auto& n : n_vector)
    {
      MultiplicityTagged<u3::U3>::vector omega_rhomax_vector =
          u3::KroneckerProduct(sigma, n);
      for (const auto& [omega, rho_max] : omega_rhomax_vector)
        states[omega].emplace_back(n, rho_max);

      // Create U3Subspace
      for (const auto& [omega, n_rho_vector] : states)
        PushSubspace(U3Subspace(omega,n_rho_vector));
    }
  }


  std::string U3Subspace::DebugStr() const
    {
      std::string ss;

      // print subspace labels
      u3::U3 omega = std::get<0>(labels());
      ss +=fmt::format(" subspace {}\n",omega);

      // enumerate state labels within subspace
      for (int i_state=0; i_state<size(); ++i_state)
      {
        const auto& state = GetState(i_state);
        ss += fmt::format("  {} : [{},{}]\n",i_state,state.n(),state.rho_max());
      }

      return ss;
    }

   std::string U3BosonSpace::DebugStr() const
   {
    std::string ss;

    // print space labels
    ss += fmt::format("space sigma {} Nn_max {}\n",sigma(),Nn_max());
    
    // iterate over subspaces
    for (int i_subspace=0; i_subspace<size(); ++i_subspace)
    {
      const auto& subspace = GetSubspace(i_subspace);

        // generate debug information for subspace
      ss += subspace.DebugStr();
    }

    return ss;
  }



  double BosonCreationRME(const u3::U3& np, const u3::U3& n)
  //  SU(3) Reduced matrix element of a^\dagger boson creation operator
  {
    double rme=0;
    const auto n1=int(n.f1());
    const auto n2=int(n.f2());
    const auto n3=int(n.f3());
    const auto n1p=int(np.f1());
    const auto n2p=int(np.f2());
    const auto n3p=int(np.f3());

    if((n1p==(n1+2))&&(n2p==n2)&&(n3p==n3))
      rme=std::sqrt((n1+4)*(n1-n2+2)*(n1-n3+3)/(2.*(n1-n2+3)*(n1-n3+4)));

    else if ((n1p==n1)&&(n2p==(n2+2))&&(n3p==n3))
      rme=std::sqrt((n2+3)*(n1-n2)*(n2-n3+2)/(2.*(n1-n2-1)*(n2-n3+3)));

    else if ((n1p==n1)&&(n2p==n2)&&(n3p==(n3+2)))
      rme=std::sqrt((n3+2)*(n2-n3)*(n1-n3+1)/(2.*(n1-n3)*(n2-n3-1)));

    return rme;

  }

  double U3BosonCreationRME(
    const u3::U3& sigmap, const u3::U3& np, unsigned int rhop, const u3::U3& omegap,
    const u3::U3& sigma,  const u3::U3& n,  unsigned int rho,  const u3::U3& omega
  )
  {
    unsigned int rho0_max=u3::OuterMultiplicity(omega.SU3(),{2,0},omegap.SU3());
    unsigned int rhon_max=u3::OuterMultiplicity(n.SU3(),{2,0},np.SU3());
    
    bool allowed = sigma==sigmap;
    allowed &= u3::OuterMultiplicity(omega.SU3(),{2,0},omegap.SU3());
    allowed &= u3::OuterMultiplicity(n.SU3(),{2,0},np.SU3());

    double rme=0.0;
    if (allowed)
      {
        rme = ParitySign(u3::ConjugationGrade(omegap)+u3::ConjugationGrade(omega))
          *u3::U(u3::SU3(2,0),n.SU3(),omegap.SU3(),sigma.SU3(),np.SU3(),1,rhop,omega.SU3(),rho,1)
          *BosonCreationRME(np,n); 
      }
    
    return rme;
  }

  double U3BosonCreationRME(
    const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
    const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
    )
  {
    const auto& [n,rho] = n_rho;
    const auto& [np,rhop] = np_rhop;
    return U3BosonCreationRME(sigmap,np,rhop,omegap,sigma,n,rho,omega);
  } 

}  //  namespace
