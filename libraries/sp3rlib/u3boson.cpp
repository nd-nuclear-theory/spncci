/****************************************************************
  u3boson.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "sp3rlib/u3boson.h"

#include <numeric>
#include "fmt/format.h"
#include <Eigen/Eigenvalues>
#include "sp3rlib/u3coef.h"
#include "mcutils/eigen.h"
#include "cppitertools/itertools.hpp"

namespace u3boson{
  std::vector<u3::U3> RaisingPolynomialLabels(int Nn_max)
  {
    std::vector<u3::U3> poly_labels;

  // prime with Nn=0 entry
    if (Nn_max>=0)
      poly_labels.push_back(u3::U3(0,0,0));

  // append remaining entries
    for (int N=2; N<=Nn_max; N+=2)
      for (int a=N-2; a>=0; a-=2)
        for (int b=2*(a/4); b>=std::max((2*a-N),0); b-=2)
         poly_labels.push_back(u3::U3(N-a,a-b,b));

       return poly_labels;

     }

  U3BosonSpace::U3BosonSpace(const u3::U3& sigma, const int Nn_max)
      : sigma_{sigma}, Nn_max_(Nn_max)
  {
    // find all raising polynomials
    std::vector<u3::U3> n_vector = u3boson::RaisingPolynomialLabels(Nn_max);

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
    }

    // Create U3Subspace
    for (const auto& [omega, n_rho_vector] : states)
      PushSubspace(U3Subspace(omega,n_rho_vector));
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

    U3BosonSectors::U3BosonSectors(
      const U3BosonSpace& space,
      const u3::U3& omega0
    )
      : BaseSectors{space,space}, omega0_(omega0)
    {
      // omega0_ = omega0;
      for(int bra_index=0; bra_index<space.size(); ++bra_index)
        for(int ket_index=0; ket_index<space.size(); ++ket_index)
          {
            const u3::U3& omega_bra = space.GetSubspace(bra_index).omega();
            const u3::U3& omega_ket = space.GetSubspace(ket_index).omega();
            int multiplicity = u3::OuterMultiplicity(omega_ket,omega0,omega_bra);

            if(multiplicity>0)
              PushSector(bra_index,ket_index,multiplicity);
          }
    }

  // Eq. 8 of G. Rosensteel and D. J. Rowe, J. Math. Phys. 24 (1983) 2461.
  // for N0 = 3/2.
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

  // Eq. 9 of G. Rosensteel and D. J. Rowe, J. Math. Phys. 24 (1983) 2461.
  double U3BosonCreationRME(
    const u3::U3& sigmap, const u3::U3& np, unsigned int rhop, const u3::U3& omegap,
    const u3::U3& sigma,  const u3::U3& n,  unsigned int rho,  const u3::U3& omega
  )
  {

    bool allowed = sigma==sigmap;
    allowed &= u3::OuterMultiplicity(omega.SU3(),{2u,0u},omegap.SU3());
    allowed &= u3::OuterMultiplicity(n.SU3(),{2u,0u},np.SU3());

    double rme=0.0;
    if (allowed)
      {
        rme = ParitySign(u3::ConjugationGrade(omegap)+u3::ConjugationGrade(omega))
          *u3::U({2u,0u},n.SU3(),omegap.SU3(),sigma.SU3(),np.SU3(),1,rhop,omega.SU3(),rho,1)
          *u3boson::BosonCreationRME(np,n);
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


  double U3BosonAnnihilationRME(
    const u3::U3& sigmap, const u3::U3& np, unsigned int rhop, const u3::U3& omegap,
    const u3::U3& sigma,  const u3::U3& n,  unsigned int rho,  const u3::U3& omega
  )
  {
    // Eq. 7 of G. Rosensteel and D. J. Rowe, J. Math. Phys. 24 (1983) 2461.
    double rme
      =ParitySign(u3::ConjugationGrade(omegap)+u3::ConjugationGrade(omega))
        *std::sqrt(1.*u3::dim(omega.SU3())/u3::dim(omegap.SU3()))
        *U3BosonCreationRME(sigma,n,rho,omega,sigmap,np,rhop,omegap);
    return rme;
  }


  basis::OperatorBlock<double>
  U3BosonCreationOperator(
    const u3::U3& sigma,
    const u3boson::U3Subspace& bra_subspace,
    const u3boson::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    )
  {
    const auto& omega_bra = bra_subspace.omega();
    const auto& omega_ket = ket_subspace.omega();

    basis::OperatorBlock<double> operator_block;
    operator_block = basis::OperatorBlock<double>::Zero(bra_subspace.dimension(),ket_subspace.dimension());

    for(int bra_state_index=0; bra_state_index<bra_subspace.size(); ++bra_state_index)
      for(int ket_state_index=0; ket_state_index<ket_subspace.size(); ++ket_state_index)
        {
          const auto& bra_state = bra_subspace.GetState(bra_state_index);
          const auto& ket_state = ket_subspace.GetState(ket_state_index);

          const auto& n_bra = bra_state.n();
          const auto& n_ket = ket_state.n();
          const auto& rho_max_bra = bra_state.rho_max();
          const auto& rho_max_ket = ket_state.rho_max();
          for(int rho_bra=1; rho_bra<=rho_max_bra; ++rho_bra)
            for(int rho_ket=1; rho_ket<=rho_max_ket; ++rho_ket)
              {
                int row = bra_subspace.GetStateOffset(bra_state_index,rho_bra);
                int col = ket_subspace.GetStateOffset(ket_state_index,rho_ket);

                operator_block(row,col)
                  =U3BosonCreationRME(
                      sigma, n_bra,rho_bra,omega_bra,
                      sigma, n_ket,rho_ket,omega_ket
                    );
              }
        }

    return operator_block;
  }

  basis::OperatorBlock<double>
  U3BosonAnnihilationOperator(
    const u3::U3& sigma,
    const u3boson::U3Subspace& bra_subspace,
    const u3boson::U3Subspace& ket_subspace,
    u3::UCoefCache& u_coef_cache
    )
  {
    const auto& x_bra=bra_subspace.omega().SU3();
    const auto& x_ket=ket_subspace.omega().SU3();

    double conjugation_factor
      =ParitySign(u3::ConjugationGrade(x_bra)+u3::ConjugationGrade(x_ket))
        *std::sqrt(1.*u3::dim(x_ket)/u3::dim(x_bra));

    return conjugation_factor*U3BosonCreationOperator(sigma,ket_subspace,bra_subspace,u_coef_cache).transpose();

  }


}  //  namespace
