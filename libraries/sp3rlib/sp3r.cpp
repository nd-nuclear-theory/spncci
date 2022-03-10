/****************************************************************
  sp3r.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>


#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3boson.h"
#include "sp3rlib/vcs.h"
  namespace sp3r
  {

  // ////////////////////////////////////////////////////////////////
  // // Sp(3,R) raising polynomial
  // ////////////////////////////////////////////////////////////////
    // TauConjugate computes the first two rows of the partition conjugate 
    // of tau = [f1-f3,f2-f3] = [lambda+mu,mu], defined in 
    // jpa-18-1985-939-Rowe. Note, in ref., tau is referred to as lambda. 
    std::array<int,2> TauConjugate(const u3::U3& sigma)
    {
      std::array<int,2> tau_tilde ={0};
      const auto& lambda = sigma.SU3().lambda();
      const auto& mu = sigma.SU3().mu();

      if((lambda+mu)>=1)
        tau_tilde[0]++;

      if((lambda+mu)>=2)
        tau_tilde[1]++;

      if(mu>=1)
        tau_tilde[0]++;

      if(mu>=2)
        tau_tilde[1]++;

      return tau_tilde;
    }

  bool IsUnitary(const u3::U3& sigma)
  // Based on criterion for unitary irrep in jpa-18-1985-939-Rowe
  // equations (2.7) and (2.8).  Note, first criterion of 
  // conjugate{tau}_1<=2 always met for U(3) irrep, so we only check the 
  // second criterion conjugate{tau}_1+conjugate{tau}_2 <=2*f3.
  {
    HalfInt f3=sigma.f3();
    auto tau_conjugate = TauConjugate(sigma);
    return (tau_conjugate[0]+tau_conjugate[1])<=2*f3;
  }


  bool ModifySp3RBranching(const u3::U3& sigma)
  // From section 5 of jpa-18-1985-939-Rowe,
  // there are U(3) irreps in the U(3) boson basis
  // obtained via laddering which have no counter-part
  // in the Sp(3,R) basis and thus must be removed.
  // Such cases occur when f3<3.
  //
  // Basis must be restricted if f3<3 and 
  // conjugate{tau}_2 > (2*f3-3)
  {
    const HalfInt& f3 = sigma.f3();
    if(f3<3)
      return true;
    else
      return false;
  }


  ////////////////////////////////////////////////////////////////
  // space and subspace indexing
  ////////////////////////////////////////////////////////////////

       U3Subspace::U3Subspace(const u3::U3& omega, int upsilon_max)
        : BaseSubspace{{omega}}, upsilon_max_{upsilon_max}
       {}

       void U3Subspace::Init(const MultiplicityTagged<u3::U3>::vector& multiplicities)
      {
        // for each (n,rho_max) entry belonging to this omega subspace
        for (auto n_rho : multiplicities)
          PushStateLabels(n_rho);
      }


     std::string U3Subspace::DebugStr() const
     {
      std::ostringstream ss;

    // print subspace labels
      u3::U3 omega = labels();
      ss << fmt::format("subspace {}\nupsilon_max {}  dimension {}",labels().Str(),upsilon_max(),dimension())<<std::endl;

    // enumerate state labels within subspace
      for (int i_state=0; i_state<size(); ++i_state)
      {
        MultiplicityTagged<u3::U3> n_rho = GetStateLabels(i_state);
        ss << "  " << i_state << " " << n_rho.Str() << std::endl;
      }

      return ss.str();
    }


    double zero_threshold = 1e-6;
    Sp3RSpace::Sp3RSpace(const u3::U3& sigma, int Nn_max)
    {
      // Make sure sigma is a unitary irrep
      assert(sp3r::IsUnitary(sigma));
      bool modify_sp3r_to_u3_branching = sp3r::ModifySp3RBranching(sigma);
      assert(modify_sp3r_to_u3_branching==false);
      // set space labels
      sigma_ = sigma;
      Nn_max_ = Nn_max;


      // set up container for states
      std::map<u3::U3,MultiplicityTagged<u3::U3>::vector> states;

      // find all raising polynomials
      std::vector<u3::U3> n_vec = u3boson::RaisingPolynomialLabels(Nn_max);

      // enumerate states
      //
      // for each raising polynomial n
      //   obtain all allowed couplings omega (sigma x n -> omega)
      //     (with their multiplicities rho_max)
      //   for each allowed coupling omega
      //      store to states multimap as key value pair
      //        omega -> (n,rho_max)
      for (const auto& n : n_vec)
        {
          MultiplicityTagged<u3::U3>::vector omega_tagged_vec = KroneckerProduct(sigma,n);
          for(const auto& [omega,rho_max] : omega_tagged_vec)
            for(int rho=1; rho<=rho_max; ++rho)
              states[omega].push_back({n,rho});
        }

      // scan through spanakopita for subspaces and add subspace
     for(const auto& [omega,multiplicities] : states)
        EmplaceSubspace(omega,multiplicities.size(),multiplicities);
   }


  ////////////////////////////////////////////////////////////////////////
  // Construction from U3BosonSpace
  ////////////////////////////////////////////////////////////////////////
  void U3Subspace::Init(const u3boson::U3Subspace& u3boson_subspace)
    {
      for(int i=0; i<u3boson_subspace.size(); ++i)
        {
          const auto& u3boson_state = u3boson_subspace.GetState(i);
          const u3::U3& n = u3boson_state.n();
          for(int rho=1; rho<=u3boson_state.rho_max(); ++rho)
            {
              PushStateLabels({n,rho});
            }
        }
    }


  Sp3RSpace::Sp3RSpace(
    const u3::U3& sigma, const int Nn_max,
    const u3boson::U3BosonSpace& u3boson_space,
    const bool subspace_labels_only
    )
    {
      Nn_max_ = Nn_max;
      sigma_=sigma;

      // Check that sigma is an LGI of a unitary Sp(3,R) irrep
      assert(IsUnitary(sigma));

      // If constructing the full space, then compute K matrices and
      // use K matrices to get upsilon_max.  If subspace labels only,
      // K matrices only need to be computed if branching must be restricted.
      vcs::KmatrixMap K_matrices;
      bool get_upsilon_from_K=false;
      if(subspace_labels_only==false || sp3r::ModifySp3RBranching(sigma))
        {
          double zero_threshold = 1e-12;
          K_matrices = vcs::GenerateKmatrices(sigma,u3boson_space,zero_threshold);
          get_upsilon_from_K=true;
        }

      for(const auto& u3boson_subspace : u3boson_space)
        {
          const u3::U3& omega = u3boson_subspace.omega();
          int upsilon_max = get_upsilon_from_K?K_matrices[omega][0].cols():u3boson_subspace.dimension();

          if(upsilon_max==0)
            continue;

          // Labels only
          if(subspace_labels_only)
            PushSubspace(U3Subspace(omega,upsilon_max));
          // Full space
          else
            PushSubspace(U3Subspace(
                omega,upsilon_max,
                u3boson_subspace,
                std::move(K_matrices[omega][0]),
                std::move(K_matrices[omega][1])
              ));
        }
    }

////////////////////////////////////////////////////////////////////////



   std::string Sp3RSpace::DebugStr() const
   {
    std::ostringstream ss;

    // print space labels
    ss << "space sigma " << sigma_.Str() << " Nn_max " << Nn_max_ << std::endl;

    // iterate over subspaces
    for (int i_subspace=0; i_subspace<size(); ++i_subspace)
    {
        // set up reference to subspace of interest
        //
        // Using a reference avoids copying the U3Subspace object (and
        // all its lookup tables).
      const sp3r::U3Subspace& subspace = GetSubspace(i_subspace);

        // generate debug information for subspace
      ss << subspace.DebugStr();
    }

    return ss.str();
  }


  // std::vector<int> PartitionIrrepByNn(const sp3r::Sp3RSpace& irrep, const int Nmax)
  // {
  //   // partition irreps by Nn
  //   HalfInt Ns=irrep.GetSubspace(0).labels().N();
  //   int Nn_last=-1;
  //   std::vector<int> IrrepPartionN;
  //   for(int i=0; i<irrep.size(); i++ )
  //     {
  //       u3::U3 omega=irrep.GetSubspace(i).labels();

  //       if ( Nn_last!=int(omega.N()-Ns) )
  //         {
  //           IrrepPartionN.push_back(i);
  //           Nn_last=int(omega.N()-Ns);
  //         }
  //     }
  //   return IrrepPartionN;
  // }


}  // namespace
