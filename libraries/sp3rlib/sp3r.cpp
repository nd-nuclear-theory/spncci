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
#include "sp3rlib/vcs.h"
  namespace sp3r
  {

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) raising polynomial
  ////////////////////////////////////////////////////////////////

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

      void U3Subspace::Init(const SpanakopitaRangeType& state_range)
      {
    // for each (n,rho_max) entry belonging to this omega subspace
        for (auto it = state_range.first; it != state_range.second; ++it)
        {
	// extract state labels
	//   from omega -> (n,rho_max) entry
         MultiplicityTagged<u3::U3> n_rho_max = it->second;
         u3::U3 n = n_rho_max.irrep;
         int rho_max = n_rho_max.tag;

	// push state (n,rho) labels into subspace indexing
	//   enumerating all multiplicity indices
         for (int rho=1; rho<=rho_max; ++rho)
         {
           PushStateLabels(MultiplicityTagged<u3::U3>(n,rho));
         }
       }
     }


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
      ss << "subspace " << omega.Str() << std::endl;

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
      // set space labels
      sigma_ = sigma;
      Nn_max_ = Nn_max;


      // set up container for states
      std::map<u3::U3,MultiplicityTagged<u3::U3>::vector> states;

      // find all raising polynomials
      std::vector<u3::U3> n_vec = RaisingPolynomialLabels(Nn_max);

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

      bool modify_sp3r_to_u3_branching = sp3r::ModifySp3RBranching(sigma);
      if(modify_sp3r_to_u3_branching)
        {
          std::map<u3::U3,MultiplicityTagged<u3::U3>::vector> keep_states;
          for(const auto&[omega1,n_rho_vec1] : states)
            {
              if (omega1==sigma)
                {
                  keep_states[omega1].push_back({{0,{0,0}},1});
                  continue;
                }

              for(const auto& [n1,rho1] : n_rho_vec1)
                {
                  bool keep_irrep=false;
                  auto omega2_vec = u3::KroneckerProduct(omega1, u3::U3(-2,{0,2}));
                  for(const auto& [omega2,dummy] : omega2_vec)
                    {
                      assert(dummy==1);
                      for(const auto& [n2,rho2] : states[omega2])
                        {
                          double DeltaOmega = vcs::Omega(n1,omega1)-vcs::Omega(n2,omega2);
                          // if(DeltaOmega<=0)
                          // {
                          //   double u3boson_rme = vcs::U3BosonCreationRME(sigma,n1,rho1,omega1,sigma,n2,rho2,omega2);
                          //   fmt::print("{} {} {}  {} {} {}\n",omega1,n1,rho1,omega2,n2,rho2);
                          //   fmt::print("{}  {}\n",DeltaOmega,u3boson_rme);
                          // }
                          // assert(DeltaOmega>=0);
                          if(DeltaOmega>0)
                            {
                              double u3boson_rme = vcs::U3BosonCreationRME(sigma,n1,rho1,omega1,sigma,n2,rho2,omega2);
                              if (u3boson_rme>zero_threshold)
                                {
                                  keep_irrep=true;
                                  break;
                                }
                            }
                        }
                    }
                  if(keep_irrep)
                    keep_states[omega1].emplace_back(n1,rho1);
                }
            }
          // Overwrite list of states
          states=keep_states;
        }

      // scan through spanakopita for subspaces and add subspace
     for(const auto& [omega,multiplicities] : states)
        EmplaceSubspace(omega,multiplicities.size(),multiplicities);
   }


  Sp3RSpace::Sp3RSpace(
    const u3::U3& sigma, int Nn_max,
    const RestrictedSpanakopitaType& spanakopita
  )
  {
    // set space labels
    sigma_ = sigma;
    Nn_max_ = Nn_max;

    // scan through spanakopita for subspaces
    for(auto it=spanakopita.begin(); it!=spanakopita.end(); ++it)
      {
        // retrieve omega key of this group of states
        u3::U3 omega=it->first.irrep;
        int upsilon_max=it->first.tag;
        const auto& states = it->second;
        EmplaceSubspace(omega,upsilon_max,states);
      }
  }


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


  std::vector<int> PartitionIrrepByNn(const sp3r::Sp3RSpace& irrep, const int Nmax)
  {
    // partition irreps by Nn
    HalfInt Ns=irrep.GetSubspace(0).labels().N();
    int Nn_last=-1;
    std::vector<int> IrrepPartionN;
    for(int i=0; i<irrep.size(); i++ )
      {
        u3::U3 omega=irrep.GetSubspace(i).labels();

        if ( Nn_last!=int(omega.N()-Ns) )
          {
            IrrepPartionN.push_back(i);
            Nn_last=int(omega.N()-Ns);
          }
      }
    return IrrepPartionN;
  }


}  // namespace
