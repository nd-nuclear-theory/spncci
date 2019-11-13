/****************************************************************
  sp3r.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>

#include "gsl/gsl_sf.h"


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

    // append remaining entries
      for (int N=0; N<=Nn_max; N+=2)
        for (int a=N-2; a>=0; a-=2)
          for (int b=2*(a/4); b>=std::max((2*a-N),0); b-=2)
           poly_labels.push_back(u3::U3(N-a,a-b,b));

         return poly_labels;

       }

  ////////////////////////////////////////////////////////////////
  // space and subspace indexing
  ////////////////////////////////////////////////////////////////

       U3Subspace::U3Subspace(const u3::U3& omega, int upsilon_max)
       {
    // set subspace labels
        labels_ = omega;
        upsilon_max_=upsilon_max;
      }

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

      // void U3Subspace::Init(const std::vector<MultiplicityTagged<u3::U3>>& state_set)
      // {
      //   // for each (n,rho_max) entry belonging to this omega subspace
      //   for (auto state : state_set)
      //     PushStateLabels(state);
      // }

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



    Sp3RSpace::Sp3RSpace(const u3::U3& sigma, int Nn_max, bool restrict_sp3r_to_u3_branching)
    {

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
      for (auto n_iter = n_vec.begin(); n_iter != n_vec.end(); ++n_iter)
      {
       u3::U3 n = (*n_iter);
       MultiplicityTagged<u3::U3>::vector omega_tagged_vec = KroneckerProduct(sigma,n);
       for (
        auto omega_tagged_iter = omega_tagged_vec.begin();
        omega_tagged_iter != omega_tagged_vec.end();
        ++omega_tagged_iter
        )
       {
        // convert (omega,rho_max) to omega -> (n,rho_max)
        MultiplicityTagged<u3::U3> omega_tagged = (*omega_tagged_iter);
        u3::U3 omega = omega_tagged.irrep;
        int rho_max = omega_tagged.tag;
        for(int rho=1; rho<=rho_max; ++rho)
          states[omega].push_back(MultiplicityTagged<u3::U3>(n,rho));
       }
     }

    // scan through spanakopita for subspaces
     for(auto it=states.begin(); it!=states.end(); ++it)
     {
      // retrieve omega key of this group of states
      u3::U3 omega = it->first;

      // retrieve upsilon multiplicity states
      const MultiplicityTagged<u3::U3>::vector& multiplicities=it->second;

      // emplace subspace into space
      EmplaceSubspace(omega,multiplicities.size(),multiplicities);
     }

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
