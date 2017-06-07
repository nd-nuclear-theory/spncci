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

  // CCoef(u3::U3& n, int q, u3:: U3& np)
  // //Don't need at present.  Will come back and finish later (aem). 

  // // Coefficient of fractional parentage for expanding the raising polynomial into a
  // // polynomial of a single Jacobi coordinate and a polynomial of all other coordinates 
  // {
  //   double coef=0;
  //   int N=int(n.N());
  //   if (n.SU3()==u3::SU3(N,0))&&(np.SU3()==u3::SU3(N-q,0))
  //     coef=sqrt(
  //       gsl_sf_fact(N/2)*gsl_sf_fact(q)
  //       /(pow(2.,q/2)*gsl_sf_fact(N/2-q/2)*gsl_sf_fact(q/2))
  //       );
  //   else
  //   {
  //     if ((n.f1-2)>=n.f2)
  //       u3::U3 nb(n.f1-2,n.f2,n.f3);
  //     else if ((n.f2-2)>=f.n3)
  //       u3::U3 nb(n.f1,n.f2-2,n.f3);
  //     else if ((n.f3-2)>=0)
  //       u3::U3 nb(n.f1,n.f2,n.f3-2);
  //     // if none of the above returns true, coef is zero. 
  //     else
  //       return coef;
  //     // otherwise, calculate coef
  //  }
  //}
  
  void GenerateBCoefCache(BCoefCache& cache, int Nmax)
    {
      // TODO: Finish debugging
      std::vector<u3::U3> poly_vector=RaisingPolynomialLabels(Nmax);
      for(auto n1 : poly_vector)
        {
          MultiplicityTagged<u3::SU3>::vector n1p_set=u3::KroneckerProduct(n1.SU3(),u3::SU3(0,2));
          int N1p=int(n1.N()-2);
          u3::U3 n1p;
          bool continue_flag=true;

          for(int i=0; i<n1p_set.size(); ++i)
            {
              n1p=u3::U3(N1p,n1p_set[i].irrep);
              if((n1p.Valid())&&(n1p.SU3().lambda()%2==0)&&(n1p.SU3().mu()%2==0))
                {    
                  continue_flag=false;
                  break;
                }
            }
          if(continue_flag)
            continue;
          double coef1=vcs::BosonCreationRME(n1,n1p);
          for(auto n2 : poly_vector)
            {
              MultiplicityTagged<u3::SU3>::vector n3_set=u3::KroneckerProduct(n1.SU3(),n2.SU3());
              MultiplicityTagged<u3::SU3>::vector n3p_set=u3::KroneckerProduct(n1p.SU3(),n2.SU3());
              int N3=int(n1.N()+n2.N());
              int N3p=N3-2;
              for(auto n3_tagged :n3_set)
                {
                  if(not u3::U3::ValidLabels(N3, n3_tagged.irrep))
                    continue;
                  u3::U3 n3(N3,n3_tagged.irrep);
                  if((n3.SU3().lambda()%2!=0)||(n3.SU3().mu()%2!=0))
                    continue;
                  // If n1 or n2 are identity, the coef is 1
                  if((n1==u3::U3(0,0,0))||(n2==u3::U3(0,0,0)))
                    {
                      cache[BCoefLabels(n1,n2,n3,1)]=1;
                      continue;
                    }
                  int rho_max(n3_tagged.tag);
                  for(int rho=1; rho<=rho_max; ++rho)
                    {
                  // std::cout<<n1.Str()<<"  "<<n2.Str()<<"  "<<n3.Str()<<"  "<<rho<<std::endl;
                      BCoefLabels labels(n1,n2,n3,rho);
                      double coef2=0;
                      for(auto n3p_tagged : n3p_set)
                        {
                          u3::U3 n3p(N3p,n3p_tagged.irrep);
                          if(not n3p.Valid())
                            continue;
                          if((n3p.SU3().lambda()%2!=0)||(n3p.SU3().mu()%2!=0))
                            continue;
                          if(u3::OuterMultiplicity(n3p.SU3(),u3::SU3(2,0),n3.SU3())==0)
                            continue;
                          double coef3=vcs::BosonCreationRME(n3,n3p);
                          int rhop_max=n3p_tagged.tag;
                          for(int rhop=1; rhop<=rhop_max; ++rhop)
                            {
                              BCoefLabels labelsp(n1p,n2,n3p,rhop);
                              double coef=coef3/coef1*cache[labelsp]
                                *u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1);
                              // std::cout<<"  "<<n1p.Str()<<"  "<<n2.Str()
                              // <<"  "<<n3p.Str()<<"  "<<rhop<<"  "<<coef1<<"  "<<coef3<<"  "<<cache[labelsp]
                              // <<"  "<<u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1)
                              // <<"  "<<coef
                              // <<std::endl;
                              assert(cache.count(labelsp));
                              cache[labels]
                              +=coef3/coef1*cache[labelsp]
                                *u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1);
                              // std::cout<<"    "<<cache[labels]<<std::endl;
                              // coef2+=cache[labelsp]
                              // *u3::U(u3::SU3(2,0),n1p.SU3(),n3.SU3(),n2.SU3(),n1.SU3(),1,rho,n3p.SU3(),rhop,1);
                            }
                        }
                    }
                }
            }
        }
    }
   


  ////////////////////////////////////////////////////////////////
  // space and subspace indexing
  ////////////////////////////////////////////////////////////////

       U3Subspace::U3Subspace(const u3::U3& omega)
       {
    // set subspace labels
        labels_ = omega;
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

    Sp3RSpace::Sp3RSpace(const u3::U3& sigma, int Nn_max)
    {

    // set space labels
      sigma_ = sigma;
      Nn_max_ = Nn_max;

    // set up container for states
      SpanakopitaType states;

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
         states.insert(SpanakopitaType::value_type(omega,MultiplicityTagged<u3::U3>(n,rho_max)));
         
       }
     }

    // scan through spanakopita for subspaces

     SpanakopitaType::iterator it = states.begin();
     while (it != states.end())
     {
	// retrieve omega key of this group of states
       u3::U3 omega = it->first;

	// find range of this group of states
       SpanakopitaRangeType state_range = states.equal_range(omega);

	// set up subspace
	//
	// Note: Creation followed by push_back means U3Subspace is
	// *copied* into vector.  
	//
	// Option 1: Define a lightweight
	// constructor for U3Subspace, which just saves omega.  Push
	// the (empty) subspace into the Sp3RSpace.  Then use an Init
	// method to actually populate the subspace.  This runs up
	// against "private" protections in BaseSubspace (which could
	// be removed.)
	//
	// Option 2: Use C++11 emplace_back!  However, this requires
	// also modifying BaseSubspace to provide an "EmplaceSubspace"
	// with appropriate arguments, which would be difficult.
       PushSubspace(U3Subspace(omega));
       subspaces_.back().Init(state_range);
	// GetSubspace(size()).Init(state_range); //ILLEGAL since GetSubspace yields const reference

	// move to start of next group of states
       it = state_range.second;
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
