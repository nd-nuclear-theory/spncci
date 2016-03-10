/****************************************************************
  sp3r.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "sp3rlib/sp3r.h"

namespace sp3r
{

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) raising polynomial
  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  // space and subspace indexing
  ////////////////////////////////////////////////////////////////

  U3Subspace::U3Subspace(const u3::U3& w)
  {
    // set subspace labels
    labels_ = w;
  }

  void U3Subspace::Init(const SpanakopitaRangeType& state_range)
  {
    // set up indexing of states
    for (auto it = state_range.first; it != state_range.second; ++it)
      {
	// extract state labels
	MultiplicityTagged<u3::U3> n_rho_max = it->second;
	u3::U3 n = n_rho_max.irrep;
	int rho_max = n_rho_max.tag;

	// push states
	//   enumerating all multiplicity indices
	for (int rho=1; rho<=rho_max; ++rho)
	  {
	    PushStateLabels(MultiplicityTagged<u3::U3>(n,rho));
	  }
      }
  }

  Sp3RSpace::Sp3RSpace(const u3::U3& sigma, int Nn_max)
  {

    // set space labels
    sigma_ = sigma;
    Nn_max_ = Nn_max;

    // set up container for states
    SpanakopitaType states;

    // find all raising polynomials
    std::vector<u3::U3> n_vec;
    // std::vector<u3::U3> n_vec = RaisingPolynomialLabels(Nn_max);

    // enumerate states
    //
    // for each raising polynomial n
    //   obtain all allowed couplings w (sigma x n -> w) 
    //     (with their multiplicities rho_max)
    //   for each allowed coupling w
    //      store to states multimap as key value pair
    //        w -> (n,rho_max)
    for (auto n_iter = n_vec.begin(); n_iter != n_vec.end(); ++n_iter)
      {
	u3::U3 n = (*n_iter);
	MultiplicityTagged<u3::U3>::vector w_tagged_vec = KroneckerProduct(sigma,n);
	for (
	     auto w_tagged_iter = w_tagged_vec.begin();
	     w_tagged_iter != w_tagged_vec.end();
	     ++w_tagged_iter
	     )
	  {
	    MultiplicityTagged<u3::U3> w_tagged = (*w_tagged_iter);
	    u3::U3 w = w_tagged.irrep;
	    int rho_max = w_tagged.tag;
	    states.insert(SpanakopitaType::value_type(w,MultiplicityTagged<u3::U3>(n,rho_max)));
	  
	  }
      }

    // scan through spanakopita for subspaces

    SpanakopitaType::iterator it = states.begin();
    while (it != states.end())
      {
	// retrieve w key of this group of states
	u3::U3 w = it->first;

	// find range of this group of states
	SpanakopitaRangeType state_range = states.equal_range(w);

	// set up subspace
	//
	// Note: Creation followed by push_back means U3Subspace is
	// *copied* into vector.  
	//
	// Option 1: Eefine a lightweight
	// constructor for U3Subspace, which just saves w.  Push
	// the (empty) subspace into the Sp3RSpace.  Then use an Init
	// method to actually populate the subspace.  This runs up
	// against "private" protections in BaseSubspace (which could
	// be removed.)
	//
	// Option 2: Use C++11 emplace_back!  However, this requires
	// also modifying BaseSubspace to provide an "EmplaceSubspace"
	// with appropriate arguments, which would be difficult.
	PushSubspace(U3Subspace(w));
	subspaces_.back().Init(state_range);
	// GetSubspace(size()).Init(state_range); //ILLEGAL since GetSubspace yields const reference

	// move to start of next group of states
	it = state_range.second;
      }

  }

}  // namespace
