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

  U3Subspace::U3Subspace(const u3::U3& w, const SpanakopitaRangeType& state_range)
  {

    // set subspace labels
    labels_ = w;

    // set up indexing of states
    for (SpanakopitaType::iterator it = state_range.first; it != state_range.second; ++it)
      {
	// need to expand rho_max
	// need comparison function for MultiplicityTagged
	//PushStateLabels(it->second); 
      }

  }


  Sp3RSpace::Sp3RSpace(const u3::U3& sigma, int Nn_max)
  {

    // set space labels
    sigma_ = sigma;
    Nn_max_ = Nn_max;

    // make spanakopita
    
    SpanakopitaType states;

    // scan through spanakopita
    //   by w key group (subspace)

    SpanakopitaType::iterator it = states.begin();
    while (it != states.end())
      {
	// retrieve w key of this group of states
	u3::U3 w = it->first;

	// find range of this group of states
	SpanakopitaRangeType state_range = states.equal_range(w);

	// set up subspace
	//
	// Note: Creation followed by push-back means U3Subspace is
	// *copied* into vector.  If this overhead is troublesome, can
	// replace U3Subspace constructor with lightweight constructor
	// which just saves w.  Push the (empty) subspace into the
	// Sp3RSpace.  Then use an Init method to actually populate
	// the space once it is in the vector.

	

      }
  }

}  // namespace
