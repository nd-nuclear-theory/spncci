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

  U3Subspace::U3Subspace(const u3::U3& w, const SpanakopitaType& states)
  {

    // set subspace labels
    labels_ = w;

    // set up indexing of states
    SpanakopitaRangeType state_range = states.equal_range(w);
    for (SpanakopitaType::iterator it = state_range.first; it != state_range.second; ++it)
      PushStateLabels(it->second); 

  }

}  // namespace
