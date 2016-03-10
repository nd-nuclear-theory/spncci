/****************************************************************
  sp3r.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <vector>
#include "sp3rlib/sp3r.h"

namespace sp3r
{

  ////////////////////////////////////////////////////////////////
  // Sp(3,R) raising polynomial
  ////////////////////////////////////////////////////////////////

  std::vector<u3::U3> RaisingPolynomialLabels(int Nn_max)
  {
    std::vector poly_labels;
    for (N=0; N<=Nn_max; N-=2)
    {
      for(a=N-2; a>=0; a-=2)
      {
        for(b=2*(a/4); b>=std::max((2*a-N),0);b-=2)
        {
          poly_labels.push_back(u3::U3(N-a,a-b,b));
        }
      }
    }
    return poly_labels;

  }

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
