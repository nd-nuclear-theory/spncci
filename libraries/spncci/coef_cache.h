/****************************************************************
  coef_cache.h

  Caching of U(3) and perhaps other coefficients.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  4/21/16 (aem,mac): Created.
****************************************************************/

#ifndef COEF_CACHE_H_
#define COEF_CACHE_H_

#include <unordered_map>

#include "sp3rlib/u3coef.h"
#include "sp3rlib/sp3r.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // U coefficient cache type
  ////////////////////////////////////////////////////////////////

  typedef std::unordered_map<
    u3::UCoefLabels,
    u3::UCoefBlock,
    boost::hash<u3::UCoefLabels> > UCoefCache;


} //namespace 

#endif
