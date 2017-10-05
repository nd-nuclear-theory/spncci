/****************************************************************
  sp3rcoef.h

  Sp(3,R) labeling and branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/27/17  (aem,mac): Created. Code extracted from sp3r
****************************************************************/

#ifndef SP3RCOEF_H_
#define SP3RCOEF_H_
    
#include "sp3rlib/sp3r.h"

namespace sp3r 
{
  double CCoef(u3::U3& n, int q, u3::U3& np);
  // Coefficient of fractional parentage for expanding the raising polynomial into a
  // polynomial of a single Jacobi coordinate and a polynomial of all other coordinates 
  // To be finished later if needed

  double BCoef(u3::U3& n1, u3::U3& n2, u3::U3& n3, int rho);
  // Polynomial product coefficients 
  // To be finished later if needed 

  typedef std::tuple<u3::U3,u3::U3,u3::U3,int> BCoefLabels;
  typedef std::map<BCoefLabels,double> BCoefCache;

  void GenerateBCoefCache(BCoefCache& cache, int Nmax);
  // Generates a cache of B coefficients needed for constructing 
  // symplectic raising polynomials 

}
#endif