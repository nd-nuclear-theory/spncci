/****************************************************************
  u3.h

  U(3) and SU(3) labeling, branching, and Kronecker product.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created based on prototype u3states.py, u3.py, and so3.py.
  3/8/16 (aem,mac): Add U3ST structure and rename U3S structure.
  3/9/16 (aem,mac): Add KeyType typedefs.  Extract MultiplicityTagged.

****************************************************************/

#ifndef MOSHINSKY_H_
#define MOSHINSKY_H_

#include "am/halfint.h"
#include "u3.h"
//#include "utilities/utilities.h"
//#include "utilities/multiplicity_tagged.h"

namespace u3 
{
	double MoshinskyCoefficient(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& xr,const u3::SU3& xc,const u3::SU3& xr);
	// SU(3) Moshinsky Coefficient which is equivalent to a Wigner little d function evaluated at pi/2

	double MoshinskyCoefficient(const u3::U3& w1, const u3::U3& w2, const u3::U3& wr,const u3::U3& wc,const u3::U3& wr);
	//Overleading for U3 arguements 

	double MoshinskyCoefficient(int r1, int r2, r, R, const u3::U3& wr);
	//Overloading Moshinsky to take integer arguements for two-body and relative-center of mass arguments
	// and SU(3) total symmetry (lambda,mu)

	double MoshinskyCoefficient(int r1, int r2, int r, int R, const u3::U3& w)
	// Overloading Moshinsky to take integers and U3 for total symmetry


} //namespace









