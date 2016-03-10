/****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy on 3/9/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py.

****************************************************************/

#ifndef VCS_H_
#define VCS_H_

#include <cmath>
#include "sp3rlib/u3.h"
  


namespace vcs
{
	inline double Omega(const u3::U3& n, const u3::U3& w)

	// Omega returns the Omega factor used in Kmatrix calculations
	{
		double omega;
		omega += (2*sqr(w.f1)-sqr(n.f1)+8*(w.f1-n.f1)-2*(2*w.f1-n.f1));
		omega += (2*sqr(w.f2)-sqr(n.f2)+8*(w.f2-n.f2)-4*(2*w.f2-n.f2));
		omega += (2*sqr(w.f3)-sqr(n.f3)+8*(w.f3-n.f3)-6*(2*w.f3-n.f3));
		return omega;

	}

	double BosonCreationRME(const u3::U3& np, const u3::U3& n)
	// SU(3) Reduced matrix element of a^\dagger boson creation operator
	// 
	// Returns:
	//    rme: double reduced matrix element of boson creation operator. 

}  //  namespace

#endif