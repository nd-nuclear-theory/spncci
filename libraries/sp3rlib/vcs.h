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
#include "sp3rlib/u3coef.h"  
  


namespace vcs
{
	inline double Omega(const u3::U3& n, const u3::U3& w)

	// Omega returns the Omega factor used in Kmatrix calculations
	//
	// Based on protopye vcs.py and equation give in
	//   D. J. Rowe, J. Math Phys. 25 (1984) 2662. 
	//
	// Returns: (double) omega 
	{
		double omega=0;
		omega += (2*sqr(double(w.f1))-sqr(double(n.f1))+8*(double(w.f1)-double(n.f1))-2*(2*double(w.f1)-double(n.f1)));
		omega += (2*sqr(double(w.f2))-sqr(double(n.f2))+8*(double(w.f2)-double(n.f2))-4*(2*double(w.f2)-double(n.f2)));
		omega += (2*sqr(double(w.f3))-sqr(double(n.f3))+8*(double(w.f3)-double(n.f3))-6*(2*double(w.f3)-double(n.f3)));
		return omega/4.;

	}

	double BosonCreationRME(const u3::U3& np, const u3::U3& n);
	// SU(3) Reduced matrix element of a^\dagger boson creation operator
	// 
	// Based on protoype u3boson.py  Formula is given by:
	//   G. Rosensteel and D. J. Rowe. J. Math Phys. 24 (1983) 2461. 
	//
	// Returns:
	//    rme: (double) reduced matrix element of boson creation operator. 

	double SMatrix(const u3::U3& s, const u3::U3& w, MultiplicityTagged<u3::U3>& nr1,MultiplicityTagged<u3::U3>& nr2);
	// Calculate the K^2 matrix elements

	Eigen::MatrixXd KMatrix(const u3::U3& s, const u3::U3& w);
	//Calculates the K matrix 	

}  //  namespace

#endif
