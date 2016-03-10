/****************************************************************
  u3coef.h

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.
****************************************************************/

#ifndef U3COEF_H_
#define U3COEF_H_

#include "am/halfint.h"
#include "utilities/utilities.h"
#include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"


namespace u3
{
	namespace su3lib
	{

		const size_t MAX_K = 9;

		// Subroutines of original Fortran SU(3) library
		extern "C" 
		{ 
			extern void wu3r3w_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[MAX_K][MAX_K][MAX_K][MAX_K]);
			extern void wru3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
			extern void wzu3optimized_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&);
			//extern void xewu3_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], int[], int[], int[], const int&, const int&, const int&);
			//extern void xwu3_(const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, const int&, double[], const int&, const int&, double[], int[], int[], int[], int[], const int&, const int&, int[], const int&, const int&, const int&, const int&);
			extern void wu39lm_(const int&, const int& , const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, double[], const int&);
			extern void blocks_(void);
		}

	} //namespace






	inline double U(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3, const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23)
	{
		int r12max=u3::OuterMultiplicity(x1,x2,x12);
		int r12_3max=u3::OuterMultiplicity(x12,x3,x);
		int r23max=u3::OuterMultiplicity(x2,x2,x23);
		int r1_23max=u3::OuterMultiplicity(x1,x23,x);
		int rmax=r12max*r12_3max*r23max*r1_23max;
		int index=r12+r12max*(r12_3-1)+r12max*r12_3max*(r23-1)+r12max*r12_3max*r23max*(r1_23-1)-1;
		double* u_array = static_cast<double*>(malloc(sizeof(double)*rmax));
		su3lib::wru3optimized_(x1.lambda, x1.mu, x2.lambda, x2.mu, x.lambda, x.mu, x3.lambda, x3.mu, x12.lambda, x12.mu, r12, r12_3, x23.lambda, x23.mu, r23, r1_23, u_array, rmax);
		double ucoef=u_array[index];
		free(u_array);
		return ucoef;
	}
} //namespace 


//coef=su3lib.wru3optimized(w1.su3.lbda,w1.su3.mu,w2.su3.lbda,w2.su3.mu,w.su3.lbda,w.su3.mu,w3.su3.lbda,w3.su3.mu,w12.su3.lbda,w12.su3.mu,w23.su3.lbda,w23.su3.mu,r12max,r12_3max,r23max,r1_23max,rmax)
//Coef=coef[index]



#endif



