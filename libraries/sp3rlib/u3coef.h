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

#include <vector>

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
      extern void wu39lm_(const int&, const int& , const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, const int&, const int&, const int& , const int& , const int& , const int&, const int&, const int&, double[], const int&);
      extern void blocks_(void);
    }

  } //namespace

  inline void blocks()
  {
    su3lib::blocks_();
  } 

  double W(const u3::SU3& x1, int k1, int L1, const u3::SU3& x2, int k2, int L2, const u3::SU3& x3, int k3, int L3, int r0);
  // SU(3) reduced coupling coefficient, as referred to as Wigner coefficient 
  //
  // Provides wrapper for su3lib function wu3r3w_



  double U(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3, const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23);
  // SU(3) Racah recoupling coefficient for recoupling from (1x2)x3 to 1x(2x3). 
  //
  // Provides wrapper for su3lib function wru3optimized


  double Z(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3, const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23);
  // SU(3) Racah recoupling coefficient for recoupling from (1x2)x3 to 2x(1x3). 
  //
  // Provides wrapper for su3lib function wzu3optimized

  double Unitary9LambdaMu(
    const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x12, int r12,
    const u3::SU3& x3,  const u3::SU3& x4,  const u3::SU3& x34, int r34,
    const u3::SU3& x13, const u3::SU3& x24, const u3::SU3& x,   int r13_24,
    int r13,     int r24,     int r12_34
    );
  //SU(3) unitary 9-(lambda,mu) symbol.
  //
  // Provides wrapper for su3lib function wu39lm_


} //namespace 






#endif
