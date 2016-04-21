/****************************************************************
  u3coef.h

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.

  u3coef recquires initilization of the coeffcients by calling the function 
  u3::U3CoefInit() at the start of main

                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.


****************************************************************/

#ifndef U3COEF_H_
#define U3COEF_H_

#include <cassert>
#include <map>
#include <tuple>

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

  class UCoefLabel
  {
  public:
    typedef std::tuple<u3::SU3,u3::SU3,u3::SU3,u3::SU3,u3::SU3,u3::SU3> KeyType;

  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  inline UCoefLabel(const u3::SU3 x1, const u3::SU3 x2, const u3::SU3 x,
   const u3::SU3 x3, const u3::SU3 x12, const u3::SU3 x23)
      :x1_(x1), x2_(x2), x_(x), x3_(x3), x12_(x12), x23_(x23){}
  ////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////
 
  inline KeyType Key() const
  {
    return KeyType(x1_, x2_, x_, x3_, x12_, x23_);
  }

  ////////////////////////////////////////////////////////////////
  // hashing
  ////////////////////////////////////////////////////////////////

  inline friend std::size_t hash_value(const UCoefLabel& ucoef_labels)
  {
    // TODO (Andika)

    return 0;
  }
  ////////////////////////////////////////////////////////////////
  // string conversion
  ////////////////////////////////////////////////////////////////
  
  std::string Str() const;

  ////////////////////////////////////////////////////////////////
  // labels
  ////////////////////////////////////////////////////////////////

  // Operator labels
  u3::SU3 x1_, x2_, x_, x3_, x12_, x23_;

  };

  void U3CoefInit();


  double W(const u3::SU3& x1, int k1, int L1, const u3::SU3& x2, int k2, int L2, const u3::SU3& x3, int k3, int L3, int r0);
  // SU(3) reduced coupling coefficient, as referred to as Wigner coefficient 
  //
  // Provides wrapper for su3lib function wu3r3w_

  std::vector<double> UCoefVector(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3,
    const u3::SU3& x12, int r12_max, int r12_3_max, const u3::SU3& x23, int r23_max, int r1_23_max);


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

  double Phi(const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x3, int r, int rp);
  // Phi phase factor that arrises in chainging the coupling order of SU(3) irreps 

} //namespace 






#endif
