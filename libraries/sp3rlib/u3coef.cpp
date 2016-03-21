/****************************************************************
  u3coef.cpp

  SU(3) coupling coefficient wrappers for Akiyama and Draayer su3lib.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created based on prototype u3.py and 
  T. Dytrych CSU3Master.
****************************************************************/
#include <cassert>
#include <vector>

#include "am/halfint.h"
#include "utilities/utilities.h"
#include "utilities/multiplicity_tagged.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"


namespace u3
{
  
  void U3CoefInit()
  {
    su3lib::blocks_();
  }

  double W(const u3::SU3& x1, int k1, int L1, const u3::SU3& x2, int k2, int L2, const u3::SU3& x3, int k3, int L3, int r0)
  {
    

    double w_array[su3lib::MAX_K][su3lib::MAX_K][su3lib::MAX_K][su3lib::MAX_K];
    su3lib::wu3r3w_(x1.lambda, x1.mu, x2.lambda, x2.mu, x3.lambda, x3.mu, L1 , L2, L3, r0,1,1,1, w_array);
    //Using row-major C to access column-major Fortran array
    return w_array[k3-1][k2-1][k1-1][r0-1];
  }

  double U(
    const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3,
    const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23
   )
  {
    int r12_max=u3::OuterMultiplicity(x1,x2,x12);
    int r12_3_max=u3::OuterMultiplicity(x12,x3,x);
    int r23_max=u3::OuterMultiplicity(x2,x3,x23);
    int r1_23_max=u3::OuterMultiplicity(x1,x23,x);
    int r_max=r12_max*r12_3_max*r23_max*r1_23_max;
    // std::cout<<
    // x1.lambda<<"  "<< x1.mu<<"  "<< x2.lambda<<"  "<< x2.mu<<"  "<< x.lambda<<"  "<< x.mu<<"  "<< 
    // x3.lambda<<"  "<< x3.mu<<"  "<< x12.lambda<<"  "<< x12.mu<<"  "<< x23.lambda<<"  "<< x23.mu<<"  "<<
    // r12_max<<"  "<< r12_3_max<<"  "<< r23_max<<"  "<< r1_23_max<<
    // std::endl;
    assert(
      (r_max > 0)
      &&(r12_max   >= r12)
      &&(r12_3_max >= r12_3)
      &&(r23_max   >= r23)
      &&(r1_23_max >= r1_23)
      );

    int index=r12+r12_max*(r12_3-1)+r12_max*r12_3_max*(r23-1)+r12_max*r12_3_max*r23_max*(r1_23-1)-1;
    std::vector<double> u_array;
    u_array.resize(r_max);

    su3lib::wru3optimized_(
			   x1.lambda, x1.mu, x2.lambda, x2.mu, x.lambda, x.mu, x3.lambda, x3.mu, x12.lambda, x12.mu, x23.lambda, x23.mu,
			   r12_max, r12_3_max, r23_max, r1_23_max, 
			   &u_array[0], r_max
			   );

    return u_array[index];
  }

  double Z(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x, const u3::SU3& x3, const u3::SU3& x12, int r12, int r12_3, const u3::SU3& x23, int r23, int r1_23)
  {
    int r12_max=u3::OuterMultiplicity(x1,x2,x12);
    int r12_3_max=u3::OuterMultiplicity(x12,x3,x);
    int r23_max=u3::OuterMultiplicity(x2,x2,x23);
    int r1_23_max=u3::OuterMultiplicity(x1,x23,x);
    int r_max=r12_max*r12_3_max*r23_max*r1_23_max;
    int index=r12+r12_max*(r12_3-1)+r12_max*r12_3_max*(r23-1)+r12_max*r12_3_max*r23_max*(r1_23-1)-1;
    std::vector<double> u_array;
    u_array.resize(r_max);
    su3lib::wzu3optimized_(
      x1.lambda, x1.mu, x2.lambda, x2.mu, x.lambda, x.mu, 
      x3.lambda, x3.mu, x12.lambda, x12.mu, x23.lambda, x23.mu,
      r12_max, r12_3_max, r23_max, r1_23_max, 
      &u_array[0], r_max);

    return u_array[index];
  }

  double Phi(const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x3, int r, int rp)
  // Phi phase factor that arrises in chainging the coupling order of SU(3) irreps 
  {
    return Z(x1,u3::SU3(0,0),x3,x2,x1,1,r,x2,1,rp);
  }



  double Unitary9LambdaMu(
    const u3::SU3& x1,  const u3::SU3& x2,  const u3::SU3& x12, int r12,
    const u3::SU3& x3,  const u3::SU3& x4,  const u3::SU3& x34, int r34,
    const u3::SU3& x13, const u3::SU3& x24, const u3::SU3& x,   int r13_24,
    int r13,     int r24,     int r12_34)    
    {
      int r12_max=u3::OuterMultiplicity(x1,x2,x12);
      int r34_max=u3::OuterMultiplicity(x3,x4,x34);
      int r13_max=u3::OuterMultiplicity(x1,x3,x13);
      int r24_max=u3::OuterMultiplicity(x2,x4,x24);
      int r13_24_max=u3::OuterMultiplicity(x13,x24,x);
      int r12_34_max=u3::OuterMultiplicity(x12,x34,x);
      //std::vector<double> NLM_array;
      int index=(
        r12_max*(r34-1)
        +r12_max*r34_max*(r12_34_max-1)
        +r12_max*r34_max*r12_34_max*(r13_max-1)
        +r12_max*r34_max*r12_34_max*r13_max*(r24_max-1)
        +r12_max*r34_max*r12_34_max*r13_max*r24_max*(r13_24_max-1)
        );
      int r_max=r12_max*r34_max*r13_max*r24_max*r12_34_max*r13_24_max;
      double NLM_array[r_max];
      //NLM_array.resize(r_max);
      su3lib::wu39lm_(
        x1.lambda, x1.mu, x2.lambda, x2.mu, x12.lambda, x12.mu,
        x3.lambda, x3.mu, x4.lambda, x4.mu, x34.lambda, x34.mu,
        x13.lambda, x13.mu, x24.lambda, x24.mu, x.lambda, x.mu,
        NLM_array,r_max
        );    
      return NLM_array[index];
    }


} // namespace 