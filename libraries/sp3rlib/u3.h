/****************************************************************
  u3.h

  U(3) and SU(3) labeling, branching, and Kronecker product.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac) : Created based on prototype u3states.py, u3.py, and so3.py.

****************************************************************/

#ifndef U3_H_
#define U3_H_

#include <string>
#include <utility>

#include "am/halfint.h"

namespace u3 
{

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // SU(3) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  struct SU3
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor: syntesized default constructor (?),
    // presumably leaves data uninitialized (?)

    // construction from (lambda,mu)
    inline SU3(int lambda_, int mu_) 
      : lambda(lambda_), mu(mu_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline std::pair<int,int> Key() const
    {
      return std::pair<int,int>(lambda,mu);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Elliott labels
    int lambda, mu;

  };

  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const SU3& x1, const SU3& x2)
  {
    return x1.Key() == x2.Key();
  }

  inline bool operator < (const SU3& x1, const SU3& x2)
  {
    return x1.Key() < x2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::SU3& x)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return (x.lambda+1)*(x.mu+1)*(x.lambda+x.mu+2)/2;
  }

  inline u3::SU3 Conjugate(const u3::SU3& x)
  // Conjugate irrep.
  {
    return u3::SU3(x.mu,x.lambda);
  }

  inline int ConjugationGrade(const u3::SU3& x)
  // Integer contribution to phase on conjugation.
  {
    return x.mu + x.lambda;
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  struct U3
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor: syntesized default constructor (?),
    // presumably leaves data uninitialized (?)

    // construction from f1,f2,f3
    inline U3(const HalfInt& f1_, const HalfInt& f2_, const HalfInt& f3_) 
      : f1(f1_), f2(f2_), f3(f3_) {}

    // construction from N and lm
    inline U3(const HalfInt& N_, const u3::SU3& lm_);

    ////////////////////////////////////////////////////////////////
    // validation
    ////////////////////////////////////////////////////////////////

    inline bool Valid() const
    // Checks validity of U3 labels.
    //
    // Normally there is requirement all labels nonnegative (f3>=0),
    // but we also allow conjugate representations with all labels
    // nonpositive (f1<=0).
    {
      return (f1 >= f2) && (f2 >= f3)
	&& ((f3>=0) || (f1<=0));
    }

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // access N and SU(3) part
    
    // Note: Meed to use explicit reference to u3::SU3 since name is
    // masked here by u3::U3::SU3.

    inline HalfInt N() const
    {
      return f1+f2+f3;
    }

    inline u3::SU3 SU3() const
    {
      int lambda = int(f1-f2);
      int mu = int(f2-f3);
      return u3::SU3(lambda,mu);
    }

    inline std::pair<HalfInt,u3::SU3> Key() const
    {
      return std::pair<HalfInt,u3::SU3>(N(),SU3());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Cartesian labels
    HalfInt f1, f2, f3;

  };

  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  inline U3::U3(const HalfInt& N_, const u3::SU3& lm_)
  {

      // recover f3 first
      // N - 2mu - lambda = (f1+f2+f3)-2*(f2-f3)-(f1-f2) = 3*f3
      // but since division is not defined for HalfInt, work with twice value for division purposes
      int twice_f3 = TwiceValue(N_-2*lm_.mu-lm_.lambda) / 3;
      f3 = HalfInt(twice_f3,2);
      
      // recover f2 and f1
      f2 = f3 + lm_.mu;
      f1 = f2 + lm_.lambda;

    }


  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const U3& w1, const U3& w2)
  {
    return w1.Key() == w2.Key();
  }

  inline bool operator < (const U3& w1, const U3& w2)
  {
    return w1.Key() < w2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3& w)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(w.SU3());
  }

  inline u3::U3 Conjugate(const u3::U3& w)
  // Conjugate irrep.
  {
    return u3::U3(-w.N(),Conjugate(w.SU3()));
  }

  inline int ConjugationGrade(const u3::U3& w)
  // Integer contribution to phase on conjugation.
  {
    return ConjugationGrade(w.SU3());
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) x SU(2) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  struct U3SU2
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor: syntesized default constructor (?),
    // presumably leaves data uninitialized (?)

    // construction from (w,S)
    inline U3SU2(const u3::U3& w_, const HalfInt& S_) 
      : w(w_), S(S_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const
    {
      return w;
    }

    inline std::pair<u3::U3,HalfInt> Key() const
    {
      return std::pair<u3::U3,HalfInt>(w,S);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    // Cartesian labels
    u3::U3 w;
    HalfInt S;

  };


  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const U3SU2& wS1, const U3SU2& wS2)
  {
    return wS1.Key() == wS2.Key();
  }

  inline bool operator < (const U3SU2& wS1, const U3SU2& wS2)
  {
    return wS1.Key() < wS2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3SU2& wS)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(wS.w)*(TwiceValue(wS.S)+1);  // TODO: define dimension function for am?
  }



}  // namespace

#endif
