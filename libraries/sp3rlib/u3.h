/****************************************************************
  u3.h

  U(3) and SU(3) labeling, branching, and Kronecker product.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created based on prototype u3states.py, u3.py, and so3.py.
  3/8/16 (aem,mac): Add U3ST structure and rename U3S structure.
  3/9/16 (aem,mac): Add KeyType typedefs.  Extract MultiplicityTagged.

****************************************************************/

#ifndef U3_H_
#define U3_H_

#include <string>
#include <utility>
#include <vector>

#include "am/halfint.h"
#include "utilities/utilities.h"
#include "utilities/multiplicity_tagged.h"

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
    // typedefs
    ////////////////////////////////////////////////////////////////

    typedef std::pair<int,int> KeyType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline SU3() 
      : lambda(0), mu(0) {}
    
    // construction from (lambda,mu)
    inline SU3(int lambda_, int mu_) 
      : lambda(lambda_), mu(mu_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      return KeyType(lambda,mu);
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

  inline double Casimir2( const u3::SU3& x)
  //Second order Casimir 
  {
    return 2./3*(sqr(x.lambda)+x.lambda*x.mu+sqr(x.mu)+3*x.lambda+3*x.mu);
  } 

  inline double Casimir3(const u3::SU3& x)
  //Third order Casimir
  {
    return 1./9*(x.lambda-x.mu)*(x.lambda+2*x.mu+3)*(2*x.lambda+x.mu+3);
  }
  
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  struct U3
  {

    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

    typedef std::pair<HalfInt,u3::SU3> KeyType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline U3() 
      : f1(0), f2(0), f3(0) {}


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

    inline KeyType Key() const
    {
      return KeyType(N(),SU3());
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

  struct U3S
  {

    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

    typedef std::pair<u3::U3,HalfInt> KeyType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline U3S() 
      : S(0) {}

    // construction from (w,S)
    inline U3S(const u3::U3& w_, const HalfInt& S_) 
      : w(w_), S(S_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const
    {
      return w;
    }

    inline u3::SU3 SU3() const
    {
      return w.SU3();
    }

    inline KeyType Key() const
    {
      return KeyType(w,S);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    u3::U3 w;
    HalfInt S;

  };


  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const U3S& wS1, const U3S& wS2)
  {
    return wS1.Key() == wS2.Key();
  }

  inline bool operator < (const U3S& wS1, const U3S& wS2)
  {
    return wS1.Key() < wS2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3S& wS)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(wS.w)*(TwiceValue(wS.S)+1);  // TODO: define dimension function for am?
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) x SU(2)x SU(2) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  struct U3ST
  {

    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

    typedef std::pair<u3::U3S,HalfInt> KeyType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline U3ST() 
      : S(0), T(0) {}

    // construction from (w,S,T)
    inline U3ST(const u3::U3& w_, const HalfInt& S_, const HalfInt& T_) 
      : w(w_), S(S_), T(T_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const
    {
      return w;
    }

    inline u3::SU3 SU3() const
    {
      return w.SU3();
    }

    inline KeyType Key() const
    {
      return KeyType(U3S(w,S),T);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    u3::U3 w;
    HalfInt S, T;

  };


  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const U3ST& wST1, const U3ST& wST2)
  {
    return wST1.Key() == wST2.Key();
  }

  inline bool operator < (const U3ST& wST1, const U3ST& wST2)
  {
    return wST1.Key() < wST2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3ST& wST)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(wST.w)*(TwiceValue(wST.S)+1)*(TwiceValue(wST.T)+1);  // TODO: define dimension function for am?
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // coupling
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // outer multiplicity

  int OuterMultiplicity(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x3);
  // Calculate outer multiplicity of SU(3) irrep of SU(3) Kronecker product.
  //
  // Calculates multiplicity of x3 in product x1 x x2.
  //
  // Reference: C. K. Chew and R. T. Sharp, Can. J. Phys. 44, 2789 (1966).  To be verified.
  // As in SU3LIB MULTU3 or UNU3SU3Basics SU3::mult.
  //
  // Arguments:
  //   x1, x2, x3 (u3::SU3): irreps
  //
  // Returns:
  //   (int) : multiplicity 

  MultiplicityTagged<u3::SU3>::vector KroneckerProduct(const u3::SU3& x1, const u3::SU3& x2);
  // Generate multiplicity-tagged vector of SU(3) irreps in SU(3) Kronecker product.
  //
  // Generates Kronecker product by iterating over possible
  // (lambda,mu) in product and checking multiplicity.
  //
  // As in UNU3SU3Basics SU3::Couple but adopting more restrictive
  // bounds on product (lambda3,mu3).
  //
  // Arguments:
  //   x1, x2 (u3::SU3) : irreps
  //
  // Returns:
  //   (MultiplicityTagged<u3::SU3>::vector) : vector with each irrep
  //   (of nonzero multiplicity) tagged by its multiplicity rho_max

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // branching
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // branching multiplicity

  int BranchingMultiplicity(const u3::SU3& x, int L);
  // Calculate branching multiplicity of angular momentum in SU(3) irrep.
  //
  // Ref: e.g., Harvey, ANP 1, 67 (1968).
  //
  // Arguments:
  //   x (u3::SU3): SU(3) irrep
  //   L (int) : angular momentum
  //
  // Returns:
  //   (int) : multiplicity 
  //
  // EX:
  //   BranchingMultiplicity(u3::SU3(4,3),3)
  //   returns 2

  MultiplicityTagged<int>::vector BranchingSO3(const u3::SU3& x);
  // Generate multiplicity-tagged vector of SO(3) irreps in SU(3) irrep.
  //
  // The general branching rule is:
  //
  //   mubar=min(lambda,mu)
  //   lambdabar=max(lambda,mu)
  //   K=mubar,mubar-2,...,1 or 0 
  //   L =
  //      K, K+1,...,lambdabar          if K!=0         
  //      lambdabar, lambdabar-2,...,1 or 0      if K=0
  //
  // The list of L values with multiplicities is, however, generated
  // by iterating over the allowed L values in this range and
  // calculating their multiplicities by BranchingMultiplicity.
  //
  // Args:
  //   x (u3::SU3) : SU(3) irrep
  //
  // Returns:
  //   (MultiplicityTagged<int>::vector) : vector with each L
  //   (of nonzero multiplicity) tagged by its multiplicity 
  //   kappa_max

  MultiplicityTagged<int>::vector BranchingSO3Constrained(const u3::SU3& x, const HalfInt::pair& r);
  // Generate multiplicity-tagged vector of SO(3) irreps in SU(3)
  // irrep, constrained to lie within a constrained angular momentum
  // range.
  //
  // The intended purpose is to allow branching only to those L values
  // which will couple with a given S to yield a given J.
  //
  // Args:
  //   x (u3::SU3) : SU(3) irrep
  //   r (HalfInt::pair) : allowed angular momentum range 
  // 
  // Although range is taken as HalfInt::pair, since it could come from
  // result of coupling J and S using am::ProductAngularMomentumRange, the actual values should be
  // integral.
  //
  // Returns:
  //   (MultiplicityTagged<int>::vector) : vector with each L
  //   (of nonzero multiplicity) tagged by its multiplicity 
  //   kappa_max

}  // namespace

#endif
