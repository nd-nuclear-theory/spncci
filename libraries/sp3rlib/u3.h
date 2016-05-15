/****************************************************************
  u3.h

  U(3) and SU(3) labeling, branching, and Kronecker product.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/7/16 (aem,mac): Created based on prototype u3states.py, u3.py, and so3.py.
  3/8/16 (aem,mac): Add U3ST structure and rename U3S structure.
  3/9/16 (aem,mac): Add KeyType typedefs.  Extract MultiplicityTagged.
  3/16/16 (aem): Add validity check to U(3) Kronecker product.

****************************************************************/

#ifndef U3_H_
#define U3_H_

#include <cassert>
#include <string>
#include <utility>
#include <vector>

#include <boost/functional/hash.hpp>

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

  class SU3
  {

    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::pair<int,int> KeyType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

  public:
    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline SU3() 
      : lambda_(0), mu_(0) {}
    
    // construction from (lambda,mu)
    //
    // underscore on arguments avoids name clash with accessors
    inline SU3(int lambda, int mu) 
      : lambda_(lambda), mu_(mu) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline KeyType Key() const
    {
      //assert(false);
      return KeyType(lambda(),mu());
    }

    inline int lambda() const
    {
      return lambda_;
    }

    inline int mu() const
    {
      return mu_;
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    static const int label_width = 8;
    inline friend std::size_t hash_value(const SU3& x)
    {
      int packed_labels = (x.lambda_ << label_width) | (x.mu_ << 0);

      boost::hash<int> hasher;
      return hasher(packed_labels);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////
   
  private:

    // Elliott labels
    int lambda_, mu_;

    // Could save memory by packing labels into a uint16_t:
    // uint16_t packed_labels_;
    // packed_labels_ = (lambda << label_width) | (mu << 0);

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
    return (x.lambda()+1)*(x.mu()+1)*(x.lambda()+x.mu()+2)/2;
  }

  inline u3::SU3 Conjugate(const u3::SU3& x)
  // Conjugate irrep.
  {
    return u3::SU3(x.mu(),x.lambda());
  }

  inline int ConjugationGrade(const u3::SU3& x)
  // Integer contribution to phase on conjugation.
  {
    return x.mu() + x.lambda();
  }

  inline double Casimir2( const u3::SU3& x)
  //Second order Casimir 
  {
    return 2./3*(sqr(x.lambda())+x.lambda()*x.mu()+sqr(x.mu())+3*x.lambda()+3*x.mu());
  } 

  inline double Casimir3(const u3::SU3& x)
  //Third order Casimir
  {
    return 1./9*(x.lambda()-x.mu())*(x.lambda()+2*x.mu()+3)*(2*x.lambda()+x.mu()+3);
  }
  
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  class U3
  {

    ////////////////////////////////////////////////////////////////
    // typedefs
    ////////////////////////////////////////////////////////////////

  public:
    typedef std::pair<HalfInt,u3::SU3> KeyType;
    // typedef std::tuple<HalfInt,HalfInt,HalfInt> KeyType;

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    inline U3() 
      : f1_(0), f2_(0), f3_(0)
      // default constructor
      {}


    inline U3(const HalfInt& f1, const HalfInt& f2, const HalfInt& f3)
      : f1_(f1), f2_(f2), f3_(f3)
    // Construct from Cartesian labels [f1,f2,f3].
    {
      assert(ValidLabels(f1_,f2_,f3_));
    }

    inline U3(const HalfInt& N_, const u3::SU3& x_);
    // Construct from N(lambda,mu) labels.
    //
    // Precondition: The N(lambda,mu) are assumed to be a valid U(3)
    // combination.

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
      // return (f1_ >= f2_) && (f2_ >= f3_) && ((f3_ >=0 ) || (f1_ <= 0));
      return ValidLabels(f1_,f2_,f3_);
    }

    inline static
      bool ValidLabels(const HalfInt& f1, const HalfInt& f2, const HalfInt& f3)
    // Check validity of U3 labels in Cartesian form.
    //
    // Normally there is requirement all labels nonnegative (f3>=0),
    // but we also allow conjugate representations with all labels
    // nonpositive (f1<=0).
    {
      return (f1 >= f2) && (f2 >= f3) && ((f3 >=0 ) || (f1 <= 0));
    }

    inline static
      bool ValidLabels(const HalfInt& N, const u3::SU3& x)
    // Check validity of U3 labels in N(lambda,mu) form.
    {
      int thrice_twice_f3 = TwiceValue(N-2*x.mu()-x.lambda());
      bool valid = (thrice_twice_f3%3==0);
      return valid;
    }

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // access Cartesian labels

    inline HalfInt f1() const
    {
      return f1_;
    }

    inline HalfInt f2() const
    {
      return f2_;
    }

    inline HalfInt f3() const
    {
      return f3_;
    }

    // access N and SU(3) parts
    
    // Note: Meed to use explicit reference to u3::SU3 since name is
    // masked here by u3::U3::SU3.

    inline HalfInt N() const
    {
      return f1_+f2_+f3_;
    }

    inline u3::SU3 SU3() const
    {
      int lambda = int(f1_-f2_);
      int mu = int(f2_-f3_);
      return u3::SU3(lambda,mu);
    }

    inline KeyType Key() const
    {
      return KeyType(N(),SU3());
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    static const int label_width = 12;
    inline friend std::size_t hash_value(const U3& omega)
    {
      int packed_labels =
        (TwiceValue(omega.f1_) << 2*label_width)
        | (TwiceValue(omega.f2_) << label_width)
        | (TwiceValue(omega.f3_) << 0);

      boost::hash<int> hasher;
      return hasher(packed_labels);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

  private:
    // Cartesian labels
    HalfInt f1_, f2_, f3_;

  };

  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////

  inline U3::U3(const HalfInt& N_, const u3::SU3& x_)
  {

    assert(ValidLabels(N_,x_));
      
    // recover f3 first
    // N - 2mu - lambda = (f1+f2+f3)-2*(f2-f3)-(f1-f2) = 3*f3
    // but since division is not defined for HalfInt, work with twice value for division purposes
    int twice_f3 = TwiceValue(N_-2*x_.mu()-x_.lambda()) / 3;
    f3_ = HalfInt(twice_f3,2);
      
    // recover f2 and f1
    f2_ = f3_ + x_.mu();
    f1_ = f2_ + x_.lambda();

  }

  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const U3& omega1, const U3& omega2)
  {
    return omega1.Key() == omega2.Key();
  }

  inline bool operator < (const U3& omega1, const U3& omega2)
  {
    return omega1.Key() < omega2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3& omega)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(omega.SU3());
  }

  inline u3::U3 Conjugate(const u3::U3& omega)
  // Conjugate irrep.
  {
    return u3::U3(-omega.N(),Conjugate(omega.SU3()));
  }

  inline int ConjugationGrade(const u3::U3& omega)
  // Integer contribution to phase on conjugation.
  {
    return ConjugationGrade(omega.SU3());
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) x SU(2) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // TODO (if ever used): upgrade from struct to class and add hash_value

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

    // construction from (omega,S)
    inline U3S(const u3::U3& omega_, const HalfInt& S_) 
      : omega(omega_), S(S_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const
    {
      return omega;
    }

    inline u3::SU3 SU3() const
    {
      return omega.SU3();
    }

    inline KeyType Key() const
    {
      return KeyType(omega,S);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    u3::U3 omega;
    HalfInt S;

  };


  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const U3S& omegaS1, const U3S& omegaS2)
  {
    return omegaS1.Key() == omegaS2.Key();
  }

  inline bool operator < (const U3S& omegaS1, const U3S& omegaS2)
  {
    return omegaS1.Key() < omegaS2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3S& omegaS)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(omegaS.omega)*(TwiceValue(omegaS.S)+1);  // TODO: define dimension function for am?
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) x SU(2)x SU(2) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // TODO (if ever used): upgrade from struct to class and add hash_value

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

    // construction from (omega,S,T)
    inline U3ST(const u3::U3& omega_, const HalfInt& S_, const HalfInt& T_) 
      : omega(omega_), S(S_), T(T_) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const
    {
      return omega;
    }

    inline u3::SU3 SU3() const
    {
      return omega.SU3();
    }

    inline KeyType Key() const
    {
      return KeyType(U3S(omega,S),T);
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    u3::U3 omega;
    HalfInt S, T;

  };


  ////////////////////////////////////////////////////////////////
  // relational operators
  ////////////////////////////////////////////////////////////////

  inline bool operator == (const U3ST& omegaST1, const U3ST& omegaST2)
  {
    return omegaST1.Key() == omegaST2.Key();
  }

  inline bool operator < (const U3ST& omegaST1, const U3ST& omegaST2)
  {
    return omegaST1.Key() < omegaST2.Key();
  }

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3ST& omegaST)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(omegaST.omega)*(TwiceValue(omegaST.S)+1)*(TwiceValue(omegaST.T)+1);  // TODO: define dimension function for am?
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

  MultiplicityTagged<u3::U3>::vector KroneckerProduct(const u3::U3& omega1, const u3::U3& omega2);
  // Overloaded for U3.

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
