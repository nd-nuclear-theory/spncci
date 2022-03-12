/****************************************************************
  u3.h

  U(3) and SU(3) labeling, branching, and Kronecker product.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/7/16 (aem,mac): Created based on prototype u3states.py, u3.py, and so3.py.
  3/8/16 (aem,mac): Add U3ST structure and rename U3S structure.
  3/9/16 (aem,mac): Add KeyType typedefs.  Extract MultiplicityTagged.
  3/16/16 (aem): Add validity check to U(3) Kronecker product.
  9/6/16 (mac): Upgrade U3S and U3ST from struct to class with hash function, etc.
  09/24/21 (pjf):
    - Add constexpr to SU3 and U3 where applicable.
    - Fill in missing comparison operators for SU3 and U3.
    - Add U3 overload for OuterMultiplicity.
  11/1/21 (aem): 
    - Changed switched U3 to storing N,lambda,mu instead of f1,f2,f3
    - Changed SU3 lambda,mu type to unsigned int
    - Changed multiplicity type to unsigned int
    - Defined fmt formats for SU3, U3 and U3S
****************************************************************/

#ifndef U3_H_
#define U3_H_

#include <cassert>
#include <string>
#include <vector>

#include "boost/functional/hash.hpp"
#include "fmt/format.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "am/am.h"
#include "sp3rlib/multiplicity_tagged.h"
#include "mcutils/arithmetic.h"
#include "mcutils/deprecated.h"

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
    // constructors
    ////////////////////////////////////////////////////////////////

    public:
    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    SU3()=default;
    // constexpr inline SU3()
    //   : lambda_(0), mu_(0) {}

    // construction from (lambda,mu)
    //
    // underscore on arguments avoids name clash with accessors
    constexpr inline SU3(unsigned int lambda, unsigned int mu)
      : lambda_(lambda), mu_(mu) {}

    DEPRECATED("use unsigned int for lambda and mu")
    constexpr inline SU3(int lambda, int mu)
      : SU3((unsigned int)lambda, (unsigned int)mu)
    {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    constexpr inline int lambda() const
    {
      return int(lambda_);
    }

    constexpr inline int mu() const
    {
      return int(mu_);
    }

    constexpr inline unsigned int lambda_u() const
    {
      return lambda_;
    }

    constexpr inline unsigned int mu_u() const
    {
      return mu_;
    }


    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::pair<unsigned int,unsigned int> KeyType;

    constexpr inline KeyType Key() const
    {
      return KeyType(lambda(),mu());
    }

    constexpr inline friend bool operator==(const SU3& x1, const SU3& x2)
    {
      return x1.Key() == x2.Key();
    }

    constexpr inline friend bool operator!=(const SU3& x1, const SU3& x2)
    {
      return !(x1 == x2);
    }

    constexpr inline friend bool operator<(const SU3& x1, const SU3& x2)
    {
      return x1.Key() < x2.Key();
    }

    constexpr inline friend bool operator>(const SU3& x1, const SU3& x2)
    {
      return x2 < x1;
    }

    constexpr inline friend bool operator<=(const SU3& x1, const SU3& x2)
    {
      return !(x1 > x2);
    }

    constexpr inline friend bool operator>=(const SU3& x1, const SU3& x2)
    {
      return !(x1 < x2);
    }

    // alternative: if find need to avoid hash combination functions...
    //
    //   static const int label_width = 8;
    //   int packed_labels = (x.lambda_ << label_width) | (x.mu_ << 0);
    //   boost::hash<int> hasher;
    //   return hasher(packed_labels);

    inline friend std::size_t hash_value(const SU3& v)
    {
      boost::hash<SU3::KeyType> hasher;
      return hasher(v.Key());
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
    unsigned int lambda_, mu_;

  };

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  constexpr inline int dim(const u3::SU3& x)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return (x.lambda()+1)*(x.mu()+1)*(x.lambda()+x.mu()+2)/2;
  }

  constexpr inline u3::SU3 Conjugate(const u3::SU3& x)
  // Conjugate irrep.
  {
    return u3::SU3(x.mu_u(),x.lambda_u());
  }

  constexpr inline int ConjugationGrade(const u3::SU3& x)
  // Integer contribution to phase on conjugation.
  {
    return x.mu() + x.lambda();
  }

  constexpr inline double Casimir2( const u3::SU3& x)
  // Second order Casimir
  {
    return 2./3*((x.lambda()*x.lambda())+x.lambda()*x.mu()+(x.mu()*x.mu())+3*x.lambda()+3*x.mu());
  }

  constexpr inline double Casimir3(const u3::SU3& x)
  // Third order Casimir
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
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    // copy constructor: synthesized copy constructor since only data
    // member needs copying
      U3()=default;
    // constexpr inline U3()
    //   : N_{0},x_{}
    //   // default constructor
    //   {}


    constexpr inline U3(const HalfInt& f1, const HalfInt& f2, const HalfInt& f3)
      : N_{f1+f2+f3}, x_{(unsigned int)(f1-f2), (unsigned int)(f2-f3)}
    // Construct from Cartesian labels [f1,f2,f3].
    {
      
      assert(ValidLabels(f1,f2,f3));
    }

    constexpr inline U3(const HalfInt& N, const u3::SU3& x)
      : N_(N),x_(x)
    // Construct from N(lambda,mu) labels.
    {
      assert(ValidLabels(N,x));
    }

    ////////////////////////////////////////////////////////////////
    // validation
    ////////////////////////////////////////////////////////////////
    constexpr inline static
      bool ValidLabels(const HalfInt& f1, const HalfInt& f2, const HalfInt& f3)
    // Check validity of U3 labels in Cartesian form.
    {

      return ((f1 >= f2) && (f2 >= f3) && ((f3 >=0 )||(f1<=0))) && IsInteger(f1-f2) && IsInteger(f2-f3);
    }


    constexpr inline static
      bool ValidLabels(const HalfInt& N, const u3::SU3& x)
    // Check validity of U3 labels in N(lambda,mu) form.
    {
      bool valid = TwiceValue(N-int(x.lambda()+2*x.mu()))%3==0;
      // valid &= TwiceValue(N-x.lambda()+x.mu())%3==0;
      // valid &= TwiceValue(N-x.lambda()-2*x.mu())%3 ==0;
      return valid;
    }


    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // access Cartesian labels
    constexpr inline HalfInt N() const {return N_;}
    constexpr inline u3::SU3 SU3() const{return x_;}
   
   // but since division is not defined for HalfInt, work with twice value for division purposes
    constexpr inline HalfInt f1() const
    {
      return HalfInt(TwiceValue(N_+2*x_.lambda()+x_.mu())/3,2);
    }

    constexpr inline HalfInt f2() const
    {
      return HalfInt(TwiceValue(N_-x_.lambda()+x_.mu())/3,2);
    }

    constexpr inline HalfInt f3() const
    {
      return HalfInt(TwiceValue(N_-x_.lambda()-2*x_.mu())/3,2);
    }

    constexpr inline std::tuple<HalfInt,HalfInt,HalfInt> f() const
      {
        return {f1(),f2(),f3()};
      }
    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::pair<HalfInt,u3::SU3> KeyType;

    constexpr inline KeyType Key() const
    {
      return KeyType(N_,x_);
    }

    constexpr inline friend bool operator==(const U3& omega1, const U3& omega2)
    {
      return omega1.Key() == omega2.Key();
    }

    constexpr inline friend bool operator!=(const U3& omega1, const U3& omega2)
    {
      return !(omega1 == omega2);
    }

    constexpr inline friend bool operator<(const U3& omega1, const U3& omega2)
    {
      return omega1.Key() < omega2.Key();
    }

    constexpr inline friend bool operator>(const U3& omega1, const U3& omega2)
    {
      return omega2 < omega1;
    }

    constexpr inline friend bool operator<=(const U3& omega1, const U3& omega2)
    {
      return !(omega1 > omega2);
    }

    constexpr inline friend bool operator>=(const U3& omega1, const U3& omega2)
    {
      return !(omega1 < omega2);
    }

    inline friend std::size_t hash_value(const U3& v)
    {
      boost::hash<U3::KeyType> hasher;
      return hasher(v.Key());
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
    // HalfInt f1_, f2_, f3_;
    HalfInt N_;
    u3::SU3 x_;

  };


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

  class U3S
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:
    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline U3S()
      : S_(0) {}

    // construction from (omega,S)
    inline U3S(const u3::U3& omega, const HalfInt& S)
      : omega_(omega), S_(S) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const
    {
      return omega_;
    }

    inline u3::SU3 SU3() const
    {
      return omega_.SU3();
    }

    inline HalfInt S() const
    {
      return S_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<u3::U3,HalfInt> KeyType;

    inline KeyType Key() const
    {
      return KeyType(omega_,S_);
    }

    inline friend bool operator == (const U3S& omegaS1, const U3S& omegaS2)
    {
      return omegaS1.Key() == omegaS2.Key();
    }

    inline friend bool operator < (const U3S& omegaS1, const U3S& omegaS2)
    {
      return omegaS1.Key() < omegaS2.Key();
    }

    inline friend std::size_t hash_value(const U3S& v)
    {
      boost::hash<U3S::KeyType> hasher;
      return hasher(v.Key());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////

    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    private:

    u3::U3 omega_;
    HalfInt S_;

  };

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3S& omegaS)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(omegaS.U3())*am::dim(omegaS.S());
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) x SU(2) x SU(2) irrep
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  class U3ST
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:
    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline U3ST()
      : S_(0), T_(0) {}

    // construction from (omega,S,T)
    inline U3ST(const u3::U3& omega, const HalfInt& S, const HalfInt& T)
      : omega_(omega), S_(S), T_(T) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const
    {
      return omega_;
    }

    inline u3::SU3 SU3() const
    {
      return omega_.SU3();
    }

    inline HalfInt S() const
    {
      return S_;
    }

    inline HalfInt T() const
    {
      return T_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<u3::U3,HalfInt,HalfInt> KeyType;

    inline KeyType Key() const
    {
      return KeyType(omega_,S_,T_);
    }

    inline friend bool operator == (const U3ST& omegaST1, const U3ST& omegaST2)
    {
      return omegaST1.Key() == omegaST2.Key();
    }

    inline friend bool operator < (const U3ST& omegaST1, const U3ST& omegaST2)
    {
      return omegaST1.Key() < omegaST2.Key();
    }

    inline friend std::size_t hash_value(const U3ST& v)
    {
      boost::hash<U3ST::KeyType> hasher;
      return hasher(v.Key());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////

    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    private:

    u3::U3 omega_;
    HalfInt S_,T_;

  };

  ////////////////////////////////////////////////////////////////
  // group theory functions
  ////////////////////////////////////////////////////////////////

  inline int dim(const u3::U3ST& omegaST)
  // Calculate dimension of irrep.
  //
  // Note: Use lowercase abbreviated form "dim" to match mathematical notation.
  {
    return dim(omegaST.U3())*am::dim(omegaST.S())*am::dim(omegaST.T());
  }

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // coupling
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  // outer multiplicity

  unsigned int OuterMultiplicity(const u3::SU3& x1, const u3::SU3& x2, const u3::SU3& x3);
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

  inline unsigned int OuterMultiplicity(const u3::U3& omega1, const u3::U3& omega2, const u3::U3& omega3)
  // Overloaded for U3.
  {
    if (omega1.N() + omega2.N() != omega3.N())
      return 0;
    return OuterMultiplicity(omega1.SU3(), omega2.SU3(), omega3.SU3());
  }

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

  unsigned int BranchingMultiplicitySO3(const u3::SU3& x, int L);
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

  MultiplicityTagged<unsigned int>::vector BranchingSO3(const u3::SU3& x);
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

  MultiplicityTagged<unsigned int>::vector BranchingSO3Constrained(const u3::SU3& x, const HalfInt::pair& r);
  // DEPRECATED
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

namespace std
{
  template<> struct hash<u3::U3>
  {
    inline std::size_t operator()(const u3::U3& h) const
    {
      return hash_value(h);
    }
  };

  template<> struct hash<u3::U3S>
  {
    inline std::size_t operator()(const u3::U3S& h) const noexcept
    {
      return hash_value(h);
    }
  };


// template<> struct hash<typename lgi::LGI::UpstreamLabelsType>
//   {
//     inline std::size_t operator()(const lgi::LGI::UpstreamLabelsType& h) const noexcept
//     {
//       return boost::hash<lgi::LGI::UpstreamLabelsType>{}(h);
//     }
//   };


}


namespace fmt {

template <>
struct formatter<u3::SU3> 
{
  char presentation = 'g';

  template <typename ParseContext>
  FMT_CONSTEXPR auto parse(ParseContext& ctx) -> decltype(ctx.begin()) 
  {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'd' || *it == 'g')) presentation = *it++;

    // Check if reached the end of the range:
    if (it != end && *it != '}')
      throw format_error("invalid format");

    // Return an iterator past the end of the parsed range:
    return it;
  }

  template <typename FormatContext>
  FMT_CONSTEXPR auto format(const u3::SU3& x, FormatContext& ctx) -> decltype(ctx.out()) 
  {
    return format_to(ctx.out(), "({:d},{:d})", x.lambda(),x.mu());
  }
};

template <>
struct formatter<u3::U3> 
{
  char presentation = 'g';

  template <typename ParseContext>
  FMT_CONSTEXPR auto parse(ParseContext& ctx) -> decltype(ctx.begin()) 
  {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'd' || *it == 'g' || *it == 'f')) presentation = *it++;

    // Check if reached the end of the range:
    if (it != end && *it != '}')
      throw format_error("invalid format");

    // Return an iterator past the end of the parsed range:
    return it;
  }

  template <typename FormatContext>
  FMT_CONSTEXPR auto format(const u3::U3& w, FormatContext& ctx) -> decltype(ctx.out()) 
  {
    if(presentation == 'f')
      return format_to(ctx.out(), "{:f}{:d}",w.N(),w.SU3());
    else
    return format_to(ctx.out(), "{:g}{:d}", w.N(),w.SU3());
  }
};


template <>
struct formatter<u3::U3S> 
{
  char presentation = 'g';

  template <typename ParseContext>
  FMT_CONSTEXPR auto parse(ParseContext& ctx) -> decltype(ctx.begin()) 
  {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'd' || *it == 'g' || *it == 'f')) presentation = *it++;

    // Check if reached the end of the range:
    if (it != end && *it != '}')
      throw format_error("invalid format");

    // Return an iterator past the end of the parsed range:
    return it;
  }

  template <typename FormatContext>
  FMT_CONSTEXPR auto format(const u3::U3S& wS, FormatContext& ctx) -> decltype(ctx.out()) 
  {
    if(presentation == 'f')
      return format_to(ctx.out(), "{:f}{:f}",wS.U3(),wS.S());
    else
    return format_to(ctx.out(), "{:d}{:g}", ws.U3(),wS.S());
  }
};


}  // namespace fmt

#endif
