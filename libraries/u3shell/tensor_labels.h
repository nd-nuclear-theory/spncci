/****************************************************************
  tensor_labels.h

  Tensor operator labeling.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  5/15/16 (mac): Created as tensor.h.
  5/25/16 (mac): Rename to tensor_labels.h.  Extract 
    non-isospin-scheme labels as separate class
    (OperatorLabelsU3S, etc.).
  5/27/16 (mac): Restructure operator constuctor definitions.
  9/7/16 (mac): Add OperatorLabelsU3S conversion from OperatorLabelsU3ST.
  12/6/16 (aem) : Add Conjugation functions for OperatorLabelsU3ST and 
                  RelativeUnitTensorLabelsU3ST.
****************************************************************/

#ifndef TENSOR_LABELS_H_
#define TENSOR_LABELS_H_

#include <unordered_map>

#include "sp3rlib/u3.h"

namespace u3shell
{

  ////////////////////////////////////////////////////////////////
  // U3ST-scheme relative state labels
  ////////////////////////////////////////////////////////////////

  class RelativeStateLabelsU3ST
  // U(1)xSU(3)xSxT relative state labels
  //
  // Stored labels:
  //   eta (int) : relative quanta
  //   S (HalfInt) : spin (HalfInt for consistency, although value will always be integer)
  //   T (HalfInt) : isospin (HalfInt for consistency, although value will always be integer)
  //
  // Note: The SU(3) label (eta,0) is *not* explicitly stored but is
  // available through the accessor x().
  //
  // Note: The parity grade is *not* explicitly stored but is
  // available through the accessor g().
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline RelativeStateLabelsU3ST() 
      : eta_(0), S_(0), T_(0)
      // Default constructor.
      {}
    
    inline RelativeStateLabelsU3ST(int eta, HalfInt S, HalfInt T)
      : eta_(eta), S_(S), T_(T)
    // Construct from labels.
    {}
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline int eta() const
    {
      return eta_;
    }

    inline u3::SU3 x() const
    {
      return u3::SU3(eta(),0);
    }

    inline HalfInt S() const
    {
      return S_;
    }

    inline HalfInt T() const
    {
      return T_;
    }

    inline int g() const
    {
      return eta()%2;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,HalfInt,HalfInt> KeyType;
    // eta, S, T

    inline KeyType Key() const
    {
      return KeyType(eta_,S_,T_);
    }

    inline friend bool operator == (const RelativeStateLabelsU3ST& x1, const RelativeStateLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const RelativeStateLabelsU3ST& x1, const RelativeStateLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const RelativeStateLabelsU3ST& v)
    {
      boost::hash<RelativeStateLabelsU3ST::KeyType> hasher;
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
    int eta_;
    HalfInt S_, T_;
  };

  ////////////////////////////////////////////////////////////////
  // U3ST-scheme relative x CM state labels
  ////////////////////////////////////////////////////////////////

  class RelativeCMStateLabelsU3ST
  // U(1)xSU(3)xSxT relative x CM state labels
  //
  // Stored labels:
  //   eta_r (int) : relative quanta
  //   eta_cm (int) : relative quanta
  //   x (u3::SU3) : SU3 symmetry arising from (eta_r,0)x(eta_cm,0)
  //   S (HalfInt) : spin (HalfInt for consistency, although value will always be integer)
  //   T (HalfInt) : isospin (HalfInt for consistency, although value will always be integer)
  //
  // Note: The SU(3) label (eta,0) is *not* explicitly stored but is
  // available through the accessor x().
  //
  // Note: The parity grade is *not* explicitly stored but is
  // available through the accessor g().
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline RelativeCMStateLabelsU3ST() 
      : eta_r_(0),eta_cm_(0), S_(0), T_(0)
      // Default constructor.
      {}
    
    inline RelativeCMStateLabelsU3ST(int eta_r, int eta_cm, u3::SU3 x, HalfInt S, HalfInt T)
      : eta_r_(eta_r), eta_cm_(eta_cm), x_(x), S_(S), T_(T)
    // Construct from labels.
    {}
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline int eta_r() const
    {
      return eta_r_;
    }

    inline int eta_cm() const
    {
      return eta_cm_;
    }

    inline u3::SU3 x() const
    {
      return x_;
    }

    inline HalfInt S() const
    {
      return S_;
    }

    inline HalfInt T() const
    {
      return T_;
    }

    inline int N() const
    {
      return eta_cm_+eta_r_;
    }
    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,int,u3::SU3,HalfInt,HalfInt> KeyType;
    // eta, S, T

    inline KeyType Key() const
    {
      return KeyType(eta_r_,eta_cm_,x_,S_,T_);
    }

    inline friend bool operator == (const RelativeCMStateLabelsU3ST& x1, const RelativeCMStateLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const RelativeCMStateLabelsU3ST& x1, const RelativeCMStateLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const RelativeCMStateLabelsU3ST& v)
    {
      boost::hash<RelativeCMStateLabelsU3ST::KeyType> hasher;
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
    int eta_r_, eta_cm_;
    u3::SU3 x_;
    HalfInt S_, T_;
  };
  ////////



  ////////////////////////////////////////////////////////////////
  // U3ST-scheme two-body state labels
  ////////////////////////////////////////////////////////////////

  class TwoBodyStateLabelsU3ST
  // U(1)xSU(3)xSxT two-body state labels
  //
  // Note: It is presumed that the labels represent a valid two-body
  // state.  Thus, the U(1)xSU(3) labels *could* legitimately be
  // combined into a single U(3) label.  However, to minimize likely
  // arithmetic, the combined U(3) label is only constructed if needed
  // in an accessor.
  //
  // Stored labels:
  //   eta1, eta2 (int) : oscillator quanta for two particles
  //   x (u3::SU3) : SU(3) labels
  //   S (HalfInt) : spin (HalfInt for consistency, although value will always be integer)
  //   T (HalfInt) : isospin (HalfInt for consistency, although value will always be integer)
  //
  // Note: The parity grade is *not* explicitly stored but is
  // available througha the accessor g().
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline TwoBodyStateLabelsU3ST() 
      : eta1_(0), eta2_(0), S_(0), T_(0)
      // Default constructor.
      {}
    
    inline TwoBodyStateLabelsU3ST(int eta1, int eta2, const u3::SU3& x, HalfInt S, HalfInt T)
      : eta1_(eta1), eta2_(eta2), x_(x), S_(S), T_(T)
    // Construct from labels.
    {}
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline int eta1() const
    {
      return eta1_;
    }

    inline int eta2() const
    {
      return eta2_;
    }

    // inline int N() const
    // // Note: This is total quanta (0-based), not the U(1) label (3/2-based).
    // {
    //   return eta1()+eta2();
    // }

    inline u3::SU3 x() const
    {
      return x_;
    }

    inline HalfInt S() const
    {
      return S_;
    }

    inline HalfInt T() const
    {
      return T_;
    }

    inline int g() const
    {
      return (eta1()+eta2())%2;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,int,u3::SU3,HalfInt,HalfInt> KeyType;
    // eta1, eta2, x, S, T

    inline KeyType Key() const
    {
      return KeyType(eta1_,eta2_,x_,S_,T_);
    }

    inline friend bool operator == (const TwoBodyStateLabelsU3ST& x1, const TwoBodyStateLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const TwoBodyStateLabelsU3ST& x1, const TwoBodyStateLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const TwoBodyStateLabelsU3ST& v)
    {
      boost::hash<TwoBodyStateLabelsU3ST::KeyType> hasher;
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
    int eta1_, eta2_;
    u3::SU3 x_;
    HalfInt S_, T_;
  };

  ////////////////////////////////////////////////////////////////
  // U3ST-scheme generic tensor operator labels
  ////////////////////////////////////////////////////////////////

  class OperatorLabelsU3ST
  // U(1)xSU(3)xSxT operators labels
  //
  // For use in selection rules for enumerating operator sectors.
  // Meant also for use as base class for other, more complicated
  // operator labels.
  //
  // Note: The U(1)xSU(3) labels do *not* in general constitute a
  // valid U(3) label and thus cannot be stored in a u3::U3.  E.g.,
  // operators carrying N0=0 but an SU(3) irrep other than (0,0) are
  // perfectly well possible.
  //
  // Stored labels:
  //   N0 (int) : oscillator quanta
  //   x0 (u3::SU3) : SU(3) labels
  //   S0 (HalfInt) : spin (HalfInt for consistency, although value will always be integer)
  //   T0 (HalfInt) : isospin (HalfInt for consistency, although value will always be integer)
  //   g0 (int) : parity grade (0 or 1)
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline OperatorLabelsU3ST() 
      : N0_(0), S0_(0), T0_(0), g0_(0)
      // Default constructor.
      {}

    inline OperatorLabelsU3ST(int N0, const u3::SU3& x0, HalfInt S0, HalfInt T0, int g0)
      : N0_(N0), x0_(x0), S0_(S0), T0_(T0), g0_(g0)
    // Construct from labels.
    {}
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline int N0() const
    {
      return N0_;
    }

    inline u3::SU3 x0() const
    {
      return x0_;
    }

    inline HalfInt S0() const
    {
      return S0_;
    }

    inline HalfInt T0() const
    {
      return T0_;
    }

    inline int g0() const
    {
      return g0_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,u3::SU3,HalfInt,HalfInt,int> KeyType;
    // N0, x0, S0, T0, g0

    inline KeyType Key() const
    {
      return KeyType(N0_,x0_,S0_,T0_,g0_);
    }

    inline friend bool operator == (const OperatorLabelsU3ST& x1, const OperatorLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const OperatorLabelsU3ST& x1, const OperatorLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const OperatorLabelsU3ST& v)
    {
      boost::hash<OperatorLabelsU3ST::KeyType> hasher;
      return hasher(v.Key());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    protected:  // since derived class constructors may need to calculate some of these values
    int N0_;
    u3::SU3 x0_;
    HalfInt S0_, T0_;
    int g0_;
  };

  inline u3shell::OperatorLabelsU3ST Conjugate( const u3shell::OperatorLabelsU3ST& operator_labels)
    {
      return OperatorLabelsU3ST(
        -1*operator_labels.N0(),
        u3::Conjugate(operator_labels.x0()),
        operator_labels.S0(),
        operator_labels.T0(),
        operator_labels.g0()
        );
    }


  ////////////////////////////////////////////////////////////////
  // U3S-scheme generic tensor operator labels
  ////////////////////////////////////////////////////////////////

  class OperatorLabelsU3S
  // U(1)xSU(3)xS operators labels
  //
  // For use in selection rules for enumerating operator sectors.
  // Meant also for use as base class for other, more complicated
  // operator labels.
  //
  // Note: The U(1)xSU(3) labels do *not* in general constitute a
  // valid U(3) label and thus cannot be stored in a u3::U3.  E.g.,
  // operators carrying N0=0 but an SU(3) irrep other than (0,0) are
  // perfectly well possible.
  //
  // Stored labels:
  //   N0 (int) : oscillator quanta
  //   x0 (u3::SU3) : SU(3) labels
  //   S0 (HalfInt) : spin (HalfInt for consistency, although value will always be integer)
  //   g0 (int) : parity grade (0 or 1)
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline OperatorLabelsU3S() 
      : N0_(0), S0_(0), g0_(0)
      // Default constructor.
      {}

    inline OperatorLabelsU3S(int N0, const u3::SU3& x0, HalfInt S0, int g0)
      : N0_(N0), x0_(x0), S0_(S0), g0_(g0)
    // Construct from labels.
    {}
    
    inline OperatorLabelsU3S(int N0, const u3::SU3& x0, HalfInt S0)
      : N0_(N0), x0_(x0), S0_(S0)
    // Construct from labels.
    {
      g0_=N0%2;
    }

    inline OperatorLabelsU3S(const u3shell::OperatorLabelsU3ST& tensor)
    // Construct from U3ST opertator, for use in converting between PN and isospin
    // formalism and relative to many-body calculations 
      {
        N0_=tensor.N0(); 
        x0_=tensor.x0(); 
        S0_=tensor.S0();
        g0_=tensor.g0();
      }
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline int N0() const
    {
      return N0_;
    }

    inline u3::SU3 x0() const
    {
      return x0_;
    }

    inline HalfInt S0() const
    {
      return S0_;
    }

    inline int g0() const
    {
      return g0_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,u3::SU3,HalfInt,int> KeyType;
    // N0, x0, S0, T0, g0

    inline KeyType Key() const
    {
      return KeyType(N0_,x0_,S0_,g0_);
    }

    inline friend bool operator == (const OperatorLabelsU3S& x1, const OperatorLabelsU3S& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const OperatorLabelsU3S& x1, const OperatorLabelsU3S& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const OperatorLabelsU3S& v)
    {
      boost::hash<OperatorLabelsU3S::KeyType> hasher;
      return hasher(v.Key());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    protected:  // since derived class constructors may need to calculate some of these values
    int N0_;
    u3::SU3 x0_;
    HalfInt S0_;
    int g0_;
  };

  typedef std::tuple<u3shell::OperatorLabelsU3S,int,int> IndexedOperatorLabelsU3S;


  ////////////////////////////////////////////////////////////////
  // U3ST-scheme relative unit tensor operator labels
  ////////////////////////////////////////////////////////////////

  // Note: The constructor of the daughter class (RelativeUnitTensorLabelsU3ST) must either:
  //
  //   (1) set *all* base class members of the base class
  //   (OperatorLabelsU3ST) with a base class constructor call, e.g.,
  //   with the base class default copy constructor
  //   OperatorLabelsU3ST(operator_labels), or
  //
  //   (2) set base class members individually in the body of the
  //   constructor (after default initialization of the base class).
  //   Piecewise initialization x0_(x0), etc., is not recognized.


  class RelativeIndexedTensorLabelsU3ST
    : public OperatorLabelsU3ST
  // U(1)xSU(3)xSxT relative tensor operators labels
  //
  // Stored labels:
  //   (OperatorLabelsU3ST) : tensor operator labels (N0,x0,S0,T0,g0)
  //   index labels of operator : kappa0, L0
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline RelativeIndexedTensorLabelsU3ST() 
    // RelativeTensorLabelsU3STTagged
    // Default constructor.
    : kappa0_(0), L0_(0)
    {}

    inline RelativeIndexedTensorLabelsU3ST(int N0,
        const u3::SU3& x0, HalfInt S0, HalfInt T0,
        int kappa0, int L0
      )
      : OperatorLabelsU3ST(N0,x0,S0,T0,N0%2), kappa0_(kappa0), L0_(L0)
    // Construct from labels, with operator labels set individually.
    {}

    inline RelativeIndexedTensorLabelsU3ST(
        const u3shell::OperatorLabelsU3ST operator_labels,
        int kappa0, int L0       
      )
      : OperatorLabelsU3ST(operator_labels), kappa0_(kappa0), L0_(L0)
    // Construct from labels, with operator labels set collectively.
    {}
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // Note: OperatorLabelsU3ST accessors N0(), ..., are inherited
    inline int kappa0() const 
    {
      return kappa0_;
    }

    inline int L0() const
    {
      return L0_;
    }
    inline u3shell::OperatorLabelsU3ST operator_labels() const
    {
      return OperatorLabelsU3ST(N0(),x0(),S0(),T0(),g0());
    }
    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<OperatorLabelsU3ST::KeyType,int, int> KeyType;
    // operator, bra, ket

    inline KeyType Key() const
    {
      return KeyType(OperatorLabelsU3ST::Key(),kappa0(),L0());
    }

    typedef std::tuple<int,u3::SU3,HalfInt,HalfInt,int,int> FlatKeyType;
    inline FlatKeyType FlatKey() const
    {
      return FlatKeyType(N0_,x0_,S0_,T0_,kappa0_,L0_);
    }

    inline friend bool operator == (const RelativeIndexedTensorLabelsU3ST& x1, const RelativeIndexedTensorLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const RelativeIndexedTensorLabelsU3ST& x1, const RelativeIndexedTensorLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const RelativeIndexedTensorLabelsU3ST& v)
    {
      boost::hash<RelativeIndexedTensorLabelsU3ST::KeyType> hasher;
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
      int kappa0_,L0_;
  };


  class RelativeUnitTensorLabelsU3ST
    : public OperatorLabelsU3ST
  // U(1)xSU(3)xSxT relative unit tensor operators labels where unit tensor is defined by <bra|U|ket>=1
  //
  // Stored labels:
  //   (OperatorLabelsU3ST) : tensor operator labels (N0,x0,S0,T0,g0)
  //   bra, ket (RelativeStateLabelsU3ST) : bra and ket labels (eta,S,T)
  //
  // Note: Bra labels are always treated as ket-like, i.e., retrieving
  // the bra's SU(3) labels will give ket-like labels.
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline RelativeUnitTensorLabelsU3ST() 
    // Default constructor.
    {}

    inline RelativeUnitTensorLabelsU3ST(
        const u3::SU3& x0, HalfInt S0, HalfInt T0,
        const u3shell::RelativeStateLabelsU3ST& bra,
        const u3shell::RelativeStateLabelsU3ST& ket
      )
      : bra_(bra), ket_(ket)
    // Construct from labels, with operator labels set individually.
    //
    // DEPRECATED -- as less cleanly "structured" form
    //
    // Note: Purge use from upcoupling.cpp then delete?  Structured
    // form below requires N0 to be explicitly calculated, instead of
    // implicitly from bra-ket difference.
    //
    // Redundant operator labels are set from the bra/ket labels.
    {
      N0_ = bra_.eta() - ket_.eta();
      x0_= x0;
      S0_ = S0;
      T0_ = T0;
      g0_ = (bra_.g()+ket_.g())%2;  // equivalently, N0_%2
    }

    inline RelativeUnitTensorLabelsU3ST (
        const u3shell::OperatorLabelsU3ST operator_labels,
        const u3shell::RelativeStateLabelsU3ST& bra,
        const u3shell::RelativeStateLabelsU3ST& ket
      )
      : OperatorLabelsU3ST(operator_labels), bra_(bra), ket_(ket)
    // Construct from labels, with operator labels set collectively.
    //
    // The redundant (additive) operator labels are verified against
    // the bra/ket labels through assertions, but the triangularity of
    // couplings is not verified.
    {
      assert(N0_ == bra_.eta() - ket_.eta());
      assert(g0_ == (bra_.g()+ket_.g())%2);  // equivalently, N0_%2
    }
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // Note: OperatorLabelsU3ST accessors N0(), ..., are inherited

    inline const RelativeStateLabelsU3ST& bra() const
    {
      return bra_;
    }

    inline const RelativeStateLabelsU3ST& ket() const
    {
      return ket_;
    }

    inline const OperatorLabelsU3ST operator_labels()
      {
        return OperatorLabelsU3ST(N0_,x0_,S0_,T0_,g0_);
      }
    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<OperatorLabelsU3ST::KeyType,RelativeStateLabelsU3ST::KeyType,RelativeStateLabelsU3ST::KeyType> KeyType;
    // operator, bra, ket

    inline KeyType Key() const
    {
      return KeyType(OperatorLabelsU3ST::Key(),bra().Key(),ket().Key());
    }

    typedef std::tuple<u3::SU3,HalfInt,HalfInt,int,HalfInt,HalfInt,int,HalfInt,HalfInt> FlatKeyType;
    inline FlatKeyType FlatKey() const
    {
      return FlatKeyType(x0_,S0_,T0_,bra_.eta(),bra_.S(),bra_.T(),ket_.eta(),ket_.S(),ket_.T());
    }

    inline friend bool operator == (const RelativeUnitTensorLabelsU3ST& x1, const RelativeUnitTensorLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const RelativeUnitTensorLabelsU3ST& x1, const RelativeUnitTensorLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const RelativeUnitTensorLabelsU3ST& v)
    {
      boost::hash<RelativeUnitTensorLabelsU3ST::KeyType> hasher;
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
    RelativeStateLabelsU3ST bra_, ket_;
  };


  inline RelativeUnitTensorLabelsU3ST Conjugate(RelativeUnitTensorLabelsU3ST& tensor)
  {
    return RelativeUnitTensorLabelsU3ST(Conjugate(tensor.operator_labels()),tensor.ket(),tensor.bra());
  }




  class RelativeCMUnitTensorLabelsU3ST
    : public OperatorLabelsU3ST
  // U(1)xSU(3)xSxT relative unit tensor operators labels where unit tensor is defined by <bra|U|ket>=1
  //
  // Stored labels:
  //   (OperatorLabelsU3ST) : tensor operator labels (N0,x0,S0,T0,g0)
  //   bra, ket (RelativeStateLabelsU3ST) : bra and ket labels (eta,S,T)
  //
  // Note: Bra labels are always treated as ket-like, i.e., retrieving
  // the bra's SU(3) labels will give ket-like labels.
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline RelativeCMUnitTensorLabelsU3ST() 
    // Default constructor.
    : rho0_(0)
    {}

    inline RelativeCMUnitTensorLabelsU3ST(
        const u3::SU3& x0, HalfInt S0, HalfInt T0,
        int rho0,
        const RelativeCMStateLabelsU3ST& bra,
        const RelativeCMStateLabelsU3ST& ket
      )
      : bra_(bra), ket_(ket), rho0_(rho0)
    // Construct from labels, with operator labels set individually.
    //
    // DEPRECATED -- as less cleanly "structured" form
    //
    // Redundant operator labels are set from the bra/ket labels.
    {
      N0_ = bra_.N() - ket_.N();
      x0_= x0;
      S0_ = S0;
      T0_ = T0;
      g0_ = N0_%2;
    }

    inline RelativeCMUnitTensorLabelsU3ST(
        const u3shell::OperatorLabelsU3ST operator_labels,
        int rho0,
        const RelativeCMStateLabelsU3ST& bra,
        const RelativeCMStateLabelsU3ST& ket
      )
      : OperatorLabelsU3ST(operator_labels), bra_(bra), ket_(ket), rho0_(rho0)
    // Construct from labels, with operator labels set collectively.
    //
    // The redundant (additive) operator labels are verified against
    // the bra/ket labels through assertions, but the triangularity of
    // couplings is not verified.
    {
      assert(N0_ == bra_.N() - ket_.N());
      // assert(g0_ == (bra_.g()+ket_.g())%2);  // equivalently, N0_%2
    }
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // Note: OperatorLabelsU3ST accessors N0(), ..., are inherited

    inline const RelativeCMStateLabelsU3ST& bra() const
    {
      return bra_;
    }

    inline const RelativeCMStateLabelsU3ST& ket() const
    {
      return ket_;
    }
    inline int rho0() const
    {
      return rho0_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<OperatorLabelsU3ST::KeyType,int,RelativeCMStateLabelsU3ST::KeyType,RelativeCMStateLabelsU3ST::KeyType> KeyType;
    // operator, bra, ket

    inline KeyType Key() const
    {
      return KeyType(OperatorLabelsU3ST::Key(),rho0(),bra().Key(),ket().Key());
    }

    // typedef std::tuple<u3::SU3,HalfInt,HalfInt,int,HalfInt,HalfInt,int,HalfInt,HalfInt> FlatKeyType;
    // inline FlatKeyType FlatKey() const
    // {
    //   return FlatKeyType(x0_,S0_,T0_,bra_.eta(),bra_.S(),bra_.T(),ket_.eta(),ket_.S(),ket_.T());
    // }

    inline friend bool operator == (const RelativeCMUnitTensorLabelsU3ST& x1, const RelativeCMUnitTensorLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const RelativeCMUnitTensorLabelsU3ST& x1, const RelativeCMUnitTensorLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const RelativeCMUnitTensorLabelsU3ST& v)
    {
      boost::hash<RelativeCMUnitTensorLabelsU3ST::KeyType> hasher;
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
    RelativeCMStateLabelsU3ST bra_, ket_;
    int rho0_;
  };

  typedef std::unordered_map<
            u3shell::RelativeCMUnitTensorLabelsU3ST,
            double,
            boost::hash<u3shell::RelativeCMUnitTensorLabelsU3ST>
          > RelativeCMUnitTensorCache;


  ////////////////////////////////////////////////////////////////
  // U3ST-scheme two-body unit tensor operator labels
  ////////////////////////////////////////////////////////////////

  class TwoBodyUnitTensorLabelsU3ST
    : public OperatorLabelsU3ST
  // U(1)xSU(3)xSxT two-body unit tensor operators labels
  //
  // Stored labels:
  //   (OperatorLabelsU3ST) : tensor operator labels (N0,x0,S0,T0,g0)
  //   bra, ket (TwoBodyStateLabelsU3ST) : bra and ket labels (eta1,eta2,x,S,T)
  //   rho0 (int) : outer multiplicity index
  //
  // Note: Bra labels are always treated as ket-like, i.e., retrieving
  // the bra's SU(3) labels will give ket-like labels.
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline TwoBodyUnitTensorLabelsU3ST() 
      : rho0_(0)
      // Default constructor.
      {}
      
    inline TwoBodyUnitTensorLabelsU3ST(
        const u3shell::OperatorLabelsU3ST operator_labels,
        int rho0,
        const u3shell::TwoBodyStateLabelsU3ST& bra,
        const u3shell::TwoBodyStateLabelsU3ST& ket
      )
      : OperatorLabelsU3ST(operator_labels), rho0_(rho0), bra_(bra), ket_(ket)
    // Construct from labels, with operator labels set collectively.
    //
    // The redundant (additive) operator labels are verified against
    // the bra/ket labels through assertions, but the triangularity of
    // couplings is not verified.
    {
      assert(N0_ == bra_.eta1() + bra_.eta2() - ket_.eta1() - ket.eta2());
      assert(g0_ == (bra_.g()+ket_.g())%2);  // equivalently, N0_%2
    }
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // Note: OperatorLabelsU3ST accessors N0(), ..., are inherited

    int rho0() const
    {
      return rho0_;
    }

    inline const TwoBodyStateLabelsU3ST& bra() const
    {
      return bra_;
    }

    inline const TwoBodyStateLabelsU3ST& ket() const
    {
      return ket_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<
      OperatorLabelsU3ST::KeyType,
      int,
      TwoBodyStateLabelsU3ST::KeyType,
      TwoBodyStateLabelsU3ST::KeyType
      > KeyType;
    // operator, rho0, bra, ket

    inline KeyType Key() const
    {
      return KeyType(OperatorLabelsU3ST::Key(),rho0(),bra().Key(),ket().Key());
    }

    inline friend bool operator == (const TwoBodyUnitTensorLabelsU3ST& x1, const TwoBodyUnitTensorLabelsU3ST& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const TwoBodyUnitTensorLabelsU3ST& x1, const TwoBodyUnitTensorLabelsU3ST& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const TwoBodyUnitTensorLabelsU3ST& v)
    {
      boost::hash<TwoBodyUnitTensorLabelsU3ST::KeyType> hasher;
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
    TwoBodyStateLabelsU3ST bra_, ket_;
    int rho0_;
  };


  ////////////////////////////////////////////////////////////////
  // U3S-scheme two-body state labels
  ////////////////////////////////////////////////////////////////

  class TwoBodyStateLabelsU3S
  // U(1)xSU(3)xS two-body state labels, for use in pn scheme.
  //
  // See TwoBodyStateLabelsU3ST comments for further discussion.
  // Omits isospin labels.
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline TwoBodyStateLabelsU3S() 
      : eta1_(0), eta2_(0), S_(0)
      // Default constructor.
      {}
    
    inline TwoBodyStateLabelsU3S(int eta1, int eta2, const u3::SU3& x, HalfInt S)
      : eta1_(eta1), eta2_(eta2), x_(x), S_(S)
    // Construct from labels.
    {}
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline int eta1() const
    {
      return eta1_;
    }

    inline int eta2() const
    {
      return eta2_;
    }

    // inline int N() const
    // // Note: This is total quanta (0-based), not the U(1) label (3/2-based).
    // {
    //   return eta1()+eta2();
    // }

    inline u3::SU3 x() const
    {
      return x_;
    }

    inline HalfInt S() const
    {
      return S_;
    }

    inline int g() const
    {
      return (eta1()+eta2())%2;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<int,int,u3::SU3,HalfInt> KeyType;
    // eta1, eta2, x, S

    inline KeyType Key() const
    {
      return KeyType(eta1_,eta2_,x_,S_);
    }

    inline friend bool operator == (const TwoBodyStateLabelsU3S& x1, const TwoBodyStateLabelsU3S& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const TwoBodyStateLabelsU3S& x1, const TwoBodyStateLabelsU3S& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const TwoBodyStateLabelsU3S& v)
    {
      boost::hash<TwoBodyStateLabelsU3S::KeyType> hasher;
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
    int eta1_, eta2_;
    u3::SU3 x_;
    HalfInt S_;
  };

  ////////////////////////////////////////////////////////////////
  // U3S-scheme two-body unit tensor operator labels
  ////////////////////////////////////////////////////////////////

  class TwoBodyUnitTensorLabelsU3S
    : public OperatorLabelsU3S
  // U(1)xSU(3)xS two-body unit tensor operators labels, for use in pn scheme.
  //
  // See TwoBodyUnitTensorLabelsU3ST comments for further discussion.
  // Omits isospin labels.
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:

    inline TwoBodyUnitTensorLabelsU3S() 
      : rho0_(0)
      // Default constructor.
      {}

    inline TwoBodyUnitTensorLabelsU3S(
        const u3shell::OperatorLabelsU3S operator_labels,
        int rho0,
        const u3shell::TwoBodyStateLabelsU3S& bra,
        const u3shell::TwoBodyStateLabelsU3S& ket
      )
      : OperatorLabelsU3S(operator_labels), rho0_(rho0), bra_(bra), ket_(ket)
    // Construct from labels, with operator labels set collectively.
    //
    // The redundant (additive) operator labels are verified against
    // the bra/ket labels through assertions, but the triangularity of
    // couplings is not verified.
    {
      assert(N0_ == bra_.eta1() + bra_.eta2() - ket_.eta1() - ket.eta2());
      assert(g0_ == (bra_.g()+ket_.g())%2);  // equivalently, N0_%2
    }
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // Note: OperatorLabelsU3ST accessors N0(), ..., are inherited

    int rho0() const
    {
      return rho0_;
    }

    inline const TwoBodyStateLabelsU3S& bra() const
    {
      return bra_;
    }

    inline const TwoBodyStateLabelsU3S& ket() const
    {
      return ket_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<
      OperatorLabelsU3S::KeyType,
      int,
      TwoBodyStateLabelsU3S::KeyType,
      TwoBodyStateLabelsU3S::KeyType
      > KeyType;
    // operator, rho0, bra, ket

    inline KeyType Key() const
    {
      return KeyType(OperatorLabelsU3S::Key(),rho0(),bra().Key(),ket().Key());
    }

    inline friend bool operator == (const TwoBodyUnitTensorLabelsU3S& x1, const TwoBodyUnitTensorLabelsU3S& x2)
    {
      return x1.Key() == x2.Key();
    }
    
    inline friend bool operator < (const TwoBodyUnitTensorLabelsU3S& x1, const TwoBodyUnitTensorLabelsU3S& x2)
    {
      return x1.Key() < x2.Key();
    }

    inline friend std::size_t hash_value(const TwoBodyUnitTensorLabelsU3S& v)
    {
      boost::hash<TwoBodyUnitTensorLabelsU3S::KeyType> hasher;
      return hasher(v.Key());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;

    ////////////////////////////////////////////////////////////////
    // labels
    ////////////////////////////////////////////////////////////////

    protected:
    TwoBodyStateLabelsU3S bra_, ket_;
    int rho0_;
  };


}  // namespace

#endif
