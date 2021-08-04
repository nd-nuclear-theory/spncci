/****************************************************************
  u3spn_scheme.h

  Defines generic "skeleton" bookkeeping for U(3)xSpxSnxS) subspaces,
  i.e., labeled by i.e., N(lambda,mu)SpSnS.  No detailed enumeration
  is made of the states within a subspace.  This bookkeeping is
  intended for use in keeping track of, e.g., the subspaces which
  arise within an lsu3shell SU(3) coupled basis, and operator sectors
  between these subspaces.

  Language: C++11

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  9/6/16 (mac): Created, based on a subset of u3st_scheme and u3::U3S.
  9/8/16 (mac): Add default constructors.
  2/17/16 (mac): Rename SubspaceU3SPN method Str to LabelStr.

****************************************************************/

#ifndef U3SPN_SCHEME_H_
#define U3SPN_SCHEME_H_

#include <string>

#include "basis/basis.h"
#include "basis/degenerate.h"
#include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"

namespace u3shell {

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // U(3) x Sp x Sn x S
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  class U3SPN
  {

    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:
    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline U3SPN()
      : Sp_(0), Sn_(0) {}

    // construction from (omegaS,Sp,Sn)
    inline U3SPN(const u3::U3S& omegaS, const HalfInt& Sp, const HalfInt& Sn)
      : omegaS_(omegaS), Sp_(Sp), Sn_(Sn) {}

    // alternative construction
    inline U3SPN(const u3::U3& omega, const HalfInt& Sp, const HalfInt& Sn, const HalfInt& S)
      : Sp_(Sp), Sn_(Sn)
      {
        omegaS_=u3::U3S(omega,S);
      }

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3S U3S() const
    {
      return omegaS_;
    }

    inline u3::U3 U3() const
    {
      return omegaS_.U3();
    }

    inline HalfInt N() const
    {
      return omegaS_.U3().N();
    }

    inline u3::SU3 SU3() const
    {
      return omegaS_.SU3();
    }

    inline HalfInt S() const
    {
      return omegaS_.S();
    }

    inline HalfInt Sp() const
    {
      return Sp_;
    }

    inline HalfInt Sn() const
    {
      return Sn_;
    }

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<u3::U3S,HalfInt,HalfInt> KeyType;

    inline KeyType Key() const
    {
      return KeyType(omegaS_,Sp_,Sn_);
    }

    inline friend bool operator == (const U3SPN& omegaS1, const U3SPN& omegaS2)
    {
      return omegaS1.Key() == omegaS2.Key();
    }

    inline friend bool operator < (const U3SPN& omegaS1, const U3SPN& omegaS2)
    {
      return omegaS1.Key() < omegaS2.Key();
    }

    inline friend std::size_t hash_value(const U3SPN& v)
    {
      boost::hash<U3SPN::KeyType> hasher;
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

    u3::U3S omegaS_;
    HalfInt Sp_, Sn_;

  };


  ////////////////////////////////////////////////////////////////
  // basis indexing in U3SPN scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (omega,S) = U3S
  //
  // The parity label g can implicitly be deduced from omega if needed.
  // We do not explicitly store it.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // We will not actually enumerate states within a subspace.  We
  // provide a "dummy" int state label where a StateLabelType is
  // required by the basis template library.
  //
  // However, we keep track of the subspace dimension.  The subspace
  // dimensions must be provided to the constructor via a mapping U3S
  // -> dimension.  Then we must explicitly set the dimension_ member
  // variable (inherited from BaseSubspace).
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (omega,S).
  //
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class SubspaceU3SPN
    : public basis::BaseDegenerateSubspace<SubspaceU3SPN,u3shell::U3SPN,basis::BaseState<SubspaceU3SPN>,int>
  // Subspace class for two-body states of given U(3)xS.
  //
  // SubspaceLabelsType (u3shell::U3SPN): (omega,S)
  // StateLabelsType (int): 1 (just a place holder)
  {
    public:

    // constructor

    SubspaceU3SPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SubspaceU3SPN (const u3shell::U3SPN& omegaSPN, int dimension);

    // accessors
    u3shell::U3SPN U3SPN() const {return labels();}
    u3::U3S U3S() const {return U3SPN().U3S();}
    u3::U3 U3() const {return U3SPN().U3();}
    u3::SU3 SU3() const {return U3SPN().SU3();}
    HalfInt N() const {return U3SPN().U3().N();}
    HalfInt S() const {return U3SPN().S();}
    HalfInt Sp() const {return U3SPN().Sp();}
    HalfInt Sn() const {return U3SPN().Sn();}

    // diagnostic output
    std::string LabelStr() const;
    // std::string Str() const {return LabelStr();};  // DEPRECATED

    private:
  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class SpaceU3SPN
    : public basis::BaseSpace<SubspaceU3SPN>
  // Space class for two-body states of given U(3)xS.
  {

  public:

    // constructor

    SpaceU3SPN() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceU3SPN(const std::map<u3shell::U3SPN,int>& subspace_dimensions);

    // diagnostic output
    std::string Str() const;

  private:

  };


  ////////////////////////////////////////////////////////////////
  // sectors
  ////////////////////////////////////////////////////////////////

  class SectorsU3SPN
    : public basis::BaseSectors<SpaceU3SPN>
  // U3SPN-scheme operator sectors.
  //
  // Sectors are enumerated in lexicographical order by (bra)(ket)(rho).
  {

  public:

    // constructor

    SectorsU3SPN() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SectorsU3SPN(const SpaceU3SPN& space, const OperatorLabelsU3S& operator_labels,
                             bool spin_scalar);
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).
    //
    // Arguments:
    //   space (Space3SPN): the U3SPN space
    //   operator_labels (OperatorLabelsU3S): NxSU(3)xS character of operator
    //   spin_scalar (bool): whether or not to additionally impose
    //     Sp0=0 and Sn0=0 (in which case it is assumed S0=0)

    SectorsU3SPN(
      const SpaceU3SPN& space_bra, const SpaceU3SPN& space_ket,
      const OperatorLabelsU3S& operator_labels,bool spin_scalar
    );
    // See above for description.  Allows for possibility of different bra and ket spaces
  };


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
