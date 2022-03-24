/****************************************************************
  u3st.h

  U(3)xSU(2)xSU(2) labeling

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  3/24/22 (aem): Extracted from u3.h
****************************************************************/

#ifndef U3ST_H_
#define U3ST_H_

#include <cassert>
#include <string>
#include <vector>

#include "boost/functional/hash.hpp"
#include "fmt/format.h"
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "am/am.h"
#include "sp3rlib/u3.h"

namespace u3
{
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


 }  // u3 namespace


// Defining std hash functions for U3ST class. 
namespace std
{

  template<> struct hash<u3::U3ST>
  {
    inline std::size_t operator()(const u3::U3ST& h) const noexcept
    {
      return hash_value(h);
    }
  };

}


// Defining formatters for fmt format for U3 and SU3 classes. 
namespace fmt {



template <>
struct formatter<u3::U3ST> 
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
  FMT_CONSTEXPR auto format(const u3::U3ST& wST, FormatContext& ctx) -> decltype(ctx.out()) 
  {
    if(presentation == 'f')
      return format_to(ctx.out(), "{:f}{:f}{:f}",wST.U3(),wST.S(),wST.T());
    else
    return format_to(ctx.out(), "{:d}{:g}", wST.U3(),wST.S(),wST.T());
  }
};


}  // namespace fmt

#endif
