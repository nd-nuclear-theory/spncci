/****************************************************************
  u3s.h

  U(3)xSU(2) and SU(3)xSU(2) labeling

  Anna E. McCoy
  Institute for Nuclear Theory 

  SPDX-License-Identifier: MIT

  3/24/22 (aem): Extracted from u3.h.
****************************************************************/

#ifndef U3S_H_
#define U3S_H_

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
  class U3S
  {
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    public:
    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    // default constructor
    inline U3S()=default;

    // construction from (omega,S)
    inline U3S(const u3::U3& omega, const HalfInt& S)
      : omega_(omega), S_(S) {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    inline u3::U3 U3() const {return omega_;}

    inline u3::SU3 SU3() const {return omega_.SU3();}

    inline HalfInt S() const {return S_;}

    ////////////////////////////////////////////////////////////////
    // key tuple, comparisons, and hashing
    ////////////////////////////////////////////////////////////////

    typedef std::tuple<u3::U3,HalfInt> KeyType;

    inline KeyType Key() const {return KeyType(omega_,S_);}

    inline friend bool operator == (const U3S& omegaS1, const U3S& omegaS2)
      {return omegaS1.Key() == omegaS2.Key();}

    inline friend bool operator < (const U3S& omegaS1, const U3S& omegaS2)
      {return omegaS1.Key() < omegaS2.Key();}

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

}  // namespace


// Defining std hash functions for SU3 and U3 class. 
namespace std
{

  template<> struct hash<u3::U3S>
  {
    inline std::size_t operator()(const u3::U3S& h) const noexcept
    {
      return hash_value(h);
    }
  };

}


// Defining formatters for fmt format for U3 and SU3 classes. 
namespace fmt
{
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
      return format_to(ctx.out(), "{:d}{:g}", wS.U3(),wS.S());
    }
  };
}  // namespace fmt

#endif
