/****************************************************************
  u4.h

  Anna E. McCoy
  Argonne National Lab

  SPDX-License-Identifier: MIT

  + 01/15/25 (aem): Created.
****************************************************************/

#ifndef U4_H_
#define U4_H_

#include <string>
#include <tuple>
#include <array>
#include<math.h>
#include <iostream>
#include <cassert>
#include "fmt/format.h"
#include "sp3rlib/multiplicity_tagged.h"
#include "boost/functional/hash.hpp"

namespace u4
{
  
  class U4 
  {
    public:

    inline U4(): f1_(0),f2_(0),f3_(0),f4_(0){}
    inline U4(const unsigned int& f1, const unsigned int& f2, const unsigned int& f3,const unsigned int& f4)
        : f1_(f1), f2_(f2), f3_(f3),f4_(f4){}
    
    inline U4(const int& f1, const int& f2, const int& f3,const int& f4)
        : f1_(static_cast<unsigned int>(f1)), f2_(static_cast<unsigned int>(f2)), 
        f3_(static_cast<unsigned int>(f3)),f4_(static_cast<unsigned int>(f4)){
          assert(f1>0 && f2>0 && f3>0 && f4>0);
        }
    
    inline U4(const unsigned int& N, const std::tuple<unsigned int,unsigned int,unsigned int> su4_labels)
    {
      const auto& [t1,t2,t3] = su4_labels;
      int n((N-t1-t2-t3)/4); 
      assert(n>0);
      assert((N-t1-t2-t3)%4==0);

      f1_=t1+n;
      f2_=t2+n;
      f3_=t3+n;
      f4_=n;
    }

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // access Cartesian labels
    typedef std::tuple<unsigned int,unsigned int,unsigned int,unsigned int> KeyType;
    inline KeyType Key() const
    {
      return KeyType(f1_,f2_,f3_,f4_);
    }

    inline unsigned int one() const{return f1_;}
    inline unsigned int two() const{return f2_;}
    inline unsigned int three() const{return f3_;}
    inline unsigned int four() const{return f4_;}
    inline unsigned int N() const{return f1_+f2_+f3_+f4_;}
    inline KeyType labels() const {return Key();}

    /// Checks if given SU(4) irrep has Sn conjugate which is a valid U(N) irreps
    inline bool HasSnConjugate(unsigned int N) const {return one()<= N;}

    inline std::array<unsigned int, 4> as_array() const{
      std::array<unsigned int, 4> f={f1_,f2_,f3_,f4_};
      return f;
    }

    inline friend bool operator == (const U4& a, const U4& b)
    {
      return a.Key() == b.Key();
    }

    inline friend bool operator < (const U4& a, const U4& b)
    {
      return a.Key() < b.Key();
    }

    inline friend std::size_t hash_value(const U4& v)
    {
      boost::hash<U4::KeyType> hasher;
      return hasher(v.Key());
    }



    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    inline std::string Str() const
    {
      std::ostringstream ss;
      ss << "[" << f1_ << "," << f2_ << "," << f3_ <<"," << f4_ << "]";
      return ss.str();
    }

  private:
    unsigned int f1_,f2_,f3_,f4_;

  };


  /// Flip rows and columns of U(4) tableau f 
  /// In resulting conjugated tableau each row can have 0,..,4 tiles
  /// conjugate_tableau[0] is number of rows with 4, 
  /// conjugate_tableau[1] is number of rows with 3,
  /// etc.
  class U4Conjugate
  {
    public:

    inline U4Conjugate(): n4_(0),n3_(0),n2_(0),n1_(0),n0_(0){}
    inline U4Conjugate(
      const unsigned int num4, 
      const unsigned int num3,
      const unsigned int num2,
      const unsigned int num1,
      const unsigned int num0
    )
        : n0_(num0),n1_(num1), n2_(num2), n3_(num3),n4_(num4){}





    inline U4Conjugate(const u4::U4& f, const unsigned int N)     
    {
      assert(f.HasSnConjugate(N));
      const std::array<unsigned int,4>& f_tuple=f.as_array();
      std::vector<int> temp(f.one(),1);
      for(int i=1; i<4; i++)
        for (int j=0; j<f_tuple[i]; ++j)
          temp[j]++;

      n4_=std::count(temp.begin(), temp.end(), 4);
      n3_=std::count(temp.begin(), temp.end(), 3);
      n2_=std::count(temp.begin(), temp.end(), 2);
      n1_=std::count(temp.begin(), temp.end(), 1);
      n0_=N-f.one();
    }
    
    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////

    // access Cartesian labels
    typedef std::tuple<unsigned int,unsigned int,unsigned int,unsigned int,unsigned int> KeyType;
    inline KeyType Key() const
    {
      return KeyType(n4_,n3_,n2_,n1_,n0_);
    }

    inline unsigned int num0() const{return n0_;}
    inline unsigned int num1() const{return n1_;}
    inline unsigned int num2() const{return n2_;}
    inline unsigned int num3() const{return n3_;}
    inline unsigned int num4() const{return n4_;}
    

    inline friend bool operator == (const U4Conjugate& a, const U4Conjugate& b)
    {
      return a.Key() == b.Key();
    }

    inline friend bool operator < (const U4Conjugate& a, const U4Conjugate& b)
    {
      return a.Key() < b.Key();
    }

    inline friend std::size_t hash_value(const U4Conjugate& v)
    {
      boost::hash<U4Conjugate::KeyType> hasher;
      return hasher(v.Key());
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    inline std::string Str() const
    {
      std::ostringstream ss;
      ss << "[" << n4_ << "," << n3_ << "," << n2_ <<"," << n1_ << "," << n0_ << "]";
      return ss.str();
    }

  private:
    unsigned int n0_,n1_,n2_,n3_,n4_;

  };

  ///Based on equation 3.14 in jmp-39-1998-5631-Pan.
  unsigned int OuterMultiplicity(const u4::U4& lmbda, const u4::U4& mu, const u4::U4& nu);


  MultiplicityTagged<u4::U4>::vector KroneckerProduct(const u4::U4& f, const u4::U4& g);
  

  inline double SU4Casimir(const u4::U4& f)
  {
    unsigned int lambda1=f.one()-f.two();
    unsigned int lambda2=f.two()-f.three();
    unsigned int lambda3=f.three()-f.four();
    double cas = 3*lambda1
      +(3*lambda1*lambda1)/4
      +4*lambda2
      +lambda1*lambda2
      +lambda2*lambda2
      +3*lambda3
      +lambda1*lambda3/2
      +lambda2*lambda3+3*lambda3*lambda3/4;
    return cas;
  }


}


namespace std
{
  template<> struct hash<u4::U4>
  {
    inline std::size_t operator()(const u4::U4& h) const
    {
      return hash_value(h);
    }
  };
}


namespace fmt
{
  template<> struct formatter<u4::U4>
  {
    char presentation = 'g';

    template<typename ParseContext>
    FMT_CONSTEXPR auto parse(ParseContext& ctx) -> decltype(ctx.begin())
    {
      auto it = ctx.begin(), end = ctx.end();
      if (it != end && (*it == 'd' || *it == 'g' || *it == 'f'))
        presentation = *it++;

      // Check if reached the end of the range:
      if (it != end && *it != '}')
        throw format_error("invalid format");

      // Return an iterator past the end of the parsed range:
      return it;
    }

    template<typename FormatContext>
    FMT_CONSTEXPR auto format(const u4::U4& labels, FormatContext& ctx)
        -> decltype(ctx.out())
    {
      // if (presentation == 'f')
        return fmt::format_to(ctx.out(), "[{:d},{:d},{:d},{:d}]", labels.one(), labels.two(), labels.three(), labels.four());
      // else
      //   return fmt::format_to(ctx.out(), "{:g}{:d}", w.N(), w.SU3());
    }
  };
}  // namespace fmt

#endif