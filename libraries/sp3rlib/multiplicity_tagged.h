/****************************************************************
  multiplicity_tagged.h

  Container classes for irrep labels with multiplicity.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/9/16 (aem,mac): Extracted from u3.h.
  6/25/17 (mac): Restore and fix template specialization for
    MultiplicityTagged<int>::Str().

****************************************************************/

#ifndef MULTIPLICITY_TAGGED_H_
#define MULTIPLICITY_TAGGED_H_

#include <sstream>
#include <string>
// #include <utility>
#include <vector>
 
#include <boost/functional/hash.hpp>

template <typename tIrrep>
struct MultiplicityTagged
// Container classes for irrep labels with multiplicity.
//
// EX:
//   MultiplicityTagged<u3::SU3> xrho;
//   xrho.irrep = u3::SU3(2,1);
//   xrho.multiplicity = 4;
{

  ////////////////////////////////////////////////////////////////
  // convenience typedefs
  ////////////////////////////////////////////////////////////////

  // convenience typedef for container class
  typedef std::vector<MultiplicityTagged<tIrrep> > vector;
      
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////
      
  // copy constructor: synthesized copy constructor since only data
  // member needs copying

  // default constructor
  inline MultiplicityTagged() 
    : tag(0) {}

  // construct by (irrep, tag)
  inline MultiplicityTagged(const tIrrep& irrep_, int tag_) 
    : irrep(irrep_), tag(tag_) {}

  ////////////////////////////////////////////////////////////////
  // key tuple, comparisons, and hashing
  ////////////////////////////////////////////////////////////////

  typedef std::pair<tIrrep,int> KeyType;

  inline KeyType Key() const
  {
    return KeyType(irrep,tag);
  }

  inline friend bool operator == (const MultiplicityTagged<tIrrep>& x1, const MultiplicityTagged<tIrrep>& x2)
  {
    return x1.Key() == x2.Key();
  }
  
  inline friend bool operator < (const MultiplicityTagged<tIrrep>& x1, const MultiplicityTagged<tIrrep>& x2)
  {
    return x1.Key() < x2.Key();
  }

  inline friend std::size_t hash_value(const MultiplicityTagged<tIrrep>& v)
  {
    boost::hash<MultiplicityTagged<tIrrep>::KeyType> hasher;
    return hasher(v.Key());
  }

  ////////////////////////////////////////////////////////////////
  // string conversion
  ////////////////////////////////////////////////////////////////
    
  std::string Str() const;

  ////////////////////////////////////////////////////////////////
  // labels
  ////////////////////////////////////////////////////////////////
      
  tIrrep irrep;
  int tag;
};

template <typename tIrrep>
std::string MultiplicityTagged<tIrrep>::Str() const
// Generate string output relying on Str() method of irrep.
//
// Note: Will fail if irrep type does not have Str() method, e.g., if
// the irrep is just and int.  This problem may be overcome via
// template specialization (see <int> specialization below).
{
  std::ostringstream ss;
	
  ss << "(" << irrep.Str() << "," << tag << ")";
  return ss.str();
}

template <>
inline
std::string MultiplicityTagged<int>::Str() const
// Template specialization for tIrrep->int.
//
// Programming note: Beware that a fully specialized template no
// longer is treated as a template under the "one definition rule" and
// must therefore be explicitly inlined to avoid link-time errors
// "multiple definition of `MultiplicityTagged<int>::Str() const".
// 
//   https://stackoverflow.com/questions/4445654/multiple-definition-of-template-specialization-when-using-different-objects
{
  std::ostringstream ss;
	
  ss << "(" << irrep << "," << tag << ")";
  return ss.str();
}


#endif
