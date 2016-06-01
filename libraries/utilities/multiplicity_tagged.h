/****************************************************************
  multiplicity_tagged.h

  Container classes for irrep labels with multiplicity.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/9/16 (aem,mac): Extracted from u3.h.

****************************************************************/

#ifndef MULTIPLICITY_TAGGED_H_
#define MULTIPLICITY_TAGGED_H_

#include <sstream>
#include <string>
#include <utility>
#include <vector>
 

template <typename IRREP>
struct MultiplicityTagged
// Container classes for irrep labels with multiplicity.
//
// EX:
//   MultiplicityTagged<u3::SU3> xrho;
//   xrho.irrep = u3::SU3(2,1);
//   xrho.multiplicity = 4;
{

  ////////////////////////////////////////////////////////////////
  // typedefs
  ////////////////////////////////////////////////////////////////

  // typedef for sort key
  typedef std::pair<IRREP,int> KeyType;

  // convenience typedef for container class
  typedef std::vector<MultiplicityTagged<IRREP> > vector;
      
  ////////////////////////////////////////////////////////////////
  // constructors
  ////////////////////////////////////////////////////////////////
      
  // copy constructor: synthesized copy constructor since only data
  // member needs copying

  // default constructor
  inline MultiplicityTagged() 
    : tag(0) {}

  // construct by (irrep, tag)
  inline MultiplicityTagged(const IRREP& irrep_, int tag_) 
    : irrep(irrep_), tag(tag_) {}

  ////////////////////////////////////////////////////////////////
  // accessors
  ////////////////////////////////////////////////////////////////

  inline KeyType Key() const
  {
    return KeyType(irrep,tag);
  }

  ////////////////////////////////////////////////////////////////
  // string conversion
  ////////////////////////////////////////////////////////////////
    
  std::string Str() const;

  ////////////////////////////////////////////////////////////////
  // labels
  ////////////////////////////////////////////////////////////////
      
  IRREP irrep;
  int tag;
};

template <typename IRREP>
std::string MultiplicityTagged<IRREP>::Str() const
// Generate string output relying on Str() method of irrep.
//
// Note: Will fail if irrep type does not have Str() method, e.g.,
// if the irrep is just and int.  This problem may, in principle,
// be overcome via template specialization, but this leads to
// link-time errors (gcc 4.5).
{
  std::ostringstream ss;
	
  ss << "(" << irrep.Str() << "," << tag << ")";
  return ss.str();
}

// template <>
//   std::string MultiplicityTagged<int>::Str() const
//   // Template specialization for IRREP->int.
//   {
//     std::ostringstream ss;
// 	
//     ss << "(" << irrep << "," << tag << ")";
//     return ss.str();
//   }

////////////////////////////////////////////////////////////////
// relational operators
////////////////////////////////////////////////////////////////

template <typename IRREP>
inline bool operator == (const MultiplicityTagged<IRREP>& x1, const MultiplicityTagged<IRREP>& x2)
{
  return x1.Key() == x2.Key();
}

template <typename IRREP>
inline bool operator < (const MultiplicityTagged<IRREP>& x1, const MultiplicityTagged<IRREP>& x2)
{
  return x1.Key() < x2.Key();
}


#endif
