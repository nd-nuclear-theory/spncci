/****************************************************************
  lgi.h

  Interface for lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/8/16 (aem,mac): Created.
  1/31/17 (mac): Rename LGIVector to MultiplicityTaggedLGIVector.
****************************************************************/
#ifndef LGI_SOLVER_H_
#define LGI_SOLVER_H_

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "u3shell/u3spn_scheme.h"  

namespace lgi
{
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) LGI enumeration
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // LGI labels
  ////////////////////////////////////////////////////////////////

  // class LGI
  //
  // U3SPN labels "dressed" with easy access to the Nex quantum number
  //
  // U3SPN (parent): omega x Sp x Sn x S labels
  // Nex (int): N relative to Nsigma0 for nucleus

  class LGI
    : public u3shell::U3SPN
  {
    public:
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying
    // null constructor
    
    inline LGI():Nex_(-999){};
    inline
      LGI(const u3shell::U3SPN& sigmaSPN, int Nex)
      : U3SPN(sigmaSPN), Nex_(Nex)
    {}
    //////////////////////////////////////////////////////////////
    //accessors
    //////////////////////////////////////////////////////////////

    // Note: inherits all U3SPN accessors

    int Nex() const {return Nex_;}
    u3::U3 sigma() const {return U3();}
    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;
    //////////////////////////////////////////////////////////////
    //key tuple, comparisons and hashing
    //////////////////////////////////////////////////////////////

    // basic ordering and hashing functions inherited from U3SPN

    // pseudo-key for std::tie access to LGI properties
    //
    // This Key is *not* the key actually used in sorting and hashing,
    // which is the parent U3SPN type's key.

    typedef std::tuple<int,u3::U3,HalfInt,HalfInt,HalfInt> KeyType;

    inline KeyType Key() const
    {
      return KeyType(Nex(),sigma(),Sp(),Sn(),S());
    }


    ////////////////////////////////////////////////////////////////
    // data
    ////////////////////////////////////////////////////////////////
    private:

    // Note: U3SPN labels are inherited.

    // quick-reference information
    int Nex_;
  };

  ////////////////////////////////////////////////////////////////
  // enumeration of LGI set based on input table
  ////////////////////////////////////////////////////////////////

  // LGI input tabulation format:
  //
  //   Nex 2N lambda mu 2Sp 2Sn 2S count
  //   ...
  //
  // Each input table line results in multiple stored LGIs based on
  // the given count.
  //
  // LGI container convenience type
  //
  // STYLE: maybe LGI::vector would be more consistent
  typedef MultiplicityTagged<lgi::LGI>::vector MultiplicityTaggedLGIVector;

  void ReadLGISet(MultiplicityTaggedLGIVector& lgi_vector, const std::string& lgi_filename);

  // Generates vector of LGIs based on LGI input tabulation.
  //
  // Arguments:
  //   lgi_vector (MultiplicityTaggedLGIVector) : container for LGI list (OUTPUT)
  //   filename (string) : filename for LGI table file
}

#endif
