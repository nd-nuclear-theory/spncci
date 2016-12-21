/****************************************************************
  lgi.h

  Interface for lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/8/16 (aem,mac): Created.
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
  // LGI data structure
  ////////////////////////////////////////////////////////////////

  // Stored for each LGI:
  //   U3SPN labels -- inherited, so as to inherit accesors, etc.
  //   Nex -- useful when using LGI in recurrence relation formulas
  //   pointer to irrep indexing
  //
  // Note: It might be useful to store the running index (0-based) in
  // the LGI as well, so that it can appear in, e.g., the debugging
  // string output.
  //
  // Note: Sp and Sn are spectators for most conceivable calculations
  // but are retained for informational value.
  //
  // Note: Although sigma x S could be encoded together as a single
  // u3::U3S object, it is not clear that there is any benefit
  // (conceptual or practical) to adding an extra layer of packaging
  // at this stage in the code.

  // class LGI  
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
  typedef MultiplicityTagged<lgi::LGI>::vector LGIVector;

  void ReadLGISet(LGIVector& lgi_vector, const std::string& lgi_filename);

  // Generates vector of LGIs based on LGI input tabulation.
  //
  // Arguments:
  //   lgi_vector (LGIVector) : container for LGI list (OUTPUT)
  //   filename (string) : filename for LGI table file
}

#endif