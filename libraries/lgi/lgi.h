/****************************************************************
  lgi.h

  Interface for lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/8/16 (aem,mac): Created.
  1/31/17 (mac): Rename LGIVector to MultiplicityTaggedLGIVector.
  2/17/17 (mac): Extract WriteLGILabels from lgi_solver and change
    to accept std::ostream for output and extract to lgi.
  10/11/17 (aem) : Extract Nsigma0ForNuclide from spncci_basis
  1/15/18 (aem) : Updated Read LGI to be consistant with write functions
****************************************************************/
#ifndef LGI_SOLVER_H_
#define LGI_SOLVER_H_

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "u3shell/u3spn_scheme.h"  
#include "lsu3shell/lsu3shell_rme.h"
namespace lgi
{

  ////////////////////////////////////////////////////////////////
  // Calculation of Nsigma0
  ////////////////////////////////////////////////////////////////

  typedef std::array<int,2> NuclideType;

  HalfInt Nsigma0ForNuclide(const NuclideType& nuclide, bool intrinsic=false);
  // Calculate Nsigma0 for nuclide.
  //
  // This may be thought of as the dimensionless "oscillator energy in
  // the lowest Pauli-allowed configuration", including zero-point
  // energy.
  //
  // Example:
  //
  //   spncci::Nsigma0ForNuclide({3,3});
  //
  //      => returns 11
  //
  // Note that this is shorthand for
  //
  //   spncci::Nsigma0ForNuclide(spncci::NuclideType({3,3}))
  //
  // made possible by automatic conversion from initializer list to
  // array.
  //
  // Arguments:
  //   nuclide (input): (N,Z) for nucleus
  //
  // Returns:
  //   Nsigma0


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

  void
    WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_families,std::ostream& os);

  void
    WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_families, const std::string& filename);


  void 
  WriteLGIExpansion(
    int Z, int N, int Nmax,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    lsu3shell::OperatorBlocks&lgi_expansions,
    const std::string& filename
  );


  void 
  ReadLGISet(
    const std::string& lgi_filename, 
    const HalfInt& Nsigma0,
    MultiplicityTaggedLGIVector& lgi_vector
  );
  // Read in LGI from file and create vector of LGIs tagged by 
  // gamma_max from tabulation in file .
  //
  // Arguments:
  //   filename (string) : filename for LGI table file
  //   Nsigma0 : minimum number of oscillator quanta for the given system of nucleons.
  //    Can be obtained from lgi::Nsigma0ForNuclide. 
  //   lgi_families (MultiplicityTaggedLGIVector) : container for LGI list (OUTPUT)  

 void ReadLGILookUpTable(std::vector<int>& lgi_full_space_lookup_table, int num_irrep_families);
  // Reading in and filling out table of lgi indices in basis and lgi indices in full space by 
  // which the seed files are labeled. 


}

#endif
