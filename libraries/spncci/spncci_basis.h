/****************************************************************
  spncci_basis.h

  SpNCCI basis storage at Sp(3,R) irrep level.

  Includes branched dimension calculation functions.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created (sp_basis).
  3/11/16 (aem,mac): Implement basis irrep construction and traversal.
  3/17/16 (aem,mac): Remove superfluous data members from LGI.
  9/8/16  (aem): Rename LGI to SpIrrep and moved LGI read function to lgi.h
  12/5/16 (aem): Change SpIrrepVector from vector of SpIrreps to 
                 vector of MultiplicityTagged SpIrreps where tag is 
                 number of times sigma Sp Sn S occur in basis.
  1/26/17 (mac): Change SigmaIrrepMap from map to unordered_map.
  1/31/17 (mac):
    - Rename to spncci_basis.
    - Rename and restructure SpNCCI irrep family containers.
    - Split off branched basis definitions.
    
****************************************************************/

#ifndef SPNCCI_BASIS_H_
#define SPNCCI_BASIS_H_

#include <unordered_map>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"  
#include "u3shell/upcoupling.h"
#include "lgi/lgi.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // Sp irrep truncation scheme definition
  ////////////////////////////////////////////////////////////////

  // The Nn_max for each sigma might be determined by any number of
  // creative truncation schemes, and we put in place a general mechanism for this through the TruncatorInterface.
  //
  // The truncation is most simply defined based on an Nmax cutoff:
  //
  //   Nn_max = Nmax - Nex
  //          = Nmax - (Nsigma - Nsigma_0)

  class TruncatorInterface
  // Define generic interface for defining truncations by (sigma,S).
  {
    public:
    virtual int Nn_max(const u3::U3& sigma) const
      = 0; //pure virtual
    // Calculate Nn_max to use for given SpIrrep (sigma,S) labels.
    //
    // Arguments:
    //   sigmaS (u3::U3S) : (sigma,S) labels of the SpIrrep
    //
    // Returns:
    //   (int) : the truncation Nn_max to use for this SpIrrep
  };

  class NmaxTruncator 
    : public TruncatorInterface
  {
    public:
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////
    
    // construct by given Nmax
    inline NmaxTruncator(const HalfInt& Nsigma_0, int Nmax) 
      : Nsigma_0_(Nsigma_0), Nmax_(Nmax) {}

    ////////////////////////////////////////////////////////////////
    // truncator calculation
    ////////////////////////////////////////////////////////////////

    virtual int Nn_max(const u3::U3& sigma) const;

    ////////////////////////////////////////////////////////////////
    // truncation information
    ////////////////////////////////////////////////////////////////
    
    private:
    HalfInt Nsigma_0_;
    int Nmax_;

  };

  ////////////////////////////////////////////////////////////////
  // SpNCCI space as set of irrep families
  ////////////////////////////////////////////////////////////////

  class SpNCCIIrrepFamily
    : public lgi::LGI
  // Storage of family of Sp irreps (with common LGI labels).
  //
  // Data members:
  //   (lgi::LGI): inherited Nex-tagged U3SPN labels for irrep family
  //   gamma_max (int): multiplicity of these LGI labels
  //
  // Note: Inheritance from lgi::LGI is evil use of inheritance, since
  // it fails the "is a type of" test.  I.e., a gamma-labeled family
  // Sp(3,R) irreps is not a "type of" lowest grade SU(3) irrep.
  {
    public:
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying
    // null constructor
    
    // inline SpIrrep():Nex_(-999){};
    // inline SpIrrep():{};
    inline
      SpNCCIIrrepFamily(const lgi::LGI& lgi, int gamma_max)
      : LGI(lgi), gamma_max_(gamma_max), sp_space_ptr_(NULL)
    {
    }

    ////////////////////////////////////////////////////////////////
    // initialization
    ////////////////////////////////////////////////////////////////

    void SaveSpaceReference(const sp3r::Sp3RSpace& sp_space)
    // Save information on Sp3RSpace associated with this SpIrrep's sigma
    // for quick reference, i.e., without requiring a map lookup.
    //
    // Caution: Since only a *reference* to the Sp3RSpace is
    // mainatined here, it is important that the Sp3RSpace still be
    // safely stored elsewhere, e.g., in a sigma_irrep_map, without
    // going "out of scope" and being destroyed.
    {
      sp_space_ptr_ = &sp_space;
    }

    //////////////////////////////////////////////////////////////
    //accessors
    //////////////////////////////////////////////////////////////

    // Note: inherits all LGI (and thus U3SPN) accessors

    // alias for U3 content
    u3::U3 sigma() const {return U3();}

    // family size
    //
    // i.e., multiplicity of LGIs with givensigma x Sp x Sn x S
    int gamma_max() const {return gamma_max_;};
    
    // irrep branching info
    const sp3r::Sp3RSpace& Sp3RSpace() const
    {
      return *sp_space_ptr_;
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;
    std::string DebugStr() const;

    //////////////////////////////////////////////////////////////
    //key tuple, comparisons and hashing
    //////////////////////////////////////////////////////////////

    // basic ordering and hashing functions inherited from U3SPN

    // pseudo-key for std::tie access to SpIrrep properties
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

    // Note: LGI labels are inherited.

    // irrep family size
    int gamma_max_;
    
    // quick-reference information
    const sp3r::Sp3RSpace* sp_space_ptr_;
  };

  // SpNCCISpace
  //
  // full set of Sp irrep families
  typedef std::vector<spncci::SpNCCIIrrepFamily> SpNCCISpace;

  // Sp(3,R) branching container
  //
  // Necessary to hold Sp->U(3) branching info containers
  // (sp3r::Sp3RSpace) so that they are not destructed.
  typedef std::unordered_map<u3::U3,sp3r::Sp3RSpace,boost::hash<u3::U3>> SigmaIrrepMap;
  // typedef std::map<u3::U3,sp3r::Sp3RSpace> SigmaIrrepMap;

  void GenerateSpNCCISpace(
      const lgi::MultiplicityTaggedLGIVector& lgi_vector,
      const TruncatorInterface& truncator,
      SpNCCISpace& spncci_space,
      SigmaIrrepMap& sigma_irrep_map
    );
  // Generate Sp(3,R) irrep branching information required for given
  // set of LGIs.
  //
  // Persistent storage of the Sp3RSpace structures is in
  // sigma_irrep_map, but all the relevant information for access is
  // stored directly with the SpIrrep.
  //
  // Arguments:
  //   lgi_vector (MultiplicityTaggedLGIVector): vector of LGIs for which we
  //     must generate branchings (INPUT)
  //   truncator (TruncatorInterface): "truncator" object, which provides a 
  //     function to map each given sigma to its desired Nn_max truncation
  //   spncci_space(SpNCCISpace): vector of Sp(3,R)xSU(2) irreps with pointer to 
  //     SigmaIrrepMap which contains the branched reduced states (OUTPUT)
  //   sigma_irrep_map (SigmaIrrepMap): 
  //     mappings from sigma to its irrep indexing (OUTPUT)

  ////////////////////////////////////////////////////////////////
  // dimension counting by basis traveral
  ////////////////////////////////////////////////////////////////

  // Traversal schemes for iteration over Sp(3,R) -> U(3) subspaces &
  // states
  //
  // To traverse Sp(3,R) -> U(3) subspaces:
  //
  //   for each SpIrrep (in SpIrrep_vector)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //
  // To traverse Sp(3,R) -> U(3) states:
  //
  //   for each SpIrrep (in spncci_space)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //       for each state index
  //
  // To traverse Sp(3,R) -> U(3) -> (L,S) states:
  //
  //   A triangularity constraint may be placed on L to restrict only
  //   to those L contributing to a desired final J.
  //
  //   for each SpIrrep (in spncci_space)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //       for each state index
  //         for each L  (optionally subject to triangularity constraint for J)
  //
  // To traverse Sp(3,R) -> U(3) -> (L,S) -> J states:
  //
  //   A triangularity constraint may be placed on L to restrict only
  //   to those L contributing to a desired final J.
  //
  //   for each SpIrrep (in spncci_space)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //       for each state index
  //         for each L
  //           for each J

  // total dimension counting functions

  int TotalU3Subspaces(const SpNCCISpace& spncci_space);
  // Compute total number of U(3) subspaces from given SpIrrep set:
  //
  //    gamma sigma Sp Sn S -> omega
  //
  // That is, we are counting the number of Sp(3,R) -> U(3) subspaces
  // in the branching, as defined under "Traversal schemes" above, not
  // just the number of distinct final U(3) labels.
  //
  // It might make more sense just to count the number of spnnci
  // "leaves", blind to gamma.

  int TotalDimensionU3S(const SpNCCISpace& spncci_space);
  // Compute total number of U(3)xS reduced states from the given SpIrrep
  // set.

  int TotalDimensionU3LS(const SpNCCISpace& spncci_space);
  // Compute total number of LxS branched states from the given SpIrrep
  // set.

  int TotalDimensionU3LSJConstrained(const SpNCCISpace& spncci_space, const HalfInt& J);
  // Compute total number of LxS->J branched states of a given J, from
  // the given SpIrrep set.
  //
  // From the point of view of counting and enumerating, this is more
  // simply the number of LxS branched states contributing to the
  // given given J, from the given SpIrrep set.

  int TotalDimensionU3LSJAll(const SpNCCISpace& spncci_space);
  // Compute total number of LxS->J branched states, over all J, from
  // the given SpIrrep set.
  //
  // This is equivalent to the M-space dimension for the universal
  // donor value of M (M=0 or 1/2).

  ////////////////////////////////////////////////////////////////
  // SpNCCI irrep family pair enumeration
  ////////////////////////////////////////////////////////////////

  std::vector< std::pair<int,int> >
    GenerateSpNCCIIrrepFamilyPairs(spncci::SpNCCISpace spncci_space);
  // Enumerate pairs of Sp NCCI irrep families connected under
  // two-body allowed spin selection rules.
  //
  // Given a vector of SpIrrep's, apply angular momentum selection
  // rules to return a list of SpIrrep pairs which will have non-zero
  // matrix elements between states in their irreps.  Mapping is all
  // to all.
  //
  // Selection rules are based on fact that a two-body operator (or
  // relative operator) can carry at most 2 units of spin, by species
  // or in total.
  //
  // Selection rules: abs(Si-Sf)<=2 for total spin, neutron spin and proton spin.
  //
  // DEPRECATED (but still used in some test code)

}  // namespace

#endif