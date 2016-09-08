/****************************************************************
  sp_basis.h

  Sp(3,R) basis construction, indexing, and branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.
  3/11/16 (aem,mac): Implement basis irrep construction and traversal.
  3/17/16 (aem,mac): Remove superfluous data members from LGI.
  9/8/16  (aem): Renamed LGI to SpIrrep and moved LGI read function to lgi.h

****************************************************************/

#ifndef SP_BASIS_H_
#define SP_BASIS_H_


#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "u3shell/u3spn_scheme.h"  
#include "lgi/lgi.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) SpIrrep enumeration
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // SpIrrep data structure
  ////////////////////////////////////////////////////////////////

  // Stored for each SpIrrep:
  //   U3SPN labels -- inherited, so as to inherit accesors, etc.
  //   Nex -- useful when using SpIrrep in recurrence relation formulas
  //   pointer to irrep indexing
  //
  // Note: It might be useful to store the running index (0-based) in
  // the SpIrrep as well, so that it can appear in, e.g., the debugging
  // string output.
  //
  // Note: Sp and Sn are spectators for most conceivable calculations
  // but are retained for informational value.
  //
  // Note: Although sigma x S could be encoded together as a single
  // u3::U3S object, it is not clear that there is any benefit
  // (conceptual or practical) to adding an extra layer of packaging
  // at this stage in the code.

  // class SpIrrep  
  class SpIrrep
    : public lgi::LGI
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
      SpIrrep(const lgi::LGI& lgi)
      : LGI(lgi), sp_space_ptr_(NULL)
    {
    }

    ////////////////////////////////////////////////////////////////
    // initialization
    ////////////////////////////////////////////////////////////////

    void SaveSubspaceInfo(const sp3r::Sp3RSpace& sp_space)
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

    // Note: inherits all U3SPN accessors

    // int Nex() const {return Nex_;}
    u3::U3 sigma() const {return U3();}

    const sp3r::Sp3RSpace& Sp3RSpace() const
    {
      return *sp_space_ptr_;
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;
    std::string DebugString() const;

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

    // quick-reference information
    const sp3r::Sp3RSpace* sp_space_ptr_;
  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) SpIrrep -> U(3) branching
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // truncation scheme definition
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
  // storage of Sp(3,R) -> U(3) branchings
  ////////////////////////////////////////////////////////////////

  // Sp(3,R) container convenience type
  typedef std::map<u3::U3,sp3r::Sp3RSpace> SigmaIrrepMap;
  typedef std::vector<spncci::SpIrrep> SpIrrepVector;

  void GenerateSp3RIrreps(
      const lgi::LGIVector& lgi_vector,
      const TruncatorInterface& truncator,
      SpIrrepVector& sp_irrep_vector,
      SigmaIrrepMap& sigma_irrep_map
    );
  // Generate Sp(3,R) irrep branching information required for given set of LGIs.
  //
  // Persistent storage of the Sp3RSpace structures is in
  // sigma_irrep_map, but all the relevant information for access is
  // stored directly with the SpIrrep.
  //
  // Arguments:
  //   lgi_vector (LGIVector) : vector of LGIs for which we
  //     must generate branchings (INPUT)
  //    sp_irrep_vector(SpIrrepVector) : vector of Sp(3,R)xSU(2) irreps with pointer to 
  //      SigmaIrrepMap which contains the branched reduced states (OUTPUT)
  //   sigma_irrep_map (SigmaIrrepMap) : 
  //     mappings from sigma to its irrep indexing (OUTPUT)
  //   truncator (TruncatorInterface) : "truncator" object, which provides a 
  //     function to map each given sigma to its desired Nn_max truncation

  ////////////////////////////////////////////////////////////////
  // iteration over Sp(3,R) -> U(3) subspaces & states
  ////////////////////////////////////////////////////////////////

  // Traversal schemes
  //
  // To traverse Sp(3,R) -> U(3) subspaces:
  //
  //   for each SpIrrep (in SpIrrep_vector)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //
  // To traverse Sp(3,R) -> U(3) states:
  //
  //   for each SpIrrep (in sp_irrep_vector)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //       for each state index
  //
  // To traverse Sp(3,R) -> U(3) -> (L,S) states:
  //
  //   A triangularity constraint may be placed on L to restrict only
  //   to those L contributing to a desired final J.
  //
  //   for each SpIrrep (in sp_irrep_vector)
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
  //   for each SpIrrep (in sp_irrep_vector)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //       for each state index
  //         for each L
  //           for each J





  // total dimension counting functions

  int TotalU3Subspaces(const SpIrrepVector& sp_irrep_vector);
  // Compute total number of U(3) subspaces from given SpIrrep set.
  //
  // That is, we are counting the number of Sp(3,R) -> U(3) subspaces
  // in the branching, as defined under "Traversal schemes" above, not
  // just the number of distinct final U(3) labels.

  int TotalDimensionU3(const SpIrrepVector& sp_irrep_vector);
  // Compute total number of U(3)xS reduced states from the given SpIrrep
  // set.

  int TotalDimensionU3LS(const SpIrrepVector& sp_irrep_vector);
  // Compute total number of LxS branched states from the given SpIrrep
  // set.

  int TotalDimensionU3LSJConstrained(const SpIrrepVector& sp_irrep_vector, const HalfInt& J);
  // Compute total number of LxS->J branched states of a given J, from
  // the given SpIrrep set.
  //
  // From the point of view of counting and enumerating, this is more
  // simply the number of LxS branched states contributing to the
  // given given J, from the given SpIrrep set.

  int TotalDimensionU3LSJAll(const SpIrrepVector& sp_irrep_vector);
  // Compute total number of LxS->J branched states, over all J, from
  // the given SpIrrep set.
  //
  // This is equivalent to the M-space dimension for the universal
  // donor value of M (M=0 or 1/2).


  std::vector< std::pair<int,int> > GenerateSpIrrepPairs(spncci::SpIrrepVector sp_irrep_vector);
  // Given a vector of SpIrrep's, apply angular momentum selection rules to return a list 
  // of SpIrrep pairs which will have non-zero matrix elements between states in their irreps.  

  // Selection rules: abs(Si-Sf)<=2 for total spin, neutron spin and proton spin.
  


}  // namespace

#endif
