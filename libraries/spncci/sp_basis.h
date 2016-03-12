/****************************************************************
  sp_basis.h

  Sp(3,R) basis construction, indexing, and branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.
  3/11/16 (aem,mac): Implement basis irrep construction and traversal.

****************************************************************/

#ifndef SP_BASIS_H_
#define SP_BASIS_H_

#include <vector>

#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"

namespace spncci
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
  //   Nex sigma Sp Sn S
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
  
  struct LGI
  {
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // copy constructor: synthesized copy constructor since only data
    // member needs copying

    inline
    LGI(int Nex_, const u3::U3& sigma_, const HalfInt& Sp_, const HalfInt& Sn_, const HalfInt& S_)
      : Nex(Nex_), sigma(sigma_), Sp(Sp_), Sn(Sn_), S(S_), 
        Nn_max(0), dimension(0), irrep_ptr(NULL) {}

    ////////////////////////////////////////////////////////////////
    // initialization
    ////////////////////////////////////////////////////////////////

    void SaveSubspaceInfo(int Nn_max_, int dimension_, const sp3r::Sp3RSpace& irrep_)
    // Save information on Sp3RSpace associated with this LGI's sigma
    // for quick reference, i.e., without requiring a map lookup.
    //
    // Caution: Since only a *reference* to the Sp3RSpace is
    // mainatined here, it is important that the Sp3RSpace still be
    // safely stored elsewhere, e.g., in a sigma_irrep_map, without
    // going "out of scope" and being destroyed.
    {
      Nn_max = Nn_max_; dimension = dimension_; irrep_ptr = &irrep_;
    }


    ////////////////////////////////////////////////////////////////
    // retrieval
    ////////////////////////////////////////////////////////////////

    const sp3r::Sp3RSpace& Sp3RSpace() const
    {
      return *irrep_ptr;
    }

    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////
    
    std::string Str() const;
    std::string DebugString() const;

    ////////////////////////////////////////////////////////////////
    // data
    ////////////////////////////////////////////////////////////////
    
    // labels
    int Nex;
    u3::U3 sigma;
    HalfInt Sp, Sn, S;

    // quick-reference information
    int Nn_max;
    int dimension;
    const sp3r::Sp3RSpace* irrep_ptr;
  };

  ////////////////////////////////////////////////////////////////
  // enumeration of LGI set based on input table
  ////////////////////////////////////////////////////////////////

  // LGI input tabulation format:
  //
  //   Nex lambda mu 2Sp 2Sn 2S count
  //   ...
  //
  // Each input table line results in multiple stored LGIs based on
  // the given count.
  //
  // In calculating sigma, we also need Nsigma_0 for this nucleus, since
  //
  //   Nsigma = Nsigma_0 + Nex

  // LGI container convenience type
  //
  // STYLE: maybe LGI::vector would be more consistent
  typedef std::vector<spncci::LGI> LGIVectorType;

  void GenerateLGIVector(LGIVectorType& lgi_vector, const std::string& lgi_filename, const HalfInt& Nsigma_0);
  // Generates vector of LGIs based on LGI input tabulation.
  //
  // Arguments:
  //   lgi_vector (LGIVectorType) : container for LGI list (OUTPUT)
  //   filename (string) : filename for LGI table file
  //   Nsigma_0 (HalfInt) : Nsigma_0 for nucleus

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // Sp(3,R) LGI -> U(3) branching
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
    // Calculate Nn_max to use for given LGI (sigma,S) labels.
    //
    // Arguments:
    //   sigmaS (u3::U3S) : (sigma,S) labels of the LGI
    //
    // Returns:
    //   (int) : the truncation Nn_max to use for this LGI
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
  typedef std::map<u3::U3,sp3r::Sp3RSpace> SigmaIrrepMapType;

  void GenerateSp3RIrreps(
                          LGIVectorType& lgi_vector,
                          SigmaIrrepMapType& sigma_irrep_map,
                          const TruncatorInterface& truncator
                          );
  // Generate Sp(3,R) irrep branching information required for given set of LGIs.
  //
  // Persistent storage of the Sp3RSpace structures is in
  // sigma_irrep_map, but all the relevant information for access is
  // stored directly with the LGI.
  //
  // Arguments:
  //   lgi_vector (LGIVectorType) : vector of LGIs for which we
  //     must generate branchings (INPUT/OUTPUT)
  //   sigma_irrep_map (SigmaIrrepMapType) : 
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
  //   for each LGI (in lgi_vector)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //
  // To traverse Sp(3,R) -> U(3) states:
  //
  //   for each LGI (in lgi_vector)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //       for each state index
  //
  // To traverse Sp(3,R) -> U(3) -> (L,S) states:
  //
  //   A triangularity constraint may be placed on L to restrict only
  //   to those L contributing to a desired final J.
  //
  //   for each LGI (in lgi_vector)
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
  //   for each LGI (in lgi_vector)
  //     follow reference to Sp3RSpace
  //     for each U3Subspace in Sp3RSpace
  //       for each state index
  //         for each L
  //           for each J

  // total dimension counting functions

  int TotalU3Subspaces(const LGIVectorType& lgi_vector);
  // Compute total number of U(3) subspaces from given LGI set.
  //
  // That is, we are counting the number of Sp(3,R) -> U(3) subspaces
  // in the branching, as defined under "Traversal schemes" above, not
  // just the number of distinct final U(3) labels.

  int TotalDimensionU3(const LGIVectorType& lgi_vector);
  // Compute total number of U(3)xS reduced states from the given LGI
  // set.

  int TotalDimensionU3LS(const LGIVectorType& lgi_vector);
  // Compute total number of LxS branched states from the given LGI
  // set.

  int TotalDimensionU3LSJConstrained(const LGIVectorType& lgi_vector, const HalfInt& J);
  // Compute total number of LxS->J branched states of a given J, from
  // the given LGI set.
  //
  // From the point of view of counting and enumerating, this is more
  // simply the number of LxS branched states contributing to the
  // given given J, from the given LGI set.

  int TotalDimensionU3LSJAll(const LGIVectorType& lgi_vector);
  // Compute total number of LxS->J branched states, over all J, from
  // the given LGI set.

  

}  // namespace

#endif
