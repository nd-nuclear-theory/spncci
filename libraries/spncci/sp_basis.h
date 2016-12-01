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
  // Only compute upper triangle sectors, i.e., gamma'<=gamma
  //
  // Selection rules: abs(Si-Sf)<=2 for total spin, neutron spin and proton spin.
  

  ////////////////////////////////////////////////////////////////
  // basis indexing in U3S scheme for spncci basis branching
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (omega,S) = U3S
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // States are indexed by a tuple of numbers [gamma,sigmaS,dim,index]
  // where gamma which corresponds to the index of the lgi in lgi_vector.
  // dim=nu_max is symplectic multiplicity of omega in the lgi
  // index is the position of the nu=0 state for the given lgi in the 
  // subspace list of states labeld by gamma(sigma,Sp,Sn), nu. 
  //
  // However, we keep track of the subspace dimension.  The subspace
  // dimensions must be provided to the constructor via a mapping U3S
  // -> dimension.  Then we must explicitly set the dimension_ member
  // variable (inherited from BaseSubspace).
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (omega,S).
  // 
  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class SubspaceU3S
    : public basis::BaseSubspace<u3::U3S,std::tuple<int,u3::U3,int,int>>
  // Subspace class for two-body states of given U(3)xS.
  //
  // SubspaceLabelsType (u3shell::U3S): (omega,S)
  // StateLabelsType (std::tuple<int,int,int>): (gamma, dim, index)
  {
    public:

    // constructor
    SubspaceU3S() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SubspaceU3S (const u3::U3S& omegaS,const SpIrrepVector& sp_irrep_vector);

    // accessors
    u3::U3S omegaS() const {return labels_;}
    u3::U3 omega() const {return omegaS().U3();}
    u3::SU3 x() const {return omegaS().SU3();}
    HalfInt N() const {return omegaS().U3().N();}
    HalfInt S() const {return omegaS().S();}
    int sector_dim() const {return sector_size_;}
    // diagnostic output
    std::string Str() const;

    private:
      int sector_size_;
  };
  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateU3S
    : public basis::BaseState<SubspaceU3S>
  // State class for two-body states of given U(3)xSxT.
  {
    
  public:
    // pass-through constructors
  
  StateU3S(const SubspaceType& subspace, int index)
    // Construct state by index.
    : basis::BaseState<SubspaceU3S>(subspace, index) {}

  StateU3S(
    const SubspaceType& subspace,
    const typename SubspaceType::StateLabelsType& state_labels
    )
    // Construct state by reverse lookup on labels.
    : basis::BaseState<SubspaceU3S> (subspace, state_labels) 
    {}

    // pass-through accessors
    u3::U3S omegaS() const {return Subspace().omegaS();}
    u3::U3 omega() const {return Subspace().omega();}
    HalfInt S() const {return Subspace().S();}
    HalfInt N() const {return Subspace().N();}

    // state label accessors
    int gamma() const {return std::get<0>(GetStateLabels());} 
    u3::U3 sigma() const {return std::get<1>(GetStateLabels());}

    int sector_dim() const {return std::get<2>(GetStateLabels());}
    int index() const {return std::get<3>(GetStateLabels());}

    std::string Str() const;
  private:
 
  };




  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class SpaceU3S
    : public basis::BaseSpace<SubspaceU3S>
  // Space class for two-body states of given U(3)xS.
  {
    
  public:

    // constructor
    SpaceU3S() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceU3S(SpIrrepVector& sp_irrep_vector);

    // diagnostic output
    std::string Str() const;

  private:
    int dimension_;
  };





  ////////////////////////////////////////////////////////////////
  // Sector
  // Enumerates upper triangle sectors omegaS'<=omegaS
  ////////////////////////////////////////////////////////////////
  class SectorLabelsU3S
  {
  public:
// Need N0,x0,S0,kappa0,L0, rho0
    typedef std::tuple<int,int,u3shell::OperatorLabelsU3S,int,int,int> KeyType;
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////
    //default constructor
    inline SectorLabelsU3S()
    :rho0_(0), kappa0_(0), L0_(0){}

    // construction from labels
    inline 
    SectorLabelsU3S(
      int bra_index, int ket_index, 
      const u3shell::OperatorLabelsU3S& tensor_labels,
      int kappa0, int L0, int rho0
      )
    : bra_index_( bra_index), ket_index_(ket_index), tensor_labels_(tensor_labels),kappa0_(kappa0), L0_(L0), rho0_(rho0)
    {}

    inline
    SectorLabelsU3S(
      int bra_index, int ket_index, 
      const u3shell::IndexedOperatorLabelsU3S& tensor_labels,
      int rho0
      )
    : bra_index_( bra_index), ket_index_(ket_index), rho0_(rho0)
    {
      std::tie(tensor_labels_,kappa0_,L0_)=tensor_labels; 
    }

    inline 
    SectorLabelsU3S(
      int bra_index, int ket_index, 
      const u3shell::OperatorLabelsU3ST& tensor_labels,
      int kappa0, int L0, int rho0
      )
    : bra_index_( bra_index), ket_index_(ket_index),kappa0_(kappa0), L0_(L0), rho0_(rho0)
    {
      tensor_labels_=u3shell::OperatorLabelsU3S(tensor_labels);
    }

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
    int bra_index() const {return bra_index_;}
    int ket_index() const {return ket_index_;}
    u3shell::OperatorLabelsU3S 
      operator_labels() const {return tensor_labels_;}
    int N0() const {return tensor_labels_.N0();}
    u3::SU3 x0() const {return tensor_labels_.x0();}
    HalfInt S0() const {return tensor_labels_.S0();}
    int L0() const {return L0_;}
    int kappa0() const {return kappa0_;}
    int rho0() const {return rho0_;}

    inline KeyType Key() const
    {
      return KeyType(bra_index_,ket_index_,tensor_labels_,kappa0_,L0_,rho0_);
    }

    inline friend bool operator == (const SectorLabelsU3S& sector1, const SectorLabelsU3S& sector2)
    {
      return sector1.Key() == sector2.Key();
    }

    inline friend bool operator < (const SectorLabelsU3S& sector1, const SectorLabelsU3S& sector2)
    {
      return sector1.Key() < sector2.Key();
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(const SectorLabelsU3S& sector)
    {
      boost::hash<SectorLabelsU3S::KeyType> hasher;
      return hasher(sector.Key());
    }
    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////

    std::string Str() const;

  private:
    int bra_index_, ket_index_, kappa0_, L0_,rho0_;
    u3shell::OperatorLabelsU3S tensor_labels_;
  };


  typedef std::unordered_map<spncci::SectorLabelsU3S,int,boost::hash<spncci::SectorLabelsU3S>> SectorLabelsU3SCache;

  void GetSectorsU3S(
    const spncci::SpaceU3S& space, 
    const std::vector<u3shell::IndexedOperatorLabelsU3S>& relative_tensor_labels,
    SectorLabelsU3SCache& u3s_sectors
    );
  // Generates a cache of SectorLabelsU3S from operator labels given in 
  // relative_tensor_rmes, which are U(1)xSU(3)xSU(2) tensors labeled
  // by (N0,x0,S0,kappa0,L0). 
  // 
  // space (input) : space used to define sectors
  // relative_tensor_rmes (input) : container of rme labels keys and rme values
  //                                RelativeRMEsU3ST defined in upcoupling.h
  // u3_sectors (output) : container with SectorLabelsU3S keys and index values

  ////////////////////////////////////////////////////////////////
  // basis indexing in LS scheme for spncci basis branching
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (L,S) = LS
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // States are indexed by a tuple of numbers [omegaS, index]
  // omegaS correspond to subspace in SpaceU3S
  // index is starting position in sector matrix 
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (L,S).
  // 
  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class SubspaceLS
    : public basis::BaseSubspace<std::tuple<int,HalfInt>,std::tuple<u3::U3,int,int>>
  // Subspace class for two-body states of given SO(3)xS.
  //
  {
    public:

    // constructor

    SubspaceLS() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SubspaceLS(const int& L, const HalfInt& S,const SpaceU3S& u3s_space);

    // accessors
    HalfInt S() const {return std::get<1>(GetSubspaceLabels());}
    int L() const{return std::get<0>(GetSubspaceLabels());}
    int sector_dim() const {return sector_size_;}
    // diagnostic output
    std::string Str() const;

    private:
      int sector_size_;
  };
  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateLS
    : public basis::BaseState<SubspaceLS>
  // State class for two-body states of given U(3)xSxT.
  {
    
  public:
    // pass-through constructors
  
  StateLS(const SubspaceType& subspace, int index)
    // Construct state by index.
    : basis::BaseState<SubspaceLS>(subspace, index) {}

  StateLS(
    const SubspaceType& subspace,
    const typename SubspaceType::StateLabelsType& state_labels
    )
    // Construct state by reverse lookup on labels.
    : basis::BaseState<SubspaceLS> (subspace, state_labels) 
    {}

    // pass-through accessors
    HalfInt S() const {return Subspace().S();}
    int L() const {return Subspace().L();}
    // state label accessors
    u3::U3 omega() const {return std::get<0>(GetStateLabels());}
    int kappa_max() const {return std::get<1>(GetStateLabels());}
    int index() const {return std::get<2>(GetStateLabels());}
    std::string Str() const;
  private:
 
  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class SpaceLS
    : public basis::BaseSpace<SubspaceLS>
  // Space class for two-body states of given U(3)xS.
  {
    
  public:

    // constructor

    SpaceLS() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceLS(const SpaceU3S& u3s_space, HalfInt J);

    // diagnostic output
    std::string Str() const;

  private:
    int dimension_;
  };



}  // namespace

#endif
