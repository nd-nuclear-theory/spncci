/****************************************************************
  spncci_space.h

  SpNCCI basis storage at Sp(3,R) irrep level.

  Includes branched dimension calculation functions.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  03/10/16 (aem,mac): Created (sp_basis).
  03/11/16 (aem,mac): Implement basis irrep construction and traversal.
  03/17/16 (aem,mac): Remove superfluous data members from LGI.
  09/08/16  (aem): Rename LGI to SpIrrep and moved LGI read function to lgi.h
  12/05/16 (aem): Change SpIrrepVector from vector of SpIrreps to vector of
                  MultiplicityTagged SpIrreps where tag is gamma_max.
  01/26/17 (mac): Change SigmaIrrepMap from map to unordered_map.
  01/31/17 (mac):
    - Rename to spncci_basis.
    - Rename and restructure SpNCCI irrep family containers.
    - Split off branched basis definitions.
  02/17/17 (mac):
    - Extract BabySpNCCI indexing from spncci_branching_u3s.
    - Add U3SPN accessors sigmaSPN and omegaSPN to BabySpNCCISubspace.
  02/19/17 (mac): Move in PrecomputeKMatrices from explicit.cpp.
  02/21/17 (mac):
    - Add intrinsic coordinate mode for PrecomputeKMatrices.
    - Impose explicit attribute on SpNCCI space constructor.
  06/07/17 (mac):
    - Extract PrecomputeKMatrices to vcs_cache.
    - Extract legacy GenerateSpNCCIIrrepFamilyPairs to unit_tensor_test.
  07/01/17 (aem) : Added intrinsic option for Nsigma0ForNuclide.
  09/27/17 (aem) : Removed gamma_max=0 lgi from spncci space
  10/04/17 (aem): Modified Sp3r->U(3) branching restriction
  10/11/17 (aem) : Moved Nsigma0ForNuclide to lgi.h
  01/16/18 (aem) : Added new BabySpNCCIHypersector contructor for updated recurrence scheme
  01/31/18 (aem) : Add ObservableBabySpNCCIHypersector class
  06/21/19 (aem) : Add BabySpNCCIHypersectors constructor for seed hypersectors
  02/12/21 (aem) : Add SpaceSpBasis and sectors extracted from branching2.{cpp,h}
****************************************************************/

#ifndef SPNCCI_BASIS_H_
#define SPNCCI_BASIS_H_

#include <array>
#include <unordered_map>

#include "basis/hypersector.h"
#include "basis/degenerate.h"

#include "lgi/lgi.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_common.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"
#include "u3shell/unit_tensor_space_u3s.h"
#include "u3shell/upcoupling.h"

namespace spncci
{

  typedef std::pair<int,int> NnPair;
  typedef std::pair<int,int> LGIPair;

  int ValenceShellForNuclide(const lgi::NuclideType& nuclide);
  // Calculate valence shell for nuclide.
  //
  // This is the highest oscillator shell occupied in the lowest
  // Pauli-allowed configuration (i.e., N1v).
  //
  // Example:
  //
  //   spncci::ValenceShellForNuclide({3,3});
  //
  //      => returns 1
  //
  // Note that this is shorthand for
  //
  //   spncci::ValenceShell(spncci::NuclideType({3,3}))
  //
  // made possible by automatic conversion from initializer list to
  // array.
  //
  //
  // Arguments:
  //   nuclide (input): (N,Z) for nucleus
  //
  // Returns:
  //   N1v

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
  // Alternative truncator
  ////////////////////////////////////////////////////////////////
  class NlimitTruncator
    : public TruncatorInterface
  {
    public:
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////

    // construct by given Nmax
    inline NlimitTruncator(const HalfInt& Nsigma_0, int Nmax, int Nlimit)
      : Nsigma_0_(Nsigma_0), Nmax_(Nmax), Nlimit_(Nlimit){}

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
    int Nlimit_;
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
      : LGI(lgi), gamma_max_(gamma_max), sp_space_ptr_(nullptr)
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
      SigmaIrrepMap& sigma_irrep_map,
      bool restrict_sp3r_to_u3_branching=false
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
  // linearized indexing for SpNCCI space
  //
  // baby SpNCCI subspaces -- the most finely chopped SpNCCI
  // subspaces
  //
  // These are the U(3) subspaces (or, equivalently, U3SPN subspaces)
  // distinguished by SpNCCI irrep, but unioned over SpNCCI irreps in
  // the same family, i.e., sharing the same LGI labels.
  //
  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (sigma,Sp,Sn,S,omega)
  //
  //   sigma (U3): LGI U(3) label
  //   Sp, Sn, S (HalfInt): spin labels
  //   omega (U3): subspace U(3) label
  //
  // state labels within subspace: (gamma,upsilon) -> (dummy)
  //
  //   gamma (int): index of SpNCCI irrep within family defined by LGI
  //     quantum numbers (1<=gamma<=gamma_max)
  //
  //   upsilon (int): branching multiplicity label for sigma->omega
  //     (1<=upsilon<=upsilon_max); represents composite of (n,rho),
  //     where n is the raising U(3) and rho is the sigma x n coupling
  //     multiplicity index
  //
  //   However, we suppress storage of the state labels and instead
  //   keep track of subspace dimensions.  A dummy int label is
  //   instead provided for purposes of the template argument.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  //   No states are actually defined.  The subspace multiplicity
  //   information (gamma_max,upsilon_max) are instead stored,
  //   yielding total dimension gamma_max*upsilon_max.
  //
  //   The conventional ordering of these states would be
  //   lexicographic by (gamma,upsilon).
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Subspaces are ordered by traversal of the SpNCCISpace:
  //
  //   -- by irrep family ordering (determined by LGI input order)
  //
  //   -- by U(3) subspace within irrep family's Sp3RSpace
  //     (i.e., canonically by increasing omega)
  //
  // Dimensions are calculated as gamma_max*upsilon_max.
  //
  ////////////////////////////////////////////////////////////////

  class BabySpNCCISubspace;
  class BabySpNCCISpace;

  // labels

  typedef std::tuple<u3::U3,HalfInt,HalfInt,HalfInt,u3::U3> BabySpNCCISubspaceLabels;

  // subspace

  class BabySpNCCISubspace
    : public basis::BaseDegenerateSubspace<BabySpNCCISubspace,spncci::BabySpNCCISubspaceLabels,basis::BaseState<BabySpNCCISubspace>,int>
  // SubspaceLabelsType (u3shell::U3SPN): (omega,S)
  // StateLabelsType (int): 1 (just a place holder)
  {
    public:

    // constructor

    BabySpNCCISubspace() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes

    BabySpNCCISubspace(
        const spncci::SpNCCIIrrepFamily& spncci_irrep_family,
        int irrep_family_index,
        const sp3r::U3Subspace& u3_subspace
      );
    // Construct from native description of SpNCCI space.
    //
    // Arguments:
    //   spncci_irrep_family (spncci::SpNCCIIrrepFamily): irrep family from which to take subspace
    //   u3_subspace (sp3r::U3Subspace): U(3) subspace

    // accessors

    u3::U3 sigma() const {return std::get<0>(labels());}
    HalfInt Sp() const {return std::get<1>(labels());}
    HalfInt Sn() const {return std::get<2>(labels());}
    HalfInt S() const {return std::get<3>(labels());}
    u3::U3 omega() const {return std::get<4>(labels());}
    int Nn() const {return int(omega().N()-sigma().N());}

    u3::U3S sigmaS() const {return u3::U3S(sigma(),S());}
    u3shell::U3SPN sigmaSPN() const {return u3shell::U3SPN(sigma(),Sp(),Sn(),S());}
    u3::U3S omegaS() const {return u3::U3S(omega(),S());}
    u3shell::U3SPN omegaSPN() const {return u3shell::U3SPN(omega(),Sp(),Sn(),S());}

    int irrep_family_index() const {return irrep_family_index_;}
    int upsilon_max() const {return upsilon_max_;}
    int gamma_max() const {return gamma_max_;}
    // diagnostic strings

    std::string LabelStr() const;
    // Provide string representation of subspace labels.

    std::string DebugStr() const;
    // Dump subspace contents (or, rather, dimension info).

    private:

    // dimension info
    int gamma_max_, upsilon_max_;
    // backwards lookup into SpNCCISpace
    int irrep_family_index_;
  };

  // space
  class BabySpNCCISpace
    : public basis::BaseSpace<BabySpNCCISubspace>
  {

  public:

    // constructor

    BabySpNCCISpace() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes

    explicit BabySpNCCISpace(const spncci::SpNCCISpace& spncci_space);
    // Construct from native description of SpNCCI space.
    //
    // Arguments:
    //   spncci_space (spncci::SpNCCISpace): space from which to harvest subspaces

    // diagnostic strings

    // std::string DebugStr() const;

  private:

  };

  // sectors

  class BabySpNCCISectors
    : public basis::BaseSectors<BabySpNCCISpace>
  {

  public:

    // constructor

    BabySpNCCISectors() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes

    BabySpNCCISectors(
        const spncci::BabySpNCCISpace& space,
        const u3shell::OperatorLabelsU3S& operator_labels
      );
      // Enumerate sector pairs connected by an operator of given
      // tensorial character ("constrained" sector
      // enumeration).
      //
      // Arguments:
      //   space (BabySpNCCISpace): the space
      //   operator_labels (OperatorLabelsU3S): NxSU(3)xS character of operator
  };

  class BabySpNCCIHypersectors
    : public basis::BaseHypersectors<BabySpNCCISpace,u3shell::RelativeUnitTensorSpaceU3S>
  {

  public:
      // constructor

    BabySpNCCIHypersectors() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes

    BabySpNCCIHypersectors(
        const spncci::BabySpNCCISpace& space,
        const u3shell::RelativeUnitTensorSpaceU3S& operator_space
      );
      // Enumerate sector pairs connected by u3S subspaces of
      // relative unit tensors
      //
      // Arguments:
      //   space (BabySpNCCISpace): the space
      //   operator_space (u3shell::RelativeUnitTensorSpaceU3S) : operator space


    BabySpNCCIHypersectors(
        const spncci::BabySpNCCISpace& space,
        const u3shell::ObservableSpaceU3S& operator_space,
        int irrep_family_index_1=-1, int irrep_family_index_2=-1
      );
      //Overload of notation
      // Enumerate sector pairs connected by u3S subspaces of
      // relative observables for Nnp>=Nn
      //
      // Arguments:
      //   space (BabySpNCCISpace): the space
      //   operator_space (u3shell::ObservableSpaceU3S) : observable operator space
      //   irrep_family_index=-1 means no restriction on which irrep family.

      BabySpNCCIHypersectors(
        int Nmax,
        const spncci::BabySpNCCISpace& space,
        const u3shell::RelativeUnitTensorSpaceU3S& operator_space,
        const std::map<spncci::NnPair,std::set<int>>& operator_subsets_NnpNn,
        std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
        int irrep_family_index_1=-1, int irrep_family_index_2=-1,
        bool Nn0_conjugate_hypersectors=false
      );

      //Overload of notation
      // Enumerate unit tensor hypersector connected by baby spncci subspaces
      //
      // Arguments:
      //    space : full baby spncci space
      //    operator_space : unit tensor space.  Subspaces by x0,S0,etap,eta
      //
      //    operator_subset : vector of unit_tensor_subspace indices which have
      //      may have non-zero rmes between the given irrep family pair.
      //
      //    unit_tensor_hypersector_subset: is a vector of vectors of indices
      //      for unit hypersectors organized by Nsum=Nn+Nnp.  Each previous
      //      Nsum must be computed before the recurrence can go to the next Nsum.
      //
      //    irrep_family_index : index of irrep family in SpNCCISpace.  Used to
      //      restrict hypersectors to include baby spncci subspaces for a given
      //      irrep family pair.  If irrep_family_index=-1, then there is no
      //      restriction by irrep family.
      //    Nn0_conjugate_hypersectors : Generate hypersectors for special case that
      //      Nn=0 and Nnp!=0.  irrep_family_index_bra and irrep_family_index_ket are
      //      passed as flipped, so hypersectors enumerate for Nnp'=0 and Nn'!=0.

      BabySpNCCIHypersectors(
        const lgi::MultiplicityTaggedLGIVector& lgi_families,
        const spncci::BabySpNCCISpace& space,
        const u3shell::RelativeUnitTensorSpaceU3S& operator_space,
        const std::set<int>& operator_subset,
        // std::vector<int>& unit_tensor_hypersector_subset,
        int irrep_family_index_1, int irrep_family_index_2
      );
      // Constructor for setting up hypersectors for only the LGI of each irrep family
      // for use in storing seed RMEs for unit tensor recurrence.

  };

  void PrintHypersectors(
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
    );


  class ObservableBabySpNCCIHypersectors
    : public basis::BaseHypersectors<BabySpNCCISpace,u3shell::ObservableSpaceU3S>
  {

  public:
      // constructor

    ObservableBabySpNCCIHypersectors() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes

    ObservableBabySpNCCIHypersectors(
        const spncci::BabySpNCCISpace& baby_spncci_space,
        const u3shell::ObservableSpaceU3S& observable_space,
        int irrep_family_index_1=-1, int irrep_family_index_2=-1, bool restrict_lower_triangle=false
      );
      // Enumerate hypersectors between baby spncci subsaces for a set of operators in observables space
      //
      // Arguments:
      //    baby_spncci_space : Sp3rU3S subspaces
      //    observable_space : operator space [U(1)xSU(3)xSU(2) tensors]
      //    irrep_family_index_1, irrep_family_index_2 : index of irrep family of interest
      //      + if -1, then construct for all pairs of irrep families, not just one pair
      //    restrict_lower_triangle : if  true, only construct hyperspectors for lower triangle of
      //      full matrix.  (irrep_family_index_bra>=irrep_family_index_ket)

  };


  void NumBabySpncciSubspacesInIrrepFamily(
    int num_irrep_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    std::vector<int>& irrep_family_num_subspaces
  );

  void SortSubspacesDecending(const std::vector<int>& num_subspaces,std::vector<int>&ordered_subspaces);


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Branched symplectic basis
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // convenience typedef for (L,S)
  typedef std::pair<u3::U3,int> omegaLLabels;

  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to J level
  ////////////////////////////////////////////////////////////////
  //
  //   subspace: (sigma,Sp,Sn,S) irrep_family_index
  //     state: (omega,L)
  //       substates: (kappa,gamma,upsilon)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (sigma,Sp,Sn,S) => u3shell::U3SPN
  //
  // state labels within subspace: (omega,L) => omegaLLabels
  //
  // substate labels (implied): (kappa,gamma,upsilon)
  //
  //   (See BabySpNCCI docstring in spncci_basis for definitions of
  //   these basis labels.)
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, states are ordered by first appearance of
  // omega in irrep_family (sigma,Sp,Sn,S) then by L in order of
  // appearence in u3::BranchingSO3(omega).
  //
  // Degeneracy is kappa_max*gamma_max*upsilon_max.
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by appearance in
  // babyspncci
  //
  // ALTERNATE CODE AVAILABLE BUT NOT SELECTED:
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (omega,S).

  ////////////////////////////////////////////////////////////////

  class SubspaceSpBasis;
  class StateSpBasis;
  class SpaceSpBasis;

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class SubspaceSpBasis
    : public basis::BaseDegenerateSubspace<SubspaceSpBasis,int,StateSpBasis,omegaLLabels>
    // SubspaceLabelsType (int) : irrep_family_index
    //     Formerly   SubspaceLabelsType (u3::U3SPN): (sigma,Sp,Sn,S)
    // StateLabelsType (omegaLLabels): (omega,L)
    {
      public:

      // constructors

      SubspaceSpBasis() = default;
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      SubspaceSpBasis(
        const HalfInt& J,
        const u3shell::U3SPN& sigmaSPN,
        int irrep_family_index,
        const BabySpNCCISpace& baby_spncci_space
      );

      // subspace label accessors
      // u3shell::U3SPN sigmaSPN() const {return labels();}
      u3shell::U3SPN sigmaSPN() const {return sigmaSPN_;}

      HalfInt N() const {return sigmaSPN().U3().N();}
      u3::U3 sigma() const {return sigmaSPN().U3();}
      HalfInt S() const {return sigmaSPN().S();}

      // int irrep_family_index() const {return irrep_family_index_;}
      int irrep_family_index() const {return labels();}
      int gamma_max() const {return gamma_max_;}

      // state auxiliary data accessors
      const std::vector<int>& state_kappa_max() const {return state_kappa_max_;}
      const std::vector<int>& state_upsilon_max() const {return state_upsilon_max_;}
      const std::vector<int>& state_baby_spncci_subspace_index() const {return state_baby_spncci_subspace_index_;}

      // diagnostic output
      std::string LabelStr() const;
      std::string DebugStr() const;

      private:

      int gamma_max_;
      int irrep_family_index_;
      u3shell::U3SPN sigmaSPN_;

      // state auxiliary data
      std::vector<int> state_kappa_max_;
      std::vector<int> state_upsilon_max_;
      std::vector<int> state_baby_spncci_subspace_index_;
    };

  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateSpBasis
    : public basis::BaseDegenerateState<SubspaceSpBasis>
  {

    public:

    // pass-through constructors

    StateSpBasis(const SubspaceType& subspace, int& index)
      // Construct state by index.
      : basis::BaseDegenerateState<SubspaceSpBasis>(subspace,index) {}

    StateSpBasis(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
      // Construct state by reverse lookup on labels.
      : basis::BaseDegenerateState<SubspaceSpBasis>(subspace,state_labels)
      {}

    // pass-through accessors for subspace labels
    u3shell::U3SPN sigmaSPN() const {return subspace().sigmaSPN();}
    u3::U3 sigma() const {return subspace().sigma();}
    HalfInt S() const {return subspace().S();}
    HalfInt N() const {return subspace().N();}


    // state label accessors
    omegaLLabels omegaL() const {return labels();}
    u3::U3 omega() const {return labels().first;}
    int L() const {return labels().second;}
    int Nn() const
    {
      return int(omega().N()-N());
    }

    // diagnostic output
    // std::string LabelStr() const;

    // state auxiliary data accessors
    int upsilon_max() const
    {
      return subspace().state_upsilon_max()[index()];
    }
    int kappa_max() const
    {
      return subspace().state_kappa_max()[index()];
    }

    int baby_spncci_subspace_index() const
    {
      return subspace().state_baby_spncci_subspace_index()[index()];
    }

    private:

  };

  // ////////////////////////////////////////////////////////////////
  // // space
  // ////////////////////////////////////////////////////////////////

  class SpaceSpBasis
    : public basis::BaseSpace<SubspaceSpBasis>
  {

    public:

    // constructor
    SpaceSpBasis() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J);

    //Alternate constructor that constructs a subspace of the full space base on irrep_family_subset
    SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J, std::set<int>irrep_family_subset);

    HalfInt J() const {return J_;}

    // diagnostic output
    std::string DebugStr(bool show_subspaces=false) const;

    private:
      HalfInt J_;

  };


  ////////////////////////////////////////////////////////////////
  // sectors
  ////////////////////////////////////////////////////////////////

  class SectorsSpBasis
    : public basis::BaseSectors<SpaceSpBasis>
  {

    public:

    // constructor

    SectorsSpBasis() = default;
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SectorsSpBasis(
        const SpaceSpBasis& space,
        HalfInt J0,
        basis::SectorDirection sector_direction = basis::SectorDirection::kCanonical
      );
    // Enumerate sector pairs connected by an operator of given
    // tensorial character.

  };


  void GetSpBasisOffsets(
    const spncci::SpaceSpBasis& spbasis,
    std::vector<std::vector<int>>& offsets
    );
    // Iterate through J-branched spbasis and identify starting position
    // of each irrep family in the basis listing.  Starting index stored
    // in offsets indexed by subspace index. Subspaces are (sigma,Sp,Sn,S)





}  // namespace

#endif
