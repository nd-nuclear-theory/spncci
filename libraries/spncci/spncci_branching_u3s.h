/****************************************************************
  spncci_branching_u3s.h

  U(3)xS layer of SpNCCI basis branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/31/17 (mac): Extracted from sp_basis.
    
****************************************************************/

#ifndef SPNCCI_BRANCHING_U3S_H_
#define SPNCCI_BRANCHING_U3S_H_

#include <unordered_map>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"  
#include "u3shell/upcoupling.h"
#include "lgi/lgi.h"

namespace spncci
{


  ////////////////////////////////////////////////////////////////
  // basis indexing
  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  // baby SpNCCI subspaces -- the most finely chopped SpNCCI
  // subspaces
  //
  // These are the U(3) subspaces distinguished by SpNCCI irrep, but
  // unioned over SpNCCI irreps in the same family, i.e., sharing the
  // same LGI labels.
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

  // labels

  typedef std::tuple<u3::U3,HalfInt,HalfInt,HalfInt,u3::U3> BabySpNCCISubspaceLabels;

  // subspace

  class BabySpNCCISubspace
    : public basis::BaseSubspace<spncci::BabySpNCCISubspaceLabels,int>
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
        const sp3r::U3Subspace& u3_subspace
      );
    // Construct from native description of SpNCCI space.
    //
    // Arguments:
    //   spncci_irrep_family (spncci::SpNCCIIrrepFamily): irrep family from which to take subspace
    //   u3_subspace (sp3r::U3Subspace): U(3) subspace

    // accessors

    u3::U3 sigma() const {return std::get<0>(labels_);}
    HalfInt Sp() const {return std::get<1>(labels_);}
    HalfInt Sn() const {return std::get<2>(labels_);}
    HalfInt S() const {return std::get<3>(labels_);}
    u3::U3 omega() const {return std::get<4>(labels_);}

    // diagnostic strings

    std::string LabelStr() const;
    // Provide string representation of subspace labels.

    std::string DebugStr() const;
    // Dump subspace contents (or, rather, dimension info).

    private:
    
    // dimension info
    int gamma_max_, upsilon_max_;

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

    BabySpNCCISpace(const spncci::SpNCCISpace& spncci_space);
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

    BabySpNCCISectors() {};
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

    SubspaceU3S (const u3::U3S& omegaS,const SpNCCISpace& spncci_space);

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

    SpaceU3S(SpNCCISpace& spncci_space);

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

  // typedef std::unordered_map<spncci::SectorLabelsU3S,int,boost::hash<spncci::SectorLabelsU3S>> SectorLabelsU3SCache;

  void GetSectorsU3S(
    const spncci::SpaceU3S& space, 
    const std::vector<u3shell::IndexedOperatorLabelsU3S>& relative_tensor_labels,
    std::vector<spncci::SectorLabelsU3S>& sector_vector
    );
  // Generates a cache of SectorLabelsU3S from operator labels given in 
  // relative_tensor_rmes, which are U(1)xSU(3)xSU(2) tensors labeled
  // by (N0,x0,S0,kappa0,L0). 
  // 
  // space (input) : space used to define sectors
  // relative_tensor_rmes (input) : container of rme labels keys and rme values
  //                                RelativeRMEsU3ST defined in upcoupling.h
  // u3_sectors (output) : container with SectorLabelsU3S keys and index values


}  // namespace

#endif