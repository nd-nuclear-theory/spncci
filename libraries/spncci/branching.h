/****************************************************************
  branching.h

  Basis definitions for U(3)xS, LxS, and J branchings of SpNCCI basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/6/17 (mac): Created, as reimplementation of branching_u3s and
    branching_u3lsj.
****************************************************************/

#ifndef SPNCCI_BRANCHING_H_
#define SPNCCI_BRANCHING_H_

#include "basis/degenerate.h"
// #include "am/am.h"
#include "sp3rlib/u3.h"
#include "spncci/spncci_basis.h"
// #include "spncci/unit_tensor.h"
// #include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"  
// #include "u3shell/upcoupling.h"
// #include "lgi/lgi.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // convenience typedef for (L,S)
  ////////////////////////////////////////////////////////////////

  typedef std::tuple<int,HalfInt> LSLabels;

  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to U3S level
  ////////////////////////////////////////////////////////////////  
  //
  //   subspace: (omega,S)
  //     state: (sigma,Sp,Sn)
  //       substates: (gamma,upsilon)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (omega,S) => u3::U3S
  //
  // state labels within subspace: (sigma,Sp,Sn,[S]) => u3shell::U3SPN
  //   i.e., irrep family labels (of course, S is redundant)
  //
  // substate labels (implied): (gamma,upsilon)
  //
  //   (See BabySpNCCI docstring in spncci_basis for definitions of
  //   these basis labels.)
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, states are ordered by first appearance of
  // irrep family label in traveral of the BabySpNCCI basis, which in
  // turn follows irrep family ordering in SpNCCISpace.
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by first appearance
  // of (omega,S) in traveral of the BabySpNCCI basis.
  //
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class SubspaceSpU3S
    : public basis::BaseDegenerateSubspace<u3::U3S,u3shell::U3SPN>
    // SubspaceLabelsType (u3::U3S): (omega,S)
    // StateLabelsType (u3shell::U3SPN): (sigma,Sp,Sn,S)
    {
      public:

      // constructors

      SubspaceSpU3S() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      SubspaceSpU3S(const u3::U3S& omegaS, const BabySpNCCISpace& baby_spncci_space);

      // subspace label accessors
      u3::U3S omegaS() const {return labels_;}
      u3::U3 omega() const {return omegaS().U3();}
      u3::SU3 x() const {return omegaS().SU3();}
      HalfInt N() const {return omegaS().U3().N();}
      HalfInt S() const {return omegaS().S();}

      // state auxiliary data accessors
      const std::vector<int>& state_gamma_max() const {return state_gamma_max_;}
      const std::vector<int>& state_baby_spncci_subspace_index() const {return state_baby_spncci_subspace_index_;}
      const std::vector<int>& state_irrep_family_index() const {return state_irrep_family_index_;}

      // diagnostic output
      std::string LabelStr() const;
      std::string DebugStr() const;

      private:

      // state auxiliary data
      std::vector<int> state_gamma_max_;
      std::vector<int> state_baby_spncci_subspace_index_;
      std::vector<int> state_irrep_family_index_;
    };

  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateSpU3S
    : public basis::BaseDegenerateState<SubspaceSpU3S>
  {
    
    public:

    // pass-through constructors

    StateSpU3S(const SubspaceType& subspace, int& index)
      // Construct state by index.
      : basis::BaseDegenerateState<SubspaceSpU3S>(subspace,index) {}

    StateSpU3S(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
      // Construct state by reverse lookup on labels.
      : basis::BaseDegenerateState<SubspaceSpU3S>(subspace,state_labels) 
      {}

    // pass-through accessors
    u3::U3S omegaS() const {return subspace().omegaS();}
    u3::U3 omega() const {return subspace().omega();}
    HalfInt S() const {return subspace().S();}
    HalfInt N() const {return subspace().N();}


    // state label accessors
    u3shell::U3SPN sigmaSPN() const {return labels();}
    int Nn() const 
    {
      return int(subspace().N()-sigmaSPN().U3().N());
    }

    // state auxiliary data accessors
    int gamma_max() const
    {
      return subspace().state_gamma_max()[index()];
    }
    int baby_spncci_subspace_index() const
    {
      return subspace().state_baby_spncci_subspace_index()[index()];
    }

    int irrep_family_index() const
    {
      return subspace().state_irrep_family_index()[index()];
    }


    private:
 
  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class SpaceSpU3S
    : public basis::BaseDegenerateSpace<SubspaceSpU3S>
  {
    
    public:

    // constructor
    SpaceSpU3S() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceSpU3S(const BabySpNCCISpace& baby_spncci_space);

    // diagnostic output
    std::string DebugStr(bool show_subspaces=false) const;

  };



  ////////////////////////////////////////////////////////////////
  // Sector
  // Enumerates omegaS sectors
  ////////////////////////////////////////////////////////////////

  class SectorLabelsSpU3S
  {
    public:
    // Need N0,x0,S0,kappa0,L0, rho0
    typedef std::tuple<int,int,u3shell::OperatorLabelsU3S,int,int,int> KeyType;
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////
    //default constructor
    inline SectorLabelsSpU3S()
      :rho0_(0), kappa0_(0), L0_(0){}

    // construction from labels
    inline 
      SectorLabelsSpU3S(
          int bra_index, int ket_index, 
          const u3shell::OperatorLabelsU3S& tensor_labels,
          int kappa0, int L0, int rho0
        )
      : bra_index_( bra_index), ket_index_(ket_index), tensor_labels_(tensor_labels),kappa0_(kappa0), L0_(L0), rho0_(rho0)
    {}

    inline
      SectorLabelsSpU3S(
          int bra_index, int ket_index, 
          const u3shell::IndexedOperatorLabelsU3S& tensor_labels,
          int rho0
        )
      : bra_index_( bra_index), ket_index_(ket_index), rho0_(rho0)
    {
      std::tie(tensor_labels_,kappa0_,L0_)=tensor_labels; 
    }

    inline 
      SectorLabelsSpU3S(
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

    inline friend bool operator == (const SectorLabelsSpU3S& sector1, const SectorLabelsSpU3S& sector2)
    {
      return sector1.Key() == sector2.Key();
    }

    inline friend bool operator < (const SectorLabelsSpU3S& sector1, const SectorLabelsSpU3S& sector2)
    {
      return sector1.Key() < sector2.Key();
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(const SectorLabelsSpU3S& sector)
    {
      boost::hash<SectorLabelsSpU3S::KeyType> hasher;
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

  void GetSectorsSpU3S(
      const spncci::SpaceSpU3S& space, 
      const std::vector<u3shell::IndexedOperatorLabelsU3S>& relative_tensor_labels,
      std::vector<spncci::SectorLabelsSpU3S>& sector_vector
    );





  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to LS level
  ////////////////////////////////////////////////////////////////  
  //
  //   subspace: (L,S)
  //     state: (omega,kappa,sigma,Sp,Sn)
  //       substates: (gamma,upsilon)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (L,S) => spncci::LSLabels
  //
  // state labels within subspace: (omega,kappa,(sigma,Sp,Sn,[S]))
  //   => (u3::U3,int,u3shell::U3SPN)
  //
  // substate labels (implied): (gamma,upsilon)
  //
  //   (See BabySpNCCI docstring in spncci_basis for definitions of
  //   these basis labels.)
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, states are ordered by:
  //
  //   - ordering of (omega,S) subspaces in SpaceSpU3S, and thus by
  //     first appearance of (omega,S) label in traveral of the
  //     BabySpNCCI basis
  //
  //   - increasing kappa
  //
  //   - ordering of irrep family labels within (omega,S) subspace,
  //     which follows first appearance of irrep family label in
  //     traveral of the BabySpNCCI basis, which in turn follows irrep
  //     family ordering in SpNCCISpace
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (L,S).  Only (L,S) values appearing by branching of an (omega,S)
  // subspace are included.
  //
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  typedef std::tuple<u3::U3,int,u3shell::U3SPN> StateLabelsSpLS;
  class SubspaceSpLS
    : public basis::BaseDegenerateSubspace<spncci::LSLabels,spncci::StateLabelsSpLS>
    {
      public:

      // constructors
      SubspaceSpLS() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      SubspaceSpLS(const spncci::LSLabels& ls_labels, const SpaceSpU3S& spu3s_space);

      // subspace label accessors
      int L() const {return std::get<0>(labels_);}
      HalfInt S() const {return std::get<1>(labels_);}

      // state auxiliary data accessors
      const std::vector<int>& state_gamma_max() const {return state_gamma_max_;}
      const std::vector<int>& state_baby_spncci_subspace_index() const {return state_baby_spncci_subspace_index_;}
      const std::vector<int>& state_spu3s_subspace_index() const {return state_spu3s_subspace_index_;}

      // diagnostic output
      std::string LabelStr() const;
      std::string DebugStr() const;

      private:

      // state auxiliary data
      std::vector<int> state_gamma_max_;
      std::vector<int> state_baby_spncci_subspace_index_;
      std::vector<int> state_spu3s_subspace_index_;
    };


  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateSpLS
    : public basis::BaseDegenerateState<SubspaceSpLS>
  {
    
    public:

    // pass-through constructors

    StateSpLS(const SubspaceType& subspace, int& index)
      // Construct state by index.
      : basis::BaseDegenerateState<SubspaceSpLS>(subspace,index) {}

    StateSpLS(
        const SubspaceType& subspace,
        const typename SubspaceType::StateLabelsType& state_labels
      )
      // Construct state by reverse lookup on labels.
      : basis::BaseDegenerateState<SubspaceSpLS>(subspace,state_labels) 
      {}

    // pass-through accessors
    int L() const {return subspace().L();}
    HalfInt S() const {return subspace().S();}

    // state label accessors -- (Sp,Sn)
    HalfInt Sp() const {return sigmaSPN().Sp();}
    HalfInt Sn() const {return sigmaSPN().Sn();}

    // state label accessors -- omega
    u3::U3 omega() const {return std::get<0>(labels());}
    u3::U3S omegaS() const {return u3::U3S(omega(),S());}
    u3shell::U3SPN omegaSPN() const {return u3shell::U3SPN(omega(),Sp(),Sn(),S());}
    HalfInt N() const {return omega().N();}
    // int Nex() const {return int(N()-Nsigma_0);}

    // state label accessors -- kappa
    int kappa() const {return std::get<1>(labels());}

    // state label accessors -- sigma
    u3::U3 sigma() const {return sigmaSPN().U3();}
    u3::U3S sigmaS() const {return sigmaSPN().U3S();}
    u3shell::U3SPN sigmaSPN() const {return std::get<2>(labels());}
    HalfInt Nsigma() const {return sigmaSPN().N();}
    // int Nsigmaex() const {return int(Nsigma()-Nsigma_0);}

    // state auxiliary data accessors
    int gamma_max() const
    {
      return subspace().state_gamma_max()[index()];
    }

    int baby_spncci_subspace_index() const
    // which baby_spncci_subspace it came from
    {
      return subspace().state_baby_spncci_subspace_index()[index()];
    }

    int spu3s_subspace_index() const
    // which spu3s_subspace it came from
    {
      return subspace().state_spu3s_subspace_index()[index()];
    }

    private:
 
  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class SpaceSpLS
    : public basis::BaseDegenerateSpace<SubspaceSpLS>
  {
    
    public:

    // constructor
    SpaceSpLS() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SpaceSpLS(const SpaceSpU3S& spu3s_space);
    // Construct full LS space.

    SpaceSpLS(const SpaceSpU3S& spu3s_space, HalfInt J);
    // Construct J-constrained LS space.


    // diagnostic output
    std::string DebugStr(bool show_subspaces=false) const;

  };


}  // namespace

#endif
