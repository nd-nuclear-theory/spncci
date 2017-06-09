/****************************************************************
  branching_u3lsj.h

  U(3)xLS and U(3)xLSJ layers of SpNCCI basis branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/31/17 (mac): Extracted from sp_basis (as spncci_branching_u3lsj).
  2/19/17 (mac): Rename to branching_u3lsj.
  5/27/17 (mac): Overhaul implementation of U3LS subspaces and store
    parent irrep info.
    
****************************************************************/

#ifndef SPNCCI_SPNCCI_BRANCHING_U3LSJ_H_
#define SPNCCI_SPNCCI_BRANCHING_U3LSJ_H_

#include <unordered_map>

#include "am/am.h"  
#include "lgi/lgi.h"
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "spncci/branching_u3s.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"  
#include "u3shell/upcoupling.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // basis indexing in U3LS scheme for spncci basis branching
  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (L,S) = LS
  //
  //   L (int): orbital angular momentum
  //   S (HalfInt): total spin
  //
  // state labels within subspace: (N)
  //
  //   subpace_index (int): index of U3S subspace from which state
  //     is branched
  //
  // ----------------
  //
  // The idea is that states are grouped in the hierarchy
  //
  //   subspace: (L,S)
  //     state: (omega) => actually u3s_subspace_index stored
  //       substates: (kappa,(sigma,Sp,Sn),gamma,upsilon)
  //
  // This unfortunately means we have "thrown away" genuine labels,
  // and our "states" no longer have good Sp symmetry, but rather are
  // aggregates of substates with different symmetry, i.e., from
  // different irrep families.  And it destroys all simple "Cartesian"
  // product structure to the substate indexing.  Although the kappa
  // multiplicity is defined by a single kappa_max (from omega->L),
  // now upsilon_max and gamma_max are no longer constant, but rather
  // depend on parent Sp irrep family (sigma,Sp,Sn,S).
  //
  // A more intuitive scheme, which retains this substate ordering,
  // but also retains all necessary information for decompositions by
  // symmetry labels, would be:
  //
  //   subspace: (L,S)
  //     state: (omega,kappa,sigma,Sp,Sn)
  //       substates: (gamma,upsilon)
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered lexicographically by
  // (L,S).
  //
  // NOT TRUE?
  //
  // Aren't they just ordered by first appearance of (L,S) as a
  // possible branching for an (omega,S) subspace, then by L?
  //
  // RATHER...
  //
  // Ordered by first occurrence of S in an (omega,S) subspace, then
  // by increasing L, subject to triangularity to target J value, but
  // without regard to whether or not L actually exists in the
  // branching of this (or any) (omega,S).  An empty subspace is
  // therefore possible.
  // 
  ////////////////////////////////////////////////////////////////
  
  typedef std::pair<int,HalfInt> LSPair;

  class StateLS;  // forward declaration (to permit use as "friend" of SubspaceU3S)

  class SubspaceLS
    : public basis::BaseSubspace<LSPair,std::tuple<int>>
    // SubspaceLabelsType (std::pair<int,HalfInt>): (L,S)
    // StateLabelsType (int): index into U3S space

    {
      public:

      // constructors

      SubspaceLS() {};
      // default constructor -- provided since required for certain
      // purposes by STL container classes (e.g., std::vector::resize)

      SubspaceLS(int L, HalfInt S,const SpaceU3S& u3s_space);

      // accessors
      int L() const{return std::get<0>(labels());}
      HalfInt S() const {return std::get<1>(labels());}

      int full_dimension() const {return full_dimension_;}
      // int sector_dim() const {return full_dimension();} // DEPRECATED in favor of full_dimension

      // diagnostic output
      std::string Str() const;

      // state indexing lookup -- DEPRECATED
      int sector_index(int state_index) const
      // DEPRECATED -- in favor of StateLS::starting_index()
      {
        return state_substate_offset_.at(state_index);
      }

      private:

      // state auxiliary data
      //
      // Note: These parallel arrays are getting cumbersome.  Perhaps
      // this information should be bundled into a struct.
      friend class StateLS;
      std::vector<int> state_substate_offset_;  // starting index, counting (gamma,upsilon) multiplicity
      std::vector<int> state_dimension_;  // number of substates, counting (gamma,upsilon) multiplicity
      // std::vector<int> state_gamma_max_;  // gamma_max
      // std::vector<int> state_kappa_max_;  // kappa_max
      // std::vector<u3shell::U3SPN> state_sigmaSPN_;  // Sp irrep symmetry labels
      std::vector<u3::U3> state_omega_;  // U(3) label
      
      // dimension, counting (kappa,gamma,upsilon) multiplicity of subspace
      int full_dimension_;
    };

  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class StateLS
    : public basis::BaseState<SubspaceLS>
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
    int L() const {return subspace().L();}
    HalfInt S() const {return subspace().S();}

    // supplemental data accessors
    int substate_offset() const
    // Provide offset of first substate into fully expanded listing of
    // substates in subspace.
    {
      return subspace().state_substate_offset_[index()];
    }
    int dimension() const
    // Provide number of substates of this composite state.
    {
      return subspace().state_dimension_[index()];
    }
    // int gamma_max() const
    // // Provide gamma multiplicity of Sp irrep.
    // {
    //   return subspace().state_gamma_max_[index()];
    // }
    // int kappa_max() const
    // // Provide kappa branching multiplicity of L irrep.
    // {
    //   return subspace().state_kappa_max_[index()];
    // }
    // u3shell::U3SPN sigmaSPN() const
    //   // Provide full symmetry labels (sigma,Sp,Sn,S) of Sp irrep.
    //   {
    //     return subspace().state_sigmaSPN_[index()];
    //   }
    // u3::U3 sigma() const
    //   // Extract sigma label of Sp irrep.
    //   {
    //     return sigmaSPN().U3();
    //   }
    u3::U3 omega() const
      // Provide full symmetry labels (sigma,Sp,Sn,S) of Sp irrep.
      {
        return subspace().state_omega_[index()];
      }

  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////
  class SpaceLS
    : public basis::BaseSpace<SubspaceLS>
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
  };


  ////////////////////////////////////////////////////////////////
  // Sector
  // Enumerates sectors LS
  ////////////////////////////////////////////////////////////////
  typedef std::pair<int,int> OperatorLabelsLS;  // (L0,S0)

  class SectorLabelsLS
  {
    public:
    // Need N0,x0,S0,kappa0,L0, rho0
    typedef std::tuple<int,int,int,HalfInt> KeyType;
    ////////////////////////////////////////////////////////////////
    // constructors
    ////////////////////////////////////////////////////////////////
    //default constructor
    inline SectorLabelsLS()
      :bra_index_(-1), ket_index_(-1){}

    // construction from labels
    inline 
      SectorLabelsLS(
          int bra_index, int ket_index, 
          const OperatorLabelsLS& tensor_labels
        )
      : bra_index_( bra_index), ket_index_(ket_index),tensor_labels_(tensor_labels)
    {}

    ////////////////////////////////////////////////////////////////
    // accessors
    ////////////////////////////////////////////////////////////////
    int bra_index() const {return bra_index_;}
    int ket_index() const {return ket_index_;}
    OperatorLabelsLS operator_labels() const {return tensor_labels_;}
    int L0() const {return tensor_labels_.first;}
    HalfInt S0() const {return tensor_labels_.second;}

    inline KeyType Key() const
    {
      return KeyType(bra_index_,ket_index_,tensor_labels_.first,tensor_labels_.second);
    }

    inline friend bool operator == (const SectorLabelsLS& sector1, const SectorLabelsLS& sector2)
    {
      return sector1.Key() == sector2.Key();
    }

    inline friend bool operator < (const SectorLabelsLS& sector1, const SectorLabelsLS& sector2)
    {
      return sector1.Key() < sector2.Key();
    }

    ////////////////////////////////////////////////////////////////
    // hashing
    ////////////////////////////////////////////////////////////////

    inline friend std::size_t hash_value(const SectorLabelsLS& sector)
    {
      boost::hash<SectorLabelsLS::KeyType> hasher;
      return hasher(sector.Key());
    }
    ////////////////////////////////////////////////////////////////
    // string conversion
    ////////////////////////////////////////////////////////////////

    std::string Str() const;

    private:
    int bra_index_, ket_index_, L0_;
    OperatorLabelsLS tensor_labels_;
  };


  void GenerateOperatorLabelsLS(const HalfInt& J0, std::vector<OperatorLabelsLS>& tensor_labels);
  // Generate a list of two-body tensor labels given the triangular restriction on L0,S0, and J0
  // and the two-body limit on S0<=2.


  void GetSectorsLS(
      const spncci::SpaceLS& space_bra, 
      const spncci::SpaceLS& space_ket,
      const std::vector<OperatorLabelsLS>& tensor_labels,
      std::vector<spncci::SectorLabelsLS>& sector_labels
    );
  // Generates a cache of SectorLabelsU3S from operator labels given in 
  // relative_tensor_rmes, which are U(1)xSU(3)xSU(2) tensors labeled
  // by (N0,x0,S0,kappa0,L0). 
  // 
  // space_bra (input) : space for bra used to define sectors
  // space_ket (input) : space for ket used to define sectors
  // relative_tensor_rmes (input) : container of rme labels keys and rme values
  //                                RelativeRMEsU3ST defined in upcoupling.h
  // u3_sectors (output) : container with SectorLabelsU3S keys and index values

  void 
    ContractAndRegroupLSJ(
        const HalfInt& Jp,const HalfInt& J0, const HalfInt& J,
        const spncci::SpaceU3S& u3s_space,
        const std::vector<spncci::SectorLabelsU3S>& source_sector_labels,
        const basis::MatrixVector& source_sectors,
        const spncci::SpaceLS& target_space_bra,
        const spncci::SpaceLS& target_space_ket,
        const std::vector<spncci::SectorLabelsLS>& target_sector_labels,
        basis::MatrixVector& target_sectors
      );

  // Sums over omega,omega', omega0, kappa', kappa, and kappa0 to obtain
  // LS reduced matrix elements 

  void
  ConstructOperatorMatrix(
    const spncci::SpaceLS& bra_source_space,
    const spncci::SpaceLS& ket_source_space,
    std::vector<spncci::SectorLabelsLS>& source_sector_labels,
    basis::MatrixVector& source_sectors,
    Eigen::MatrixXd& operator_matrix
    );
  

}  // namespace

#endif
