/****************************************************************
  branching_u3s.h

  U(3)xS layer of SpNCCI basis branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/31/17 (mac): Extracted from sp_basis (as spncci_branching_u3s).
  2/17/17 (mac): Extract BabySpNCCISubspace to spncci_basis.
  2/19/17 (mac): Rename to branching_u3s.
****************************************************************/

#ifndef SPNCCI_SPNCCI_BRANCHING_U3S_H_
#define SPNCCI_SPNCCI_BRANCHING_U3S_H_

#include <unordered_map>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "spncci/unit_tensor.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"  
#include "u3shell/upcoupling.h"
#include "lgi/lgi.h"

namespace spncci
{

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
  // States are indexed by a tuple of an int (to distinguish it from an 
  // alternative constructor) which corresponds to the index of the
  // corresponding baby spncci subspace. 
  // The lgi family index and dimension (nu_max*gamma_max) as well as gamma_max
  // and nu_max can be extracted from baby spncci. 
  //
  // Associated with each subspace is a look-up table which can look up
  // the starting index of the particular "state" in the U3S sector.  
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
    : public basis::BaseSubspace<u3::U3S, std::tuple<int> >
  // Subspace class for two-body states of given U(3)xS.
  //
  // SubspaceLabelsType (u3shell::U3S) : (omega,S)
  // StateLabelsType (int) : Baby Spncci index
  {
    public:
    // constructor
    SubspaceU3S() {};
    // default constructor -- provided since required for certain
    // purposes by STL container classes (e.g., std::vector::resize)

    SubspaceU3S (const u3::U3S& omegaS,const BabySpNCCISpace& baby_spncci_space);

    // accessors
    u3::U3S omegaS() const {return labels_;}
    u3::U3 omega() const {return omegaS().U3();}
    u3::SU3 x() const {return omegaS().SU3();}
    HalfInt N() const {return omegaS().U3().N();}
    HalfInt S() const {return omegaS().S();}
    int sector_dim() const {return sector_size_;}
    // diagnostic output
    std::string Str() const;

    int sector_index(int state_index) const
      {
        int sector_index=-1;
        for(auto it=sector_index_lookup_.begin(); it!=sector_index_lookup_.end(); ++it)
          {
            if(it->first==state_index)
              {
                sector_index=it->second;
                return sector_index;
              }
          }
        // if none, found, then return -1.
        return sector_index;
      }

    private:
      int sector_size_;
      // Look up table to find starting index of state in sectors
      std::map<int,int> sector_index_lookup_;

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
  
  StateU3S(const SubspaceType& subspace, int& index)
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


    // // state label accessors
    // int gamma() const {return std::get<0>(GetStateLabels());} 
    // u3::U3 sigma() const {return std::get<1>(GetStateLabels());}

    // int sector_dim() const {return std::get<2>(GetStateLabels());}


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

    SpaceU3S(const BabySpNCCISpace& baby_spncci_space);

    // diagnostic output
    std::string Str() const;

  private:
    int dimension_;
  };





  ////////////////////////////////////////////////////////////////
  // Sector
  // Enumerates omegaS sectors
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
 

  void 
  ContractAndRegroupU3S(
      int Nmax, int N1b,
      const std::vector<spncci::SectorLabelsU3S>& sector_labels_vector,
      const u3shell::RelativeRMEsU3ST& interaction_rme_cache,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      spncci::SpaceU3S& target_space,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_sector_cache,
      basis::MatrixVector& matrix_vector
    );
  // Args:
  //  Nmax (input) : Basis truncation parameter
  //  N1b (input) : Basis single particle cutoff for Nmax=0
  //  sector_labels_vector (input) : vector of sector labels 
  //  interaction_rme_cache (input) : Container holding interaction rme's keyed
  //     by RelativeUnitTensorU3ST labels
  //  space (input) : space of omegaS subspaces
  //  unit_tensor_sector_cache (input) : nested container holding unit tensor rmes
  //  matrix_vector (output) : vector of U3S sectors indexed by U3S labels and kappa0,L0


}  // namespace

#endif
