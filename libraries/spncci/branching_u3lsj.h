/****************************************************************
  branching_u3lsj.h

  U(3)xLS and U(3)xLSJ layers of SpNCCI basis branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/31/17 (mac): Extracted from sp_basis (as spncci_branching_u3lsj).
  2/19/17 (mac): Rename to branching_u3lsj.
    
****************************************************************/

#ifndef SPNCCI_SPNCCI_BRANCHING_U3LSJ_H_
#define SPNCCI_SPNCCI_BRANCHING_U3LSJ_H_

#include <unordered_map>


#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "spncci/branching_u3s.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"  
#include "u3shell/upcoupling.h"
#include "lgi/lgi.h"

namespace spncci
{

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
  // States are indexed by a tuple of an int corresponding to 
  // a U3S subspace. 
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
    : public basis::BaseSubspace< std::pair<int,HalfInt>,std::tuple<int> >
    // Subspace class for two-body states of given SO(3)xS.
    //
    {
      public:

      // constructors

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
      // Look up table to find starting index of state in 
      std::map<int,int> sector_index_lookup_;
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
        const spncci::SpaceLS& source_space,
        std::vector<spncci::SectorLabelsLS>& source_sector_labels,
        basis::MatrixVector& source_sectors,
        Eigen::MatrixXd& operator_matrix
      );
  // TODO (mac): generalize to case where bra_J (and bra_space_ls) != ket_J (and ket_space_ls)

}  // namespace

#endif
