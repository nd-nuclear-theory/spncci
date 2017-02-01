/****************************************************************
  spncci_branching_u3lsj.h

  U(3)xLS and U(3)xLSJ layers of SpNCCI basis branching.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  1/31/17 (mac): Extracted from sp_basis.
    
****************************************************************/

#ifndef SPNCCI_BRANCHING_U3LSJ_H_
#define SPNCCI_BRANCHING_U3LSJ_H_

#include <unordered_map>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "spncci/spncci_basis.h"
#include "spncci/spncci_branching_u3s.h"
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
