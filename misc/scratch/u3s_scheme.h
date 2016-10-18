/****************************************************************
  u3s_scheme.h

  Defines generic "skeleton" bookkeeping for U(3)xSU(3) subspaces, i.e.,
  labeled by i.e., N(lambda,mu)S.  No detailed enumeration is made of
  the states within a subspace.  This bookkeeping is intended for use
  in keeping track of, e.g., the U(3)xSU(2) subspaces which arise
  within an lsu3shell SU(3) coupled basis, and operator sectors
  between these subspaces.

  DEPRECATED: Written but immediately disfavored for u3spn_scheme.

  Language: C++11
                                 
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  9/6/16 (mac): Created, based on a subset of u3st_scheme.

****************************************************************/

#ifndef U3S_SCHEME_H_
#define U3S_SCHEME_H_

#include <string>

#include "basis/indexing.h"
#include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"

namespace u3shell {


  ////////////////////////////////////////////////////////////////
  // basis indexing in U3S scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (omega,S) = U3S
  //
  // The parity label g can implicitly be deduced from omega if needed.
  // We do not explicitly store it.
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // We will not actually enumerate states within a subspace.  We
  // provide a "dummy" int state label where a StateLabelType is
  // required by the basis template library.
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
  
  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class SubspaceU3S
    : public basis::BaseSubspace<u3::U3S,int>
  // Subspace class for two-body states of given U(3)xS.
  //
  // SubspaceLabelsType (u3::U3S): (omega,S)
  // StateLabelsType (int): 1 (just a place holder)
  {
    public:

    // constructor
    SubspaceU3S (const u3::U3S& omegaS, int dimension);

    // accessors
    u3::U3S U3S() const {return labels_;}
    u3::U3 omega() const {return U3S().U3();}
    u3::SU3 SU3() const {return U3S().SU3();}
    HalfInt N() const {return U3S().U3().N();}
    HalfInt S() const {return U3S().S();}

    // diagnostic output
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
    SpaceU3S(const std::map<u3::U3S,int>& subspace_dimensions);

    // diagnostic output
    std::string Str() const;

  private:

  };


  ////////////////////////////////////////////////////////////////
  // sectors
  ////////////////////////////////////////////////////////////////

  class SectorsU3S
    : public basis::BaseSectors<SpaceU3S>
  // U3S-scheme operator sectors.
  //
  // Sectors are enumerated in lexicographical order by (bra)(ket)(rho).
  {

  public:

    // constructor

    SectorsU3S(SpaceU3S& space, const OperatorLabelsU3S& operator_labels);
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).
    //
  };


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
