/****************************************************************
  indexing_u3st.h

  Defines relative, relative-cm, and two-body state indexing in U3ST
  coupling scheme.

  Language: C++11
                                 
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/11/16 (aem,mac): Created, styled after indexing_lsjt.h.
  5/12/16 (aem,mac): Restructure by U(3) subspace.

****************************************************************/

#ifndef INDEXING_U3ST_H_
#define INDEXING_U3ST_H_

#include <string>

#include "utilities/indexing_base.h"
#include "sp3rlib/u3.h"


namespace u3shell {

  ////////////////////////////////////////////////////////////////
  // two-body states in U3ST scheme
  ////////////////////////////////////////////////////////////////
  
  // subspace labels: (omega,S,T,g)
  //   where omega=(N,(lambda,mu)) is a U(3) label
  // state labels within subspace: (N1)
  //
  // Subspace labels must obey:
  //   -- N~g
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N1
  //   -- [N2 is implied by N1+N2=N]
  // and subject to:
  //   -- SU(3) coupling constraint (N1,0)x(N2,0)->(lambda,mu)
  //   -- antisymmetry constraint lambda+mu+S+T~1 if N1==N2
  //
  // This basis is for *distinguishable* particle states, subject to
  // the symmetry constraints which apply to *identical*-particle
  // states.  That is, although the antisymmetry constraint for identical
  // orbitals (omega+S+T~1) is applied to the quantum numbers when
  // enumerating the basis, the states
  //    
  //   |(N1,N2)...>  and  |(N2,N1)...>
  //
  // are still counted separately in the basis.  No lexicographical
  // ordering constraint N1<=N2 on the two single-particle states is
  // imposed.  The basis is therefore redundant.  This simplifies
  // implementation of double contractions over "particle 1" and
  // "particle 2" indices as matrix multiplication, without the
  // calling code having to worry about "swapping" single particle
  // states within the two-body state and applying the relevant phase
  // (~omega+S+T+g+1).
  //
  // Subpaces in given Nmax truncation:
  //    How do we generate the omega list?
  //    Option 1:
  //      for each N in 0..Nmax
  //        for each (N1,N2) pair with N1+N2=N  [N1=0..Nmax]
  //          enumerate (N1,0)x(N2,0) -> omega
  //      and collect set of distinct omegas
  //    Option 2:
  //      for each N in 0..Nmax
  //        for lambda in 0..N
  //          for mu in 0..floor(N/2)
  //      attempt to build subspace, but prune to 
  //      subspaces of nonzero dimension
  //    
  // Within the full space, subspaces are ordered by:
  //    -- increasing omega 
  //    -- increasing S (S=0,1)
  //    -- increasing T (T=0,1)
  //    -- [g is implied by omega]
  // 
  // Subspaces are pruned to those of nonzero dimension.

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class TwoBodySubspaceU3ST
    : public shell::BaseSubspace<std::tuple<u3::U3,int,int,int>,std::tuple<int>>
  // Subspace class for two-body states of given U(3)xSxT.
  //
  // SubspaceLabelsType (std::tuple): <omega, S, T, g>
  //   omega (u3::U3) :  U(3) label
  //   S (int) : spin
  //   T (int) : isospin
  //   g (int) : grade (parity % 2) 
  // StateLabelsType (std::tuple): <N1>
  //   N1 (int) : N1
  // Note: tuple of int (as opposed to plain int) is necessary for the two forms 
  // of the state constructor to be syntactically distinct.
  {
  public:

    // constructor
    TwoBodySubspaceU3ST (u3::U3 omega, int S, int T, int g);

    // accessors
    u3::U3 omega() const {return std::get<0>(labels_);}
    int S() const {return std::get<1>(labels_);}
    int T() const {return std::get<2>(labels_);}
    int g() const {return std::get<3>(labels_);}
    int N() const {return int(omega().N());}

    // diagnostic output
    std::string Str() const;

  private:

    // validation
    bool ValidLabels() const;

  };

  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class TwoBodyStateU3ST
    : public shell::BaseState<TwoBodySubspaceU3ST>
  // State class for two-body states of given U(3)xSxT.
  {
    
  public:

    // pass-through constructors

  TwoBodyStateU3ST(const SubspaceType& subspace, int index)
    // Construct state by index.
    : shell::BaseState<TwoBodySubspaceU3ST>(subspace, index) {}

  TwoBodyStateU3ST(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : shell::BaseState<TwoBodySubspaceU3ST> (subspace, state_labels) {}

    // pass-through accessors
    u3::U3 omega() const {return Subspace().omega();}
    int S() const {return Subspace().S();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}
    int N() const {return Subspace().N();}

    // state label accessors
    int N1() const {return std::get<0>(GetStateLabels());}
    int N2() const {return N()-N1();}

  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class TwoBodySpaceU3ST
    : public shell::BaseSpace<TwoBodySubspaceU3ST>
  // Space class for two-body states of given U(3)xSxT.
  {
    
  public:

    // constructor
    TwoBodySpaceU3ST(int Nmax);

    // accessors
    int Nmax() const {return Nmax_;}

    // diagnostic output
    std::string Str() const;

  private:
    // truncation
    int Nmax_;

  };

  ////////////////////////////////////////////////////////////////
  // sectors
  ////////////////////////////////////////////////////////////////
  
  struct OperatorLabelsU3ST
  // U(1)xSU(3)xSxT operators labels
  //
  // For use in selection rules for enumerating operator sectors.
  //
  // Note: The U(1)xSU(3) labels do *not* in general constitute a
  // valid U(3) label and thus cannot be stored in a u3::U3.  E.g.,
  // operators carrying N0=0 but an SU(3) irrep other than (0,0) are
  // perfectly well possible.
  {

  OperatorLabelsU3ST(int N0_, const u3::SU3& x0_, HalfInt S0_, HalfInt T0_, int g0_)
  : N0(N0_), x0(x0_), S0(S0_), T0(T0_), g0(g0_)
    {}


    int N0;
    u3::SU3 x0;
    HalfInt S0, T0;
    int g0;
  };

  class TwoBodySectorsU3ST
    : public shell::BaseSectors<TwoBodySpaceU3ST>
  {

  public:

    // constructor

    TwoBodySectorsU3ST(TwoBodySpaceU3ST& space, const OperatorLabelsU3ST& operator_labels);

    // Enumerates sector pairs connected by an operator of given tensorial and parity character.
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).
    // Sectors are included in "both directions", i.e., there is no
    // asumption of hermiticity or attempt to therefore only store
    // "half" the sectors.

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
