/****************************************************************
  indexing_u3st.h

  Defines relative, relative-cm, and two-body state indexing in U3ST
  coupling scheme.

  Language: C++11
                                 
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/11/16 (aem,mac): Created, styled after indexing_lsjt.h.
  5/12/16 (aem,mac): Restructure by U(3) subspace.
  5/14/16 (mac): Implement sectors.  Impose canonical ordering on two-body states.

****************************************************************/

#ifndef INDEXING_U3ST_H_
#define INDEXING_U3ST_H_

#include <string>

#include "utilities/indexing_base.h"
#include "sp3rlib/u3.h"
#include "u3shell/tensor.h"


namespace u3shell {

  ////////////////////////////////////////////////////////////////
  // two-body states in U3ST scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (omega,S,T,g)
  //   where omega=(N,(lambda,mu)) is a U(3) label
  //
  // state labels within subspace: (N1)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing omega 
  //    -- increasing S (S=0,1)
  //    -- increasing T (T=0,1)
  //    -- [g is implied by omega (N~g)]
  // and subject to:
  //   -- N~g
  // 
  // Subspaces are pruned to those of nonzero dimension.
  //
  // How do we generate the omega list under a given Nmax truncation?
  //
  //   Option 1:
  //     for each N in 0..Nmax
  //       for each (N1,N2) pair with N1+N2=N  [N1=0..Nmax]
  //         enumerate (N1,0)x(N2,0) -> omega
  //     and collect set of distinct omegas
  //
  //     This requires duplicative coding and possibly iteration over
  //     large numbers of duplicative set elements.
  //
  //   Option 2 (ADOPTED):
  //     for each N in 0..Nmax
  //       for lambda in 0..N
  //         for mu in 0..floor(N/2)
  //     attempt to build subspace, but prune to 
  //     subspaces of nonzero dimension
  //
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N1
  //   -- [N2 is implied by N1+N2=N]
  // and subject to:
  //   -- SU(3) coupling constraint (N1,0)x(N2,0)->(lambda,mu)
  //   -- antisymmetry constraint omega+S+T~1 if N1==N2
  //
  // Here, by the phase contribution from omega, we mean lambda+mu.
  //
  // This basis is for *identical* particle states:
  //   -- The labels are subject to the antisymmetry constraint (omega+S+T~1).
  //   -- A canonical (lexicographic) ordering constraint is applied to the 
  //      single-particle quantum numbers.  That is, when
  //      enumerating the basis, the states
  //    
  //        |(N1,N2)...>  and  |(N2,N1)...>
  //
  //      would be redundant, and only the first (for N1<=N2) is
  //      retained.
  //
  ////////////////////////////////////////////////////////////////  

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
