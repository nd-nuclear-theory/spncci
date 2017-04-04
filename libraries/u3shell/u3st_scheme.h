/****************************************************************
  u3st_scheme.h

  Defines relative, relative-cm, and two-body state indexing in U3ST
  coupling scheme.

  Language: C++11
                                 
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/11/16 (aem,mac): Created, styled after indexing_lsjt.h.
  5/12/16 (aem,mac): Restructure by U(3) subspace.
  5/14/16 (mac): Implement sectors.  Impose canonical ordering on two-body states.
  5/25/16 (mac): Add "all-to-all" sectors constructors.
  6/8/16 (mac): Rename to u3st_scheme.h.

****************************************************************/

#ifndef U3ST_SCHEME_H_
#define U3ST_SCHEME_H_

#include <string>

#include "basis/basis.h"
#include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"


namespace u3shell {


  ////////////////////////////////////////////////////////////////
  // relative states in U3ST scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (N,S,T,g)
  //   where the SU(3) character of the state is given by omega=(N,(N,0))
  //
  // In each case there is only one state in the subspace
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing N 
  //    -- increasing S (S=0,1)
  //    -- increasing T (T=0,1)
  //    -- [g is implied by omega (N~g)] note this is the relative g
  // and subject to:
  //   -- N~g
  // 
  // Subspaces are pruned to those of nonzero dimension.
  // Since spaces of interest are those which transform to two-body operators,
  // we restrict the basis to states which satisfy
  //    -- N+S+T odd
  //
  ////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class RelativeSubspaceU3ST
    : public basis::BaseSubspace<std::tuple<int,int,int,int>,int>
    // Subspace class for two-body states of given U(3)xSxT.
    //
    // SubspaceLabelsType (std::tuple): <N, S, T>
    //   N (int) :  U(3) label
    //   S (int) : spin
    //   T (int) : isospin
    //   g (int) : grade (parity % 2) 
    // StateLabelsType (int): 1 (just a place holder)
    {
  public:

    // constructor
    RelativeSubspaceU3ST (int N, int S, int T, int g);

    // accessors
    u3::U3 omega() const {return u3::U3(std::get<0>(labels_),0,0);}
    u3::SU3 SU3() const {return u3::SU3(std::get<0>(labels_),0);}
    int N() const {return std::get<0>(labels_);}
    int S() const {return std::get<1>(labels_);}
    int T() const {return std::get<2>(labels_);}
    int g() const {return std::get<3>(labels_);}
    std::tuple<int,int,int,int> Key() const {return labels_;}
    // diagnostic output
    std::string Str() const;

  private:

    // validation
    bool ValidLabels() const;

  };

  // ////////////////////////////////////////////////////////////////
  // // state
  // ////////////////////////////////////////////////////////////////

  // class TwoBodyStateU3ST
  //   : public basis::BaseState<TwoBodySubspaceU3ST>
  // // State class for two-body states of given U(3)xSxT.
  // {
    
  // public:

  //   // pass-through constructors

  // TwoBodyStateU3ST(const SubspaceType& subspace, int index)
  //   // Construct state by index.
  //   : basis::BaseState<TwoBodySubspaceU3ST>(subspace, index) {}

  // TwoBodyStateU3ST(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
  //   // Construct state by reverse lookup on labels.
  //   : basis::BaseState<TwoBodySubspaceU3ST> (subspace, state_labels) {}

  //   // pass-through accessors
  //   u3::U3 omega() const {return Subspace().omega();}
  //   int S() const {return Subspace().S();}
  //   int T() const {return Subspace().T();}
  //   int g() const {return Subspace().g();}
  //   int N() const {return Subspace().N();}

  //   // state label accessors
  //   int N1() const {return std::get<0>(GetStateLabels());}
  //   int N2() const {return std::get<1>(GetStateLabels());}
  //   // equivalently, int N2() const {return N()-N1();}

  // };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class RelativeSpaceU3ST
    : public basis::BaseSpace<RelativeSubspaceU3ST>
  // Space class for two-body states of given U(3)xSxT.
  {
    
  public:

    // constructor
    RelativeSpaceU3ST(int Nmax);

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

  class RelativeSectorsU3ST
    : public basis::BaseSectors<RelativeSpaceU3ST>
  // U3ST-scheme relative sectors.
  //
  // Sectors are enumerated in lexicographical order by (bra)(ket).
  {

  public:

    // constructor

    RelativeSectorsU3ST(RelativeSpaceU3ST& space);
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //

    RelativeSectorsU3ST(RelativeSpaceU3ST& space, const OperatorLabelsU3ST& operator_labels);
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).
    //
  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // Relative-CM states in U3ST scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////  
  //
  // Labeling
  //
  // subspace labels: (omega,S,T,g)
  //   where omega=(N,(lambda,mu)) is a U(3) label
  //
  // state labels within subspace: (Nr,Ncm)
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
  // See documentation for TwoBodySubspaceU3ST below. 
  ////////////////////////////////////////////////////////////////
  //
  // States
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing Nr
  //   -- [Ncm is implied by Nr+Ncm=N]
  // and subject to:
  //   -- SU(3) coupling constraint (Nr,0)x(Ncm,0)->(lambda,mu)
  //   -- antisymmetry constraint Nr+S+T~1
  //
  // This basis is for *identical* particle states, as enforced by the
  // antisymmetry constraint on Nr.
  //
  ////////////////////////////////////////////////////////////////  

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////

  class RelativeCMSubspaceU3ST
    : public basis::BaseSubspace<std::tuple<u3::U3,int,int,int>,std::tuple<int,int>>
    // Subspace class for two-body states of given U(3)xSxT.
    //
    // SubspaceLabelsType (std::tuple): <omega, S, T, g>
    //   omega (u3::U3) :  U(3) label
    //   S (int) : spin
    //   T (int) : isospin
    //   g (int) : grade (parity % 2) 
    // StateLabelsType (std::tuple): <Nr,Ncm>
    //   Nr (int) : Nr
    //   Ncm (int) : Ncm (=N-Nr, with N determined by omega) -- redundant label
    {
  public:

    // constructor
    RelativeCMSubspaceU3ST (u3::U3 omega, int S, int T, int g);

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

  class RelativeCMStateU3ST
    : public basis::BaseState<RelativeCMSubspaceU3ST>
  // State class for two-body states of given U(3)xSxT.
  {
    
  public:

    // pass-through constructors

  RelativeCMStateU3ST(const SubspaceType& subspace, int index)
    // Construct state by index.
    : basis::BaseState<RelativeCMSubspaceU3ST>(subspace, index) {}

  RelativeCMStateU3ST(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : basis::BaseState<RelativeCMSubspaceU3ST> (subspace, state_labels) {}

    // pass-through accessors
    u3::U3 omega() const {return Subspace().omega();}
    int S() const {return Subspace().S();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}
    int N() const {return Subspace().N();}

    // state label accessors
    int Nr() const {return std::get<0>(GetStateLabels());}
    int Ncm() const {return std::get<1>(GetStateLabels());}
    // equivalently, int N2() const {return N()-N1();}

  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class RelativeCMSpaceU3ST
    : public basis::BaseSpace<RelativeCMSubspaceU3ST>
  // Space class for two-body states of given U(3)xSxT.
  {
    
  public:

    // constructor
    RelativeCMSpaceU3ST(int Nmax);

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

  class RelativeCMSectorsU3ST
    : public basis::BaseSectors<RelativeCMSpaceU3ST>
  // U3ST-scheme relative-cm sectors.
  //
  // Sectors are enumerated in lexicographical order by (bra)(ket).
  // Sectors are included in "both directions", i.e., there is no
  // asumption of hermiticity or attempt to therefore only store
  // "half" the sectors.
  {

  public:

    // constructor

    RelativeCMSectorsU3ST(RelativeCMSpaceU3ST& space);
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Since no operator SU(3) is selected, there is no meaning for
    // multiplicity, and unit multiplicity is assumed.  This is useful
    // if we *subsequently* wish to iterate over allowed operators on
    // this sector.

    RelativeCMSectorsU3ST(RelativeCMSpaceU3ST& space, const OperatorLabelsU3ST& operator_labels);
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).
    //
    // Multiplicity is determined based on the given operator SU(3)
    // character, and sectors with multiplicity indices
    // 1..multiplicity are defined.  This is useful if we
    // wish to iterate over allowed sectors for a given operator.

  };
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////



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
  // state labels within subspace: (N1,N2)
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
  //   -- The labels are subject to the antisymmetry constraint (omega+S+T~1)
  //      if N1==N2.
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
    : public basis::BaseSubspace<std::tuple<u3::U3,int,int,int>,std::tuple<int,int>>
    // Subspace class for two-body states of given U(3)xSxT.
    //
    // SubspaceLabelsType (std::tuple): <omega, S, T, g>
    //   omega (u3::U3) :  U(3) label
    //   S (int) : spin
    //   T (int) : isospin
    //   g (int) : grade (parity % 2) 
    // StateLabelsType (std::tuple): <N1,N2>
    //   N1 (int) : N1
    //   N2 (int) : N2 (=N-N1, with N determined by omega) -- redundant label
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
    : public basis::BaseState<TwoBodySubspaceU3ST>
  // State class for two-body states of given U(3)xSxT.
  {
    
  public:

    // pass-through constructors

  TwoBodyStateU3ST(const SubspaceType& subspace, int index)
    // Construct state by index.
    : basis::BaseState<TwoBodySubspaceU3ST>(subspace, index) {}

  TwoBodyStateU3ST(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : basis::BaseState<TwoBodySubspaceU3ST> (subspace, state_labels) {}

    // pass-through accessors
    u3::U3 omega() const {return Subspace().omega();}
    int S() const {return Subspace().S();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}
    int N() const {return Subspace().N();}

    // state label accessors
    int N1() const {return std::get<0>(GetStateLabels());}
    int N2() const {return std::get<1>(GetStateLabels());}
    // equivalently, int N2() const {return N()-N1();}

  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class TwoBodySpaceU3ST
    : public basis::BaseSpace<TwoBodySubspaceU3ST>
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
    : public basis::BaseSectors<TwoBodySpaceU3ST>
  // U3ST-scheme two-body sectors.
  //
  // Sectors are enumerated in lexicographical order by (bra)(ket).
  // Sectors are included in "both directions", i.e., there is no
  // asumption of hermiticity or attempt to therefore only store
  // "half" the sectors.
  {

  public:

    // constructor

    TwoBodySectorsU3ST(TwoBodySpaceU3ST& space);
    // Enumerate all sector pairs ("all-to-all" sector enumeration).
    //
    // Since no operator SU(3) is selected, there is no meaning for
    // multiplicity, and unit multiplicity is assumed.  This is useful
    // if we *subsequently* wish to iterate over allowed operators on
    // this sector.

    TwoBodySectorsU3ST(TwoBodySpaceU3ST& space, const OperatorLabelsU3ST& operator_labels);
    // Enumerate sector pairs connected by an operator of given
    // tensorial and parity character ("constrained" sector
    // enumeration).
    //
    // Multiplicity is determined based on the given operator SU(3)
    // character, and sectors with multiplicity indices
    // 1..multiplicity are defined.  This is useful if we
    // wish to iterate over allowed sectors for a given operator.

  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
