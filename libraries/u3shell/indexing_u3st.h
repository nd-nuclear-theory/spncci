/****************************************************************
  indexing_u3st.h

  Defines relative, relative-cm, and two-body state indexing in U3ST
  coupling scheme.

  Language: C++11
                                 
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  5/11/16 (aem,mac): Created, styled after indexing_lsjt.h.

****************************************************************/

#ifndef INDEXING_U3ST_H_
#define INDEXING_U3ST_H_

#include "utilities/indexing_base.h"
#include "sp3rlib/u3.h"


namespace u3shell {

  ////////////////////////////////////////////////////////////////
  // relative states in U3ST scheme
  ////////////////////////////////////////////////////////////////
  
  // Note: Relative subspaces in U3ST scheme are trivial (unit
  // dimension) so indexing is not actually implemented.  But we note
  // the allowed quantum numbers here for reference.

  // subspace labels: (Nr,S,T,g)    P=(-)^g
  // state labels within subspace: ()
  //
  // Within the space, subspaces are ordered by:
  //   -- increasing Nr (Nr=0,1,...,Nmax)
  //   -- increasing S (S=0,1)
  //   -- [T forced by Nr+S+T~1]
  //   -- [g forced by g~Nr]
  // subject to:
  //   -- parity constraint Nr~g
  //   -- antisymmetry constraint Nr+S+T~1
  //
  // Note that ordering of subspaces is therefore lexicographic by (N,S).

  ////////////////////////////////////////////////////////////////
  // relative-cm states in U3ST scheme
  ////////////////////////////////////////////////////////////////

  // subspace labels: (Ncm,x,S,T,g)
  // state labels within subspace: (Nr)
  //
  // Truncation of the subspace is by the two-body Nmax.
  //
  // Subspace subject to:
  //   -- antisymmetry constraint Nr+S+T~1 on states
  //      implies Ncm+S+T+g~1 for subspace
  //
  // Within a subspace, states are ordered by:
  //   -- increasing Nr
  //        or, equivalently, increasing N=Ncm+Nr (N~g)
  // and subject to:
  //   -- coupling constraint (Ncm,0)x(Nr,0)->x
  //   -- parity constraint N=Nr+Ncm~g
  //   -- [antisymmetry constraint Nr+S+T~1]
  // 
  // A full space (and ordering of subspaces) has not been
  // implemented, as it is not needed for the Moshinsky calculations.
  //
  // Note: Dimension of relative-cm subspace appears to always be either 0 or 1.

  //subspace
  
  class RelativeCMSubspaceU3ST
    : public shell::BaseSubspace< std::tuple<int,u3::SU3,int,int,int> , std::tuple<int> >
    {
    
    public:

      // constructor

      RelativeCMSubspaceU3ST(int Ncm, const u3::SU3& x, int S, int T, int g, int Nmax);

      // accessors

      int Ncm() const {return std::get<0>(labels_);}
      u3::SU3 SU3() const {return std::get<1>(labels_);}
      int S() const {return std::get<2>(labels_);}
      int T() const {return std::get<3>(labels_);}
      int g() const {return std::get<4>(labels_);}
      int Nmax() const {return Nmax_;}

    private:

      //validation
      bool ValidLabels() const;

      // truncation
      int Nmax_;

    };

  // state

  class RelativeCMStateU3ST
    : public shell::BaseState<RelativeCMSubspaceU3ST>
  {
    
  public:

    // pass-through constructors

  RelativeCMStateU3ST(const SubspaceType& subspace, int index)
    // Construct state by index.
    : shell::BaseState<RelativeCMSubspaceU3ST>(subspace, index) {}

  RelativeCMStateU3ST(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : shell::BaseState<RelativeCMSubspaceU3ST>(subspace, state_labels) {}

    // pass-through accessors
    int Ncm() const {return Subspace().Ncm();}
    u3::SU3 SU3() const {return Subspace().SU3();}
    int S() const {return Subspace().S();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}

    // state label accessors
    int Nr() const {return std::get<0>(GetStateLabels());}
    int N() const {return  Nr()+Ncm();}

  };

  // space -- not defined

  // sectors -- not defined



  ////////////////////////////////////////////////////////////////
  // two-body states in LSJT scheme
  ////////////////////////////////////////////////////////////////
  
  // subspace labels: (L,S,J,T,g)
  // state labels within subspace: (N1,l1,N2,l2)
  //
  // Truncation of the space is by the two-body Nmax.
  //
  // Within a subspace, the states are ordered by:
  //   -- increasing N (N~g)
  //   -- lexicographically increasing (N1,l1)
  //   -- lexicographically increasing (N2,l2)
  // and subject to:
  //   -- triangularity constraint on (l1,l2,L)
  //   -- parity constraint N~g
  //   -- antisymmetry constraint L+S+T~1 if (N1,l1)==(N2,l2)
  //
  // Note that, although the antisymmetry constraint for identical
  // orbitals (L+S+T~1) is applied to the quantum numbers when
  // enumerating the basis, the states
  //    
  //   |(N1,l1)(N2,l2)...>  and  |(N2,l2)(N1,l1)...>
  //
  // are still counted separately in the basis.  That is, *no*
  // lexicographical ordering constraint (N1,l1)<=(N2,l2) on the two
  // single-particle states is imposed.  The basis is therefore
  // redundant.  This simplifies implementation of double contractions
  // over "particle 1" and "particle 2" indices as matrix
  // multiplication, without the calling code having to worry about
  // "swapping" single particle states within the two-body state and
  // applying the relevant phase (~L+S+T+g+1).
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing L (L=0,1,...,Nmax)
  //    -- increasing S (S=0,1)
  //    -- increasing J
  //    -- increasing T (T=0,1)
  //    -- increasing g (g=0,1)
  // subject to:
  //    -- triangularity of (L,S,J)
  // 
  // Subspaces are pruned to those of nonzero dimension.
  //
  // Note that ordering of subspaces is lexicographic by (L,S,J,g).

  //subspace

#if 0
  class TwoBodySubspaceLSJT
    : public shell::BaseSubspace< std::tuple<int,int,int,int,int> , std::tuple<int,int,int,int> >
  {
    
  public:

    // constructor

    TwoBodySubspaceLSJT (int L, int S, int J, int T, int g, int Nmax);
    // Set up indexing in Nmax truncation.

    // accessors

    int L() const {return std::get<0>(labels_);}
    int S() const {return std::get<1>(labels_);}
    int J() const {return std::get<2>(labels_);}
    int T() const {return std::get<3>(labels_);}
    int g() const {return std::get<4>(labels_);}
    int Nmax() const {return Nmax_;}

  private:

    //validation
    bool ValidLabels() const;

    // truncation
    int Nmax_;
  };

  // state

  class TwoBodyStateLSJT
    : public shell::BaseState<TwoBodySubspaceLSJT>
  {
    
  public:

    // pass-through constructors

  TwoBodyStateLSJT(const SubspaceType& subspace, int index)
    // Construct state by index.
    : shell::BaseState(subspace, index) {}

  TwoBodyStateLSJT(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : shell::BaseState (subspace, state_labels) {}

    // pass-through accessors
    int L() const {return Subspace().L();}
    int S() const {return Subspace().S();}
    int J() const {return Subspace().J();}
    int T() const {return Subspace().T();}
    int g() const {return Subspace().g();}

    // state label accessors
    int N1() const {return std::get<0>(GetStateLabels());}
    int l1() const {return std::get<1>(GetStateLabels());}
    int N2() const {return std::get<2>(GetStateLabels());}
    int l2() const {return std::get<3>(GetStateLabels());}

    int N() const {return  N1()+N2();}

  };

  // space

  class TwoBodySpaceLSJT
    : public shell::BaseSpace<TwoBodySubspaceLSJT>
  {
    
  public:

    // constructor

    TwoBodySpaceLSJT(int Nmax);
    // Enumerates all relative LSJT subspaces of given dimension up to a
    // given Nmax cutoff.

    // diagnostic output
    void Print(std::ostream& os) const;

    // accessors
    int Nmax() const {return Nmax_;}

  private:
    // truncation
    int Nmax_;


  };

  // sectors

  class TwoBodySectorsLSJT
    : public shell::BaseSectors<TwoBodySpaceLSJT>
  {

  public:

    // constructor

    TwoBodySectorsLSJT(TwoBodySpaceLSJT& space, int J0, int g0);

    // Enumerates sector pairs connected by an operator of given tensorial and parity character.
    //
    // Sectors are enumerated in lexicographical order by (bra)(ket).
    // Sectors are included in "both directions", i.e., there is no
    // asumption of hermiticity or attempt to therefore only store
    // "half" the sectors.

  };


#endif


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace

#endif
