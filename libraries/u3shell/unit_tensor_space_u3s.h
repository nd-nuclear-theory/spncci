/****************************************************************
  unit_tensor_space_u3s.h

  Defines unit tensor subspaces grouped by u3s operator labels

  Language: C++11

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  4/4/17 (aem,mac): Created.
  1/31/18 (aem) : Add ObservableSpaceU3S for relative observables
    decomposed by upcouling
****************************************************************/

// TODO: Rename to operator_space_u3s

#ifndef UNIT_TENSOR_SPACE_U3ST_H_
#define UNIT_TENSOR_SPACE_U3ST_H_

#include <string>

#include "basis/basis.h"
#include "basis/hypersector.h"
#include "sp3rlib/u3.h"
#include "u3shell/tensor_labels.h"

namespace u3shell {

  typedef std::tuple<u3::SU3,HalfInt,int,int> UnitTensorSubspaceLabels;
  typedef std::set<UnitTensorSubspaceLabels> UnitTensorSubspaceLabelsSet;

  ////////////////////////////////////////////////////////////////
  // relative unit tensors subspaces in U3S scheme
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (x0,S0,etap,eta)
  //
  // state labels : (T0,Sp,Tp,S,T)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing N0
  //    -- increasing etap
  //    -- increasing x0
  //    -- increasing S0 (S0=0,2)
  //    -- [g is implied by omega (N~g)] note this is the relative g
  // and subject to:
  //   -- N~g
  //
  ////////////////////////////////////////////////////////////////

  class RelativeUnitTensorSubspaceU3S;
  class RelativeUnitTensorStateU3S;
  class RelativeUnitTensorSpaceU3S;

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////
  typedef std::tuple<int,int,int,int,int> RelativeUnitTensorStateLabelsU3S;

  class RelativeUnitTensorSubspaceU3S
    : public basis::BaseSubspace<RelativeUnitTensorSubspaceU3S,
        std::tuple<u3::SU3,HalfInt,int,int>,
        RelativeUnitTensorStateU3S,
        std::tuple<int,int,int,int,int>>
    // Subspace class for two-body states of given S(3)xS.
    //
    // SubspaceLabelsType (std::tuple): <x0, S0,etap,eta>
    //   x0   (SU3) : SU(3) label
    //   S0   (HalfInt) : spin
    //   etap (int) : oscillator quanta in bra
    //   eta  (int) : oscillator quanta in ket

    // StateLabelsType (std::tuple): <T0,Sp,Tp,S,T>
    {
  public:

    // constructor
    RelativeUnitTensorSubspaceU3S (
      u3::SU3 x0, HalfInt S0, int etap, int eta,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels
     );

    // accessors
    u3::SU3 x0() const {return std::get<0>(labels());}
    HalfInt S0() const {return std::get<1>(labels());}
    int etap() const {return std::get<2>(labels());}
    int eta() const {return std::get<3>(labels());}
    int N0() const {return etap()-eta();}

    std::tuple<u3::SU3,HalfInt,int,int> Key() const {return labels();}
    // diagnostic output
    std::string Str() const;
    std::string LabelStr() const;

  private:

  };

  ////////////////////////////////////////////////////////////////
  // state
  ////////////////////////////////////////////////////////////////

  class RelativeUnitTensorStateU3S
    : public basis::BaseState<RelativeUnitTensorSubspaceU3S>
  // State class for relative unit tensor states of given U(3)xS subspaces.
  {

  public:

    // pass-through constructors

  RelativeUnitTensorStateU3S(const SubspaceType& subspace, int index)
    // Construct state by index.
    : basis::BaseState<RelativeUnitTensorSubspaceU3S>(subspace, index) {}

  RelativeUnitTensorStateU3S(const SubspaceType& subspace, const typename SubspaceType::StateLabelsType& state_labels)
    // Construct state by reverse lookup on labels.
    : basis::BaseState<RelativeUnitTensorSubspaceU3S> (subspace, state_labels) {}

    // pass-through accessors
    u3::SU3 x0() const {return subspace().x0();}
    HalfInt S0() const {return subspace().S0();}
    int etap() const {return subspace().etap();}
    int eta() const {return subspace().eta();}
    int N0() const {return etap()-eta();}

    // state label accessors
    int T0() const {return std::get<0>(labels());}
    int Sp() const {return std::get<1>(labels());}
    int Tp() const {return std::get<2>(labels());}
    int S() const {return std::get<3>(labels());}
    int T() const {return std::get<4>(labels());}

  };

  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class RelativeUnitTensorSpaceU3S
    : public basis::BaseSpace<RelativeUnitTensorSubspaceU3S>
  // Space class for relative unit tensor states of given U(3)xS.
  {

  public:

    // constructor
    RelativeUnitTensorSpaceU3S(
      int Nmax, int N1v,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels
    );

    RelativeUnitTensorSpaceU3S(
      int Nmax, int N1v,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
      const std::map< std::pair<int,int>, UnitTensorSubspaceLabelsSet>&
        NnpNn_organized_unit_tensor_subspaces
    );

    // accessors
    int Nmax() const {return Nmax_;}

    // diagnostic output
    std::string Str() const;

  private:
    // truncation
    int Nmax_;
    int N1v_;

  };





  ////////////////////////////////////////////////////////////////
  // relative tensors subspaces in U3S scheme
  ////////////////////////////////////////////////////////////////
  typedef std::tuple<int,u3::SU3,HalfInt,int,int> ObservableSubspaceLabels;

  ////////////////////////////////////////////////////////////////
  //
  // Labeling
  //
  // subspace labels: (x0,S0,kappa0,L0)
  //
  // state labels : (dummy)
  //
  ////////////////////////////////////////////////////////////////
  //
  // Subspaces
  //
  // Within the full space, subspaces are ordered by:
  //    -- increasing N0
  //    -- increasing x0
  //    -- increasing S0 (S0=0,2)
  //    -- increasing k0
  //    -- increasing L0
  //    -- [g is implied by omega (N~g)] note this is the relative g
  // and subject to:
  //   -- N~g
  //
  ////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////
  // subspace
  ////////////////////////////////////////////////////////////////
  class ObservableSubspaceU3S
    : public basis::BaseSubspace<ObservableSubspaceU3S,
    std::tuple<int,u3::SU3,HalfInt,int,int>,
    basis::BaseState<ObservableSubspaceU3S>,
    int>
    // Subspace class for two-body states of given S(3)xS.
    //
    // SubspaceLabelsType (std::tuple): <x0, S0,etap,eta>
    //   N0   (int) : U(1) label
    //   x0   (SU3) : SU(3) label
    //   S0   (HalfInt) : spin
    //   kappa0 (int) : index
    //   L0  (int) : index

    // StateLabelsType (int): dummy variable
    {
  public:

    // constructor
    ObservableSubspaceU3S (
      int N0, u3::SU3 x0, HalfInt S0, int kappa0, int L0
     );

    // accessors
    int N0() const {return std::get<0>(labels());}
    u3::SU3 x0() const {return std::get<1>(labels());}
    HalfInt S0() const {return std::get<2>(labels());}
    int kappa0() const {return std::get<3>(labels());}
    int L0() const {return std::get<4>(labels());}

    std::tuple<int, u3::SU3,HalfInt,int,int> Key() const {return labels();}
    // diagnostic output
    std::string Str() const;
    std::string LabelStr() const;

  private:

  };


  ////////////////////////////////////////////////////////////////
  // space
  ////////////////////////////////////////////////////////////////

  class ObservableSpaceU3S
    : public basis::BaseSpace<ObservableSubspaceU3S>
  // Space class for relative unit tensor states of given U(3)xS.
  {

  public:

    // constructor
    inline ObservableSpaceU3S()
      {}

    // constructor
    ObservableSpaceU3S(
      const std::vector<u3shell::IndexedOperatorLabelsU3S>& observable_labels
    );


    // diagnostic output
    std::string Str() const;

  private:

  };






} // namespace

#endif
