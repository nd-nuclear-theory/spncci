/****************************************************************
recurrence_indexing.h

    Indexing for SpNCCI recurrence

    SpNCCIIrrep
    ->  sigma
        -> omega (states (n,rho) -> upsilon)

    
    SpNCCILGISpace()
    ->SpNCCISpinSpace() [sigma]
      -> SpNCCISpinSubspace() [S]
        -> SpNCCISpinState() [SpSn] -> gamma_max

    RecurrenceSpatialIndexing()
    ->BasisSpatialIndexing()
        ->OperatorSpatialIndexing()
            ->UnitTensorSpatialIndexing()


  Anna E. McCoy and Patrick J. Fasano
  Institute for Nuclear Theory  
  and University of Notre Dame

  SPDX-License-Identifier: MIT

  06/24/21 (aem) : Created. 
****************************************************************/

#ifndef RECURRENCE_INDEXING_H_
#define RECURRENCE_INDEXING_H_

#include <array>
#include <unordered_map>

// #include "basis/hypersector.h"
#include "basis/degenerate.h"

#include "am/halfint.h"
#include "lgi/lgi.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/sp3r.h"
// #include "spncci/spncci_common.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/u3spn_scheme.h"
#include "u3shell/unit_tensor_space_u3s.h"
// #include "u3shell/upcoupling.h"

namespace spncci
{
  ///////////////////////////////////////////////////////////////////////////////////////////
  //   SpNCCILGISpace()
  // ->SpNCCISpinSpace() [sigma]
  //   -> SpNCCISpinSubspace() [S]
  //     -> SpNCCISpinState() [SpSn] -> gamma_max
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  using SpSn = std::tuple<HalfInt,HalfInt>;
  
  class SpNCCISpinSubspace
    : public basis::BaseDegenerateSubspace<std::tuple<HalfInt>,SpSn>
    {
      public:
        SpNCCISpinSubspace() = default;
        
        SpNCCISpinSubspace(
          const HalfInt& S, 
          const MultiplicityTagged<SpSn>::vector& spin_vector
          );
        
        HalfInt S() const {return std::get<0>(labels());}
      
      private:

    };


 class SpNCCISpinState
    : public basis::BaseDegenerateState<SpNCCISpinSubspace>
  {

    public:
    // pass-through constructors

      SpNCCISpinState(const SubspaceType& subspace, int index)
        // Construct state by index.
        : basis::BaseDegenerateState<SpNCCISpinSubspace> (subspace, index) {}
    
      SpNCCISpinState(
          const SubspaceType& subspace,
          const typename SubspaceType::StateLabelsType& state_labels
        )
        // Construct state by reverse lookup on labels.
        : basis::BaseDegenerateState<SpNCCISpinSubspace>(subspace,state_labels)
        {}

    // pass-through accessors for subspace labels
    HalfInt S() const {return subspace().S();}

    // state label accessors
    HalfInt Sp() const {return std::get<0>(labels());}
    HalfInt Sn() const {return std::get<1>(labels());}
    private:
  };

class SpNCCISpinSpace
    : public basis::BaseSpace<SpNCCISpinSubspace,std::tuple<u3::U3>>
  {
    public:
    SpNCCISpinSpace() = default;
    
    SpNCCISpinSpace(
      const u3::U3& sigma, 
      const std::map<HalfInt,MultiplicityTagged<SpSn>::vector>& spin_map
      );

    u3::U3 sigma() const {return std::get<0>(labels());}
    
    private:
    u3::U3 sigma_;
  };


class SpNCCILGISpace
    : public basis::BaseSpace<SpNCCISpinSpace>
  {
  
    public:
    SpNCCILGISpace()=default;
    SpNCCILGISpace(const lgi::MultiplicityTaggedLGIVector& lgi_vector);
    // private:
  };


}  // namespace
#endif