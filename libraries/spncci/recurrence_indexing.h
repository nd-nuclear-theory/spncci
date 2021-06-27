/****************************************************************
recurrence_indexing.h

    Indexing for SpNCCI recurrence

    SpNCCIIrrep
    ->  sigma
        -> omega (states (n,rho) -> upsilon)

    
    spncci::spatial::Space() []
    ->spncci::spatial::LGISpace() [sigma]
        ->spncci::spatial::U3Subspace() [omega]
          ->spncci::spatial::U3State() [n,rho] n=Nn(lambda_n,mu_n)/(nx,ny,nz)

    spin::Space() []
    ->spin::LGISpace() [sigma]
      ->spin::SpinSubspace() [S]
        ->spin::[Iso]SpinState() [SpSn]/[T] -> gamma_max

    spatial::RecurrenceSpace() []
    ->spatial::RecurrenceLGISpace() [sigma,sigma']
      ->spatial::RecurrenceNnsumSpace() [Nsum]
        ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
          ->spatial::RecurrenceOperatorSubspace() [x0] ->rho0_max
            ->spatial::RecurrenceOperatorState() [Nbar,Nbar']

    spin::RecurrenceLGISpace() [sigma,sigma']
    ->spin::RecurrenceSpinSpace() [S,S']
      -> spin::RecurrenceS[PN/T]Subspace() [Sp,Sn,Sp',Sn']/[T,T']->(gamma,gamma')
        ->spin::RecurrenceOperatorState() [operator_index]


  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory  
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

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
namespace spin
{
  ///////////////////////////////////////////////////////////////////////////////////////////
  // spin::Space() []
  // ->spin::LGISpace() [sigma]
  //   ->spin::SpinSubspace() [S]
  //     ->spin::[Iso]SpinState() [SpSn]/[T] -> gamma_max
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


 class SpinState
    : public basis::BaseDegenerateState<SpNCCISpinSubspace>
  {

    public:
    // pass-through constructors

      SpinState(const SubspaceType& subspace, int index)
        // Construct state by index.
        : basis::BaseDegenerateState<SpNCCISpinSubspace> (subspace, index) {}
    
      SpinState(
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

  class LGISpace
      : public basis::BaseSpace<SpNCCISpinSubspace,std::tuple<u3::U3>>
    {
      public:
      LGISpace() = default;
      
      LGISpace(
        const u3::U3& sigma, 
        const std::map<HalfInt,MultiplicityTagged<SpSn>::vector>& spin_map
        );

      u3::U3 sigma() const {return std::get<0>(labels());}
      
    };


  class Space
      : public basis::BaseSpace<LGISpace>
    {
    
      public:
      Space()=default;
      Space(const lgi::MultiplicityTaggedLGIVector& lgi_vector);

      HalfInt Nsigma0() const {return Nsigma0_;}

    HalfInt Nsigma0_;

    };


  }  // namespace spin
  namespace spatial
  {
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // spncci::spatial::Space() []
    // ->spncci::spatial::LGISpace() [sigma]
    //     ->spncci::spatial::U3Subspace() [omega]
    //       ->spncci::spatial::U3State() [n,rho] n=Nn(lambda_n,mu_n)/(nx,ny,nz)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////// 

  class U3Subspace
      : public basis::BaseSubspace<std::tuple<u3::U3>,MultiplicityTagged<u3::U3>>
      {
        public:
          U3Subspace() = default;
          
          U3Subspace(
            const u3::U3& omega, 
            const MultiplicityTagged<u3::U3>::vector& nrho_vector
            );
          
          U3Subspace(const sp3r::U3Subspace& u3subspace);

          u3::U3 omega() const {return std::get<0>(labels());}
        
        // private:
        
      };

  class U3State
      : public basis::BaseDegenerateState<U3Subspace>
    {

      public:
      // pass-through constructors

        U3State(const SubspaceType& subspace, int index)
          // Construct state by index.
          : basis::BaseDegenerateState<U3Subspace> (subspace, index) {}
      
        U3State(
            const SubspaceType& subspace,
            const typename SubspaceType::StateLabelsType& state_labels
          )
          // Construct state by reverse lookup on labels.
          : basis::BaseDegenerateState<U3Subspace>(subspace,state_labels) {}

        // // pass-through accessors for subspace labels
        u3::U3 n() const {return labels().irrep;}
        int rho() const {return labels().tag;}
        
      // private:

    };

    class LGISpace
        : public basis::BaseSpace<U3Subspace,std::tuple<u3::U3>>
      {
        public:
        LGISpace() = default;
        LGISpace(const u3::U3& sigma, const int Nn_max);

        u3::U3 sigma() const {return std::get<0>(labels());}
        // private:

      };

    class Space
        : public basis::BaseSpace<LGISpace>
      {
        public:
          Space()=default;

          // Construct from list of sigma
          Space(
            const std::vector<u3::U3>& sigma_vector, 
            const int N0,
            const int Nmax
          );

          // Construct from spin::Space
          Space(const spin::Space& spin_space, const int Nmax);

      };
  
  }// namespace spatial

}  // namespace spncci
#endif