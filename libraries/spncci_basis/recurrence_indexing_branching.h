/****************************************************************
recurrence_indexing_spatial.h

    Indexing for SpNCCI recurrence

    spncci::spatial::Space() []
    ->spncci::spatial::Sp3RSpace() [sigma]
        ->spncci::spatial::U3Subspace() [omega] (upsilon_max)
          ->spncci::spatial::U3State() [L] (kappa_max)

    spatial::RecurrenceSpace() []
    ->spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
      ->spatial::RecurrenceNnsumSpace() [Nnsum]
        ->spatial::RecurrenceU3Space() [omega,omega'] (upsilon,upsilon')
          ->spatial::RecurrenceOperatorSubspace() [x0] (rho0)
            ->spatial::RecurrenceOperatorState() [Nbar,Nbar']

    spatial::ContractionSpace() [J0]
    ->spatial::ContractionSp3RSpace() [sigma,sigma',parity_bar]
      ->spatial::ContractionU3Space() [omega,omega'] (upsilon,upsilon')
        ->spatial::ContractionOperatorSubspace() [L0] (kappa0) <-- J0-2<=L0<=J0+2
          ->spatial::ContractionOperatorState() [x0] (rho0)

    spatial::BranchingSpace() [J,J',J0]
    ->spatial::BranchingSp3RSpace() [sigma,sigma',parity_bar]
      ->spatial::BranchingU3Subspace() [omega,omega'] (upsilon,upsilon')
        ->spatial::BranchingState() [L,L'] (kappa,kappa')
    TODO: Find efficient way to actually store in (L,kappa,L',kappa') order


    Final indexing of branched J space
    -> J
      -> sigma S Sp Sn gamma omega upsilon L kappa

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  + 05/10/22 (aem): Created.
****************************************************************/
#ifndef RECURRENCE_INDEXING_BRANCHING_H_
#define RECURRENCE_INDEXING_BRANCHING_H_

#include <cppitertools/enumerate.hpp>  // iter
#include <memory>                      // for std::shared_ptr
#include <tuple>
#include <utility>  // for std::forward

#include "am/halfint.h"
#include "basis/basis.h"
#include "basis/degenerate.h"
#include "fmt/format.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/u3.h"
#include "spncci_basis/recurrence_indexing_spatial.h"
#include "spncci_basis/recurrence_indexing_spin.h"
#include "u3shell/operator_indexing_spatial.h"

namespace spncci
{

// Starting from output of recurrence sector with seeds
// Branch.

// Sectors between two U3Subspaces
//using U3Subspace = sp3r::U3Subspace;
//Need Sp3RSpace restricted by specific S and J.

// Construct Operator block for omega,omega' block with fixed spin quantum numbers.
// Block will be by upsilon, by L, by kappa.

}


#endif
