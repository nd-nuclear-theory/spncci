/****************************************************************
  vcs_cache.h

  Caching of K matrices for spncci code.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  6/7/17 (mac): Extracted from spncci_basis.
  10/4/17 (aem): Add option for K matrices computed in basis with
    restricted sp3r->u3 branching (A<6)
****************************************************************/

#ifndef SPNCCI_VCS_CACHE_H_
#define SPNCCI_VCS_CACHE_H_

#include "sp3rlib/vcs.h"
#include "spncci/spncci_basis.h"

namespace spncci
{

  typedef std::unordered_map<u3::U3,vcs::MatrixCache,boost::hash<u3::U3>> KMatrixCache;
  // storage for K matrices
  //
  // maps sigma -> K matrix cache for that Sp irrep (vcs::MatrixCache),
  // where then vcs::MatrixCache maps omega to K matrix
  //
  // Usage: k_matrix_cache[sigma][omega]

  void
  PrecomputeKMatrices(
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::KMatrixCache& k_matrix_cache,
      spncci::KMatrixCache& kinv_matrix_cache,
      bool restrict_sp3r_u3_branching=false
    );
  // Precompute and cache K matrices for all symplectic irreps
  // occurring in SpNCCI space.
  //
  // May be used either in SpNCCI RME recurrence or in explicit construction of states.
  //
  // Arguments:
  //   sigma_irrep_map (input): container for distinct symplectic irreps
  //   k_matrix_cache (output): container for corresponding K matrices
  //   intrinsic (input): whether to compute K matrices for intrinsic (A-1)
  //     or lab (A) coordinates


}  // namespace

#endif
