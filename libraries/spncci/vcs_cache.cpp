/****************************************************************
  vcs_cache.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/vcs_cache.h"

namespace spncci
{
  void
  PrecomputeKMatrices(
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::KMatrixCache& k_matrix_cache,
      spncci::KMatrixCache& kinv_matrix_cache,
      bool restrict_sp3r_u3_branching
    )
  {
    for (const auto& sigma_irrep_pair : sigma_irrep_map)
      {
        // extract sigma and irrep contents
        const u3::U3& sigma = sigma_irrep_pair.first;
        const sp3r::Sp3RSpace& sp_irrep = sigma_irrep_pair.second;

        // populate K matrix cache for this irrep
        if(restrict_sp3r_u3_branching)
          vcs::GenerateKMatrices(sp_irrep,k_matrix_cache[sigma],kinv_matrix_cache[sigma]);
        else
          {
            auto& k_matrices=k_matrix_cache[sigma];
            auto& kinv_matrices=kinv_matrix_cache[sigma];
            vcs::GenerateKMatrices(sp_irrep,k_matrices);
            for(auto it=k_matrix_cache[sigma].begin(); it!=k_matrix_cache[sigma].end(); ++it)
              {
                kinv_matrices[it->first]=it->second.inverse();
              }
          }
      }

  }

}  // namespace
