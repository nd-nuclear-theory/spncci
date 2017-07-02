/****************************************************************
  vcs_cache.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/vcs_cache.h"

namespace spncci
{
  void
  PrecomputeKMatrices(
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::KMatrixCache& k_matrix_cache
    )
  {
    for (const auto& sigma_irrep_pair : sigma_irrep_map)
      {
        // extract sigma and irrep contents
        const u3::U3& sigma = sigma_irrep_pair.first;
        const sp3r::Sp3RSpace& sp_irrep = sigma_irrep_pair.second;

        // populate K matrix cache for this irrep
        vcs::GenerateKMatrices(sp_irrep,k_matrix_cache[sigma]);
      }

  }

}  // namespace
