/****************************************************************
  computation_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/computation_control.h"

#include "cppformat/format.h"
#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // processing seed RMEs
  ////////////////////////////////////////////////////////////////

  void
  TransformSeedUnitTensorRMEs(
      const basis::MatrixVector& lgi_expansions,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::vector<u3shell::SectorsU3SPN>& lgi_unit_tensor_sectors,
      const std::vector<basis::MatrixVector>& lgi_unit_tensor_lsu3shell_matrices,
      std::vector<basis::MatrixVector>& lgi_unit_tensor_spncci_matrices
    )
  {
    lgi_unit_tensor_spncci_matrices.resize(lgi_unit_tensor_labels.size());
    for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
      {
        // set up aliases for current unit tensor
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];
        const u3shell::SectorsU3SPN& unit_tensor_sectors = lgi_unit_tensor_sectors[unit_tensor_index];
        const basis::MatrixVector& unit_tensor_lsu3shell_matrices = lgi_unit_tensor_lsu3shell_matrices[unit_tensor_index];
        basis::MatrixVector& unit_tensor_spncci_matrices = lgi_unit_tensor_spncci_matrices[unit_tensor_index];
      
        // transform seed rmes to SpNCCI basis (among LGIs)
        lgi::TransformOperatorToSpBasis(
            unit_tensor_sectors,lgi_expansions,
            unit_tensor_lsu3shell_matrices,unit_tensor_spncci_matrices
          );
      }
  }

  void
  StoreSeedUnitTensorRMEs(
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::vector<u3shell::SectorsU3SPN>& lgi_unit_tensor_sectors,
      const std::vector<basis::MatrixVector>& lgi_unit_tensor_spncci_matrices,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices,
      double zero_threshold
    )
  {
    for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
      {
        // set up aliases for current unit tensor
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];
        const u3shell::SectorsU3SPN& unit_tensor_sectors = lgi_unit_tensor_sectors[unit_tensor_index];
        const basis::MatrixVector& unit_tensor_spncci_matrices = lgi_unit_tensor_spncci_matrices[unit_tensor_index];
      
        // stash each sector in big souffle (i.e., by irrep family)
        for(int sector_index=0; sector_index<unit_tensor_sectors.size(); ++sector_index)
          {
            // extract U3SPN sector information
            const typename u3shell::SectorsU3SPN::SectorType& sector = unit_tensor_sectors.GetSector(sector_index);
            const int bra_subspace_index = sector.bra_subspace_index();
            const int ket_subspace_index = sector.ket_subspace_index();
            const u3::U3& bra_sigma = sector.bra_subspace().U3();
            const u3::U3& ket_sigma = sector.ket_subspace().U3();
            const int rho0 = unit_tensor_sectors.GetSector(sector_index).multiplicity_index();

            // put rme matrix into nested maps
            std::pair<int,int> irrep_family_index_pair(bra_subspace_index,ket_subspace_index);
            std::pair<int,int> Nn_pair(0,0);
            spncci::UnitTensorU3Sector unit_tensor_sector_labels(bra_sigma,ket_sigma,unit_tensor_labels,rho0);
            if(not mcutils::IsZero(unit_tensor_spncci_matrices[sector_index],zero_threshold))
              {
                unit_tensor_matrices[irrep_family_index_pair][Nn_pair][unit_tensor_sector_labels]
                  = unit_tensor_spncci_matrices[sector_index];
              }
          }
      }
  }

}  // namespace
