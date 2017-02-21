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

  void
  RecurseUnitTensors(
      int N1v, int Nmax,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::KMatrixCache k_matrix_cache,
      u3::UCoefCache& u_coef_cache,
      u3::PhiCoefCache& phi_coef_cache,
      const std::map<int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>> unit_tensor_labels,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices
    )
  {
    for(const auto& irrep_family_indices_submap_pair : unit_tensor_matrices)
      {
        std::pair<int,int> irrep_family_indices = irrep_family_indices_submap_pair.first;
        spncci::GenerateUnitTensorMatrix(
            N1v,Nmax,irrep_family_indices,spncci_space,u_coef_cache,phi_coef_cache,k_matrix_cache,
            unit_tensor_labels,unit_tensor_matrices);
      }
  }


  void
  ConstructObservableU3S(
      int Nmax, int N1v,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceU3S& space_u3s,
      const spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices,
      const u3shell::RelativeRMEsU3ST& relative_rmes,
      std::vector<spncci::SectorLabelsU3S>& sectors_u3s,
      basis::MatrixVector& matrices_u3s
    )
  {
    // identify operator symmetry and upcoupling labels
    //
    // Note: see branching_u3s_test for possible scheme for
    // reimplementation
    std::vector<u3shell::IndexedOperatorLabelsU3S> operator_u3s_list;
    u3shell::GetInteractionTensorsU3S(relative_rmes,operator_u3s_list);

    // constract and regroup at U3S level
    spncci::GetSectorsU3S(space_u3s,operator_u3s_list,sectors_u3s);
    spncci::ContractAndRegroupU3S(
        Nmax,N1v,sectors_u3s,relative_rmes,baby_spncci_space,
        space_u3s,unit_tensor_matrices,matrices_u3s
      );
  }

}  // namespace
