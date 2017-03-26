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
  ReadAndTransformSeedUnitTensorRMEs(
      const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const basis::MatrixVector& lgi_expansions,
      const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels,
      const std::string& filename,
      u3shell::SectorsU3SPN& unit_tensor_sectors,
      basis::MatrixVector& unit_tensor_spncci_matrices
    )
  {
    
    // lgi_unit_tensor_spncci_matrices.resize(lgi_unit_tensor_labels.size());
    // for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
    //   {
        // set up aliases for current unit tensor
        // const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];
        // const u3shell::SectorsU3SPN& unit_tensor_sectors = lgi_unit_tensor_sectors[unit_tensor_index];
        // const basis::MatrixVector& unit_tensor_lsu3shell_matrices = lgi_unit_tensor_lsu3shell_matrices[unit_tensor_index];
        // basis::MatrixVector& unit_tensor_spncci_matrices = lgi_unit_tensor_spncci_matrices[unit_tensor_index];
      
    const bool spin_scalar = false;

    // generate sector labels 
    unit_tensor_sectors = u3shell::SectorsU3SPN(lsu3shell_space,unit_tensor_labels,spin_scalar);
    
    // read in lsu3shell rms for unit tensor 
    basis::MatrixVector unit_tensor_lsu3shell_matrices;
    lsu3shell::ReadLSU3ShellRMEs(
        filename,
        lsu3shell_basis_table,lsu3shell_space,
        unit_tensor_labels,unit_tensor_sectors,unit_tensor_lsu3shell_matrices
      );

    // transform seed rmes to SpNCCI basis (among LGIs)
    lgi::TransformOperatorToSpBasis(
        unit_tensor_sectors,lgi_expansions,
        unit_tensor_lsu3shell_matrices,unit_tensor_spncci_matrices
      );
      // }
  }

  void
  StoreSeedUnitTensorRMEs(
      const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels,
      const u3shell::SectorsU3SPN& unit_tensor_sectors,
      const basis::MatrixVector& unit_tensor_spncci_matrices,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices,
      HalfInt Nsigma_max,
      double zero_threshold
    )
  {
    // for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
    //   {
        // set up aliases for current unit tensor
        // const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];
        // const u3shell::SectorsU3SPN& unit_tensor_sectors = lgi_unit_tensor_sectors[unit_tensor_index];
        // const basis::MatrixVector& unit_tensor_spncci_matrices = lgi_unit_tensor_spncci_matrices[unit_tensor_index];
      
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

            // // Restrict by lgi
            // if((bra_sigma.N()>Nsigma_max) || (ket_sigma.N()>Nsigma_max))
            //   continue;

            // put rme matrix into nested maps
            std::pair<int,int> irrep_family_index_pair(bra_subspace_index,ket_subspace_index);
            std::pair<int,int> Nn_pair(0,0);
            spncci::UnitTensorU3Sector unit_tensor_sector_labels(bra_sigma,ket_sigma,unit_tensor_labels,rho0);
            if(not mcutils::IsZero(unit_tensor_spncci_matrices[sector_index],zero_threshold))
              {
                unit_tensor_matrices[irrep_family_index_pair][Nn_pair][unit_tensor_sector_labels]
                  = unit_tensor_spncci_matrices[sector_index];
                // std::cout<<unit_tensor_spncci_matrices[sector_index]<<std::endl;
              }
          }
      // }
  }

  void
  RecurseUnitTensors(
      int N1v, int Nmax,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::KMatrixCache k_matrix_cache,
      u3::UCoefCache& u_coef_cache,
      u3::PhiCoefCache& phi_coef_cache,
      const std::map<int,std::vector<u3shell::RelativeUnitTensorLabelsU3ST>> unit_tensor_labels,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_matrices,
      bool verbose
    )
  {
    // verbosity control
    const int verbosity_interval = 100;
    int irrep_family_pair_count = 0;

    // diagnostics
    // TODO replace with a statstics.DebugStr()
    if (verbose)
      {
        spncci::UnitTensorMatrixStatistics statistics
          = spncci::GenerateUnitTensorMatrixStatistics(unit_tensor_matrices);
        spncci::WriteLog(
            fmt::format(
                "  Seed values before recursing LGI family pairs...\n"
                "    num_irrep_family_index_pairs    {}\n"
                "    num_Nn_pairs                    {}\n"
                "    num_unit_tensor_sectors         {}\n"
                "    num_nonzero_unit_tensor_sectors {} ({:.2e})\n"
                "    num_matrix_elements             {}\n"
                "    num_nonzero_matrix_elements     {} ({:.2e})\n",
                //irrep_family_pair_count,
                statistics.num_irrep_family_index_pairs,
                statistics.num_Nn_pairs,
                statistics.num_unit_tensor_sectors,
                statistics.num_nonzero_unit_tensor_sectors,
                double(statistics.num_nonzero_unit_tensor_sectors)/statistics.num_unit_tensor_sectors,
                statistics.num_matrix_elements,
                statistics.num_nonzero_matrix_elements,
                double(statistics.num_nonzero_matrix_elements)/statistics.num_matrix_elements
              )
          );
      }

    for(const auto& irrep_family_indices_submap_pair : unit_tensor_matrices)
      {
        std::pair<int,int> irrep_family_indices = irrep_family_indices_submap_pair.first;
        // std::cout<<"family pair "<<irrep_family_indices.first<<"  "<<irrep_family_indices.second<<std::endl;

        spncci::GenerateUnitTensorMatrix(
            N1v,Nmax,irrep_family_indices,spncci_space,u_coef_cache,phi_coef_cache,k_matrix_cache,
            unit_tensor_labels,unit_tensor_matrices);

        // diagnostics
        ++irrep_family_pair_count;
        if (verbose)
          {
            if ((irrep_family_pair_count%verbosity_interval==0)||(irrep_family_pair_count==unit_tensor_matrices.size()))
              {
                spncci::UnitTensorMatrixStatistics statistics
                  = spncci::GenerateUnitTensorMatrixStatistics(unit_tensor_matrices);

                // diagnostics
                spncci::WriteLog(
                    fmt::format(
                        "  After recursing {} LGI family pairs...\n"
                        "    num_irrep_family_index_pairs    {}\n"
                        "    num_Nn_pairs                    {}\n"
                        "    num_unit_tensor_sectors         {}\n"
                        "    num_nonzero_unit_tensor_sectors {} ({:.2e})\n"
                        "    num_matrix_elements             {}\n"
                        "    num_nonzero_matrix_elements     {} ({:.2e})\n",
                        irrep_family_pair_count,
                        statistics.num_irrep_family_index_pairs,
                        statistics.num_Nn_pairs,
                        statistics.num_unit_tensor_sectors,
                        statistics.num_nonzero_unit_tensor_sectors,
                        double(statistics.num_nonzero_unit_tensor_sectors)/statistics.num_unit_tensor_sectors,
                        statistics.num_matrix_elements,
                        statistics.num_nonzero_matrix_elements,
                        double(statistics.num_nonzero_matrix_elements)/statistics.num_matrix_elements
                      )
                  );
              }
          }
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

  void ConstructBranchedObservables(
    const spncci::SpaceU3S& space_u3s,
    const std::vector<std::vector<spncci::SectorLabelsU3S>>& observable_sectors_u3s,
    const std::vector<basis::MatrixVector>& observable_matrices_u3s,
    std::map<HalfInt,spncci::SpaceLS>& spaces_lsj,
    int num_observables,
    const std::vector<HalfInt>& J_values,
    int J0,
    std::vector<std::map<HalfInt,Eigen::MatrixXd>>& observable_matrices
    )
  {
    // populate fully-branched many-body matrices for observables
    // map: observable -> J ->  matrix
    // std::vector<std::map<HalfInt,Eigen::MatrixXd>> observable_matrices;  
    observable_matrices.resize(num_observables);
    for (int observable_index=0; observable_index<num_observables; ++observable_index)
      for (const HalfInt J : J_values)
        {
          // set up aliases (for current observable and J space)
          const std::vector<spncci::SectorLabelsU3S>& sectors_u3s = observable_sectors_u3s[observable_index];
          const basis::MatrixVector& matrices_u3s = observable_matrices_u3s[observable_index];
          const spncci::SpaceLS& space_lsj = spaces_lsj[J];
          Eigen::MatrixXd& observable_matrix = observable_matrices[observable_index][J];

          // determine set of (L0,S0) labels for this observable (triangular with J0)
          std::vector<spncci::OperatorLabelsLS> operator_labels_ls;
          // Note: to update when J0 varies by observable
          spncci::GenerateOperatorLabelsLS(J0,operator_labels_ls);

          // determine allowed LS sectors
          const spncci::SpaceLS& bra_space_lsj = space_lsj;
          const spncci::SpaceLS& ket_space_lsj = space_lsj;
          const HalfInt bra_J = J;
          const HalfInt ket_J = J;
          std::vector<spncci::SectorLabelsLS> sectors_lsj;
          spncci::GetSectorsLS(bra_space_lsj,ket_space_lsj,operator_labels_ls,sectors_lsj);

          // branch LS sectors to LSJ
          basis::MatrixVector matrices_lsj;  
          spncci::ContractAndRegroupLSJ(
              bra_J,J0,ket_J,
              space_u3s,sectors_u3s,matrices_u3s,
              bra_space_lsj,ket_space_lsj,sectors_lsj,matrices_lsj
            );

          // collect LSJ sectors into J matrix
          //
          // Note: Interface needs to be generalized to handle J_bra != J_ket.
          ConstructOperatorMatrix(
              space_lsj,sectors_lsj,matrices_lsj,
              observable_matrix
            );
        }
  }



}  // namespace
