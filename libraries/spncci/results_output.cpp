/****************************************************************
  results_output.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/results_output.h"

#include "mcutils/eigen.h"
#include "spncci/parameters.h"
#include "spncci/spncci_common.h"


namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // output utilities
  ////////////////////////////////////////////////////////////////

  void StartNewSection(std::ostream& out_stream, const std::string& title)
  {
    out_stream << std::endl;
    out_stream << fmt::format("[{}]",title) << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // output code
  ////////////////////////////////////////////////////////////////

  void WriteCodeInformation(std::ostream& out_stream, const spncci::RunParameters& run_parameters)
  {
    StartNewSection(out_stream,"Version");
    WriteKeyValue(out_stream,"results",":d",g_results_version);

    StartNewSection(out_stream,"Internals");
    WriteKeyValue(out_stream,"g_zero_tolerance",":e",g_zero_tolerance);

    // StartNewSection(out_stream,"Parallelization");

  }

  void WriteRunParameters(std::ostream& out_stream, const spncci::RunParameters& run_parameters)
  {
    StartNewSection(out_stream,"Space");
    WriteKeyValueList(out_stream,"nuclide",":d",run_parameters.nuclide);
    WriteKeyValue(out_stream,"A",":d",run_parameters.A);
    WriteKeyValue(out_stream,"Nsigma0",":.1f",float(run_parameters.Nsigma0));
    WriteKeyValue(out_stream,"Nsigmamax",":d",run_parameters.Nsigmamax);
    WriteKeyValue(out_stream,"N1v",":d",run_parameters.N1v);
    WriteKeyValue(out_stream,"Nmax",":d",run_parameters.Nmax);

    StartNewSection(out_stream,"Interaction");
    WriteKeyValue(out_stream,"interaction",":s",run_parameters.interaction_name);
    WriteKeyValue(out_stream,"use_coulomb",":d",run_parameters.use_coulomb);

    StartNewSection(out_stream,"Eigensolver");
    WriteKeyValue(out_stream,"num_eigenvalues",":d",run_parameters.num_eigenvalues);
    WriteKeyValue(out_stream,"num_convergence",":d",run_parameters.eigensolver_num_convergence);
    WriteKeyValue(out_stream,"max_iterations",":d",run_parameters.eigensolver_max_iterations);
    WriteKeyValue(out_stream,"tolerance",":e",run_parameters.eigensolver_tolerance);

    StartNewSection(out_stream,"Mesh");
    WriteKeyValueList(out_stream,"hw",":.3f",run_parameters.hw_values);

    StartNewSection(out_stream,"Branching");
    std::vector<double> J_values_double;
    for (const HalfInt J : run_parameters.J_values)
      J_values_double.push_back(double(J));  // need to convert HalfInt's to double for output by WriteKeyValueList
    WriteKeyValueList(out_stream,"J",":.1f",J_values_double);

    StartNewSection(out_stream,"Observables");
    WriteKeyValueList(out_stream,"filenames",":s",run_parameters.observable_filenames);

  }

  void WriteBasisStatistics(
      std::ostream& out_stream,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceSpU3S& spu3s_space,
      const spncci::SpaceSpLS& spls_space,
      const spncci::SpaceSpJ& spj_space
    )
  {
    // SpNCCISpace
    StartNewSection(out_stream,"SpNCCI");
    WriteKeyValue(out_stream,"irrep_families",":d",spncci_space.size());

    // BabySpNCCI
    StartNewSection(out_stream,"BabySpNCCI");
    WriteKeyValue(out_stream,"subspaces",":d",baby_spncci_space.size());
    // WriteKeyValue(out_stream,"dimension",":d",baby_spncci_space.Dimension());
    // WriteKeyValue(out_stream,"full_dimension",":d",baby_spncci_space.FullDimension());

    // SpU3S
    StartNewSection(out_stream,"SpU3S");
    WriteKeyValue(out_stream,"subspaces",":d",spu3s_space.size());
    WriteKeyValue(out_stream,"dimension",":d",spu3s_space.Dimension());
    WriteKeyValue(out_stream,"full_dimension",":d",spu3s_space.FullDimension());

    // SpLS
    StartNewSection(out_stream,"SpLS");
    WriteKeyValue(out_stream,"subspaces",":d",spls_space.size());
    WriteKeyValue(out_stream,"dimension",":d",spls_space.Dimension());
    WriteKeyValue(out_stream,"full_dimension",":d",spls_space.FullDimension());

    // SpJ
    StartNewSection(out_stream,"SpJ");
    WriteKeyValue(out_stream,"subspaces",":d",spj_space.size());
    WriteKeyValue(out_stream,"dimension",":d",spj_space.Dimension());
    WriteKeyValue(out_stream,"full_dimension",":d",spj_space.FullDimension());

    // SpJ
    StartNewSection(out_stream,"SpJ (listing)");
    out_stream << "# subspace_index J dim" << std::endl;
    for (int subspace_index=0; subspace_index<spj_space.size(); ++subspace_index)
      {
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(subspace_index);
        out_stream
          << fmt::format(
              "{:3d} {:4.1f} {:10d}",
              subspace_index,double(spj_subspace.J()),spj_subspace.full_dimension()
            )
          << std::endl;
      }

  }

  void WriteSpU3SSubspaceListing(
      std::ostream& out_stream,
      const spncci::SpaceSpU3S& spu3s_space,
      HalfInt Nsigma0
    )
  {

    // SpU3S
    StartNewSection(out_stream,"SpU3S (listing)");
    out_stream
      << "# subspace_index " << std::endl
      << "# Nex omega.N omega.lambda omega.mu" << std::endl
      << "# S " << std::endl
      << "# dim" << std::endl;
    for (int subspace_index=0; subspace_index<spu3s_space.size(); ++subspace_index)
      {
        const SubspaceSpU3S& spu3s_subspace = spu3s_space.GetSubspace(subspace_index);
        const u3::U3 omega = spu3s_subspace.omega();
        int Nex = int(spu3s_subspace.omega().N()-Nsigma0);
        out_stream
          << fmt::format(
              "{:5d}   "
              "{:2d} {:5.1f} {:3d} {:3d}   "
              "{:5.1f}   "
              "{:5d}",
              subspace_index,
              Nex,float(omega.N()),omega.SU3().lambda(),omega.SU3().mu(),
              float(spu3s_subspace.S()),
              spu3s_subspace.full_dimension()
            )
          << std::endl;
      }


  }

  void WriteBabySpNCCISubspaceListing(
      std::ostream& out_stream,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      HalfInt Nsigma0
    )
  {

    // BabySpNCCI
    StartNewSection(out_stream,"BabySpNCCI (listing)");
    out_stream
      << "# subspace_index irrep_family_index" << std::endl
      << "# Nsigmaex sigma.N sigma.lambda sigma.mu" << std::endl
      << "# Sp Sn S " << std::endl
      << "# Nex omega.N omega.lambda omega.mu" << std::endl
      << "# gamma_max upsilon_max dim" << std::endl;
    for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
      {
        const BabySpNCCISubspace& baby_spncci_subspace = baby_spncci_space.GetSubspace(subspace_index);
        const u3::U3 sigma = baby_spncci_subspace.sigma();
        const u3::U3 omega = baby_spncci_subspace.omega();
        int Nsigmaex = int(baby_spncci_subspace.sigma().N()-Nsigma0);
        int Nex = int(baby_spncci_subspace.omega().N()-Nsigma0);
        out_stream
          << fmt::format(
              "{:5d} {:5d}   "
              "{:2d} {:5.1f} {:3d} {:3d}   "
              "{:5.1f} {:5.1f} {:5.1f}   "
              "{:2d} {:5.1f} {:3d} {:3d}   "
              "{:3d} {:3d} {:4d}",
              subspace_index,baby_spncci_subspace.irrep_family_index(),
              Nsigmaex,float(sigma.N()),sigma.SU3().lambda(),sigma.SU3().mu(),
              float(baby_spncci_subspace.Sp()),float(baby_spncci_subspace.Sn()),float(baby_spncci_subspace.S()),
              Nex,float(omega.N()),omega.SU3().lambda(),omega.SU3().mu(),
              baby_spncci_subspace.gamma_max(),baby_spncci_subspace.upsilon_max(),baby_spncci_subspace.size()
            )
          << std::endl;
      }

  }


  // void WriteU3SSectorInformation(
  //     std::ostream& out_stream,
  //     const spncci::SpaceU3S& space_u3s,
  //     int num_observables, 
  //     const std::vector<std::vector<spncci::SectorLabelsU3S>>& observables_sectors_u3s
  //   )
  // {

  //   StartNewSection(out_stream,"U3S sector dimensions");
  //   out_stream
  //     << "# observable_index num_sectors max_sector_entries observable_entries"
  //     << std::endl;

  //   // for each observable, enumerate sectors 
  //   int total_entries = 0;  // total across observables, but not across hw
  //   int max_sector_entries = 0;
  //   int max_observable_entries = 0;
  //   for(int observable_index=0; observable_index<num_observables; ++observable_index) 
  //     {
  //       const std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];

  //       // do counting for this observable
  //       int observable_entries = 0;
  //       for(int sector_index=0; sector_index<sectors_u3s.size(); ++sector_index)
  //         {
  //           int rows=space_u3s.GetSubspace(sectors_u3s[sector_index].bra_index()).full_dimension();
  //           int cols=space_u3s.GetSubspace(sectors_u3s[sector_index].ket_index()).full_dimension();
  //           int sector_entries = rows*cols;
  //           max_sector_entries = std::max(max_sector_entries,sector_entries);
  //           observable_entries += sector_entries;
  //         }

  //       // write statistics for this observable
  //       out_stream
  //         << fmt::format(
  //             "{:2d} {:5d} {:10d} {:10d}",
  //             observable_index,sectors_u3s.size(),max_sector_entries,observable_entries
  //           )
  //         << std::endl;


  //       // accumulate counting from this observable -- totals across observables are not currently output
  //       max_observable_entries = std::max(max_observable_entries,observable_entries);
  //       total_entries += observable_entries;
  //     }
  // }

    void WriteU3SHypersectorSectorInformation(
      std::ostream& out_stream,
      const spncci::SpaceU3S& space_u3s,
      int num_observables, 
      const std::vector<spncci::ObservableHypersectorsU3S>& observables_hypersectors_u3s
    )
  {

    StartNewSection(out_stream,"U3S sector dimensions");
    out_stream
      << "# observable_index num_sectors max_sector_entries observable_entries"
      << std::endl;

    // for each observable, enumerate sectors 
    int total_entries = 0;  // total across observables, but not across hw
    int max_sector_entries = 0;
    long int max_observable_entries = 0;
    for(int observable_index=0; observable_index<num_observables; ++observable_index) 
      {
        const spncci::ObservableHypersectorsU3S& hypersectors_u3s
                =observables_hypersectors_u3s[observable_index];

        // do counting for this observable
        long int observable_entries = 0;
        for(int sector_index=0; sector_index<hypersectors_u3s.size(); ++sector_index)
          {
            const auto& hypersector_u3s=hypersectors_u3s.GetHypersector(sector_index);
            int rows=space_u3s.GetSubspace(hypersector_u3s.bra_subspace_index()).full_dimension();
            int cols=space_u3s.GetSubspace(hypersector_u3s.ket_subspace_index()).full_dimension();
            int sector_entries = rows*cols;
            max_sector_entries = std::max(max_sector_entries,sector_entries);
            observable_entries += sector_entries;
          }

        // write statistics for this observable
        out_stream
          << fmt::format(
              "{:2d} {:5d} {:10d} {:10d}",
              observable_index,hypersectors_u3s.size(),max_sector_entries,observable_entries
            )
          << std::endl;


        // accumulate counting from this observable -- totals across observables are not currently output
        max_observable_entries = std::max(max_observable_entries,observable_entries);
        total_entries += observable_entries;
      }
  }



  void WriteCalculationParameters(std::ostream& out_stream, double hw)
  {
    StartNewSection(out_stream,"Calculation");
    WriteKeyValue(out_stream,"hw",":.3f",hw);
  }

  void WriteEigenvalues(
      std::ostream& out_stream,
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Vector>& eigenvalues,
      int gex
    )
  {
    StartNewSection(out_stream,"Energies");
    out_stream << "# J gex i E" << std::endl;
    for (int subspace_index=0; subspace_index<spj_space.size(); ++subspace_index)
      {
        HalfInt J = spj_space.GetSubspace(subspace_index).J();
        const Eigen::VectorXd& eigenvalues_J = eigenvalues[subspace_index];

        // iterate over eigenvalues in J subspace
        for (int eigenstate_index=0; eigenstate_index<eigenvalues_J.size(); ++eigenstate_index)
          out_stream
            << fmt::format("{:4.1f} {:1d} {:3d} {:+9.4f}",double(J),gex,eigenstate_index,eigenvalues_J[eigenstate_index])
            << std::endl;
      }
  }

  void WriteDecompositions(
      std::ostream& out_stream,
      const std::string& decomposition_name,
      const std::string& format_string,
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Matrix>& decompositions,
      int gex
    )
  {
    StartNewSection(out_stream,fmt::format("Decompositions: {}",decomposition_name));

    for (int subspace_index=0; subspace_index<spj_space.size(); ++subspace_index)
      {
        // retrieve information for subspace
        HalfInt J = spj_space.GetSubspace(subspace_index).J();
        const spncci::Matrix& decompositions_J = decompositions[subspace_index];

        // short circuit empty subspace
        if (decompositions_J.cols()==0)
          continue;

        // write header comment for subspace
        out_stream
          << fmt::format(
              "# decompositions for subspace J={:.1f}, gex={:1d} ({:d}x{:d})",
              float(J),gex,decompositions_J.rows(),decompositions_J.cols()
            )
          << std::endl;

        // write decompositions
        out_stream << mcutils::FormatMatrix(decompositions_J,format_string) << std::endl;
      }
  }

  void WriteObservables(
      std::ostream& out_stream,
      const std::vector<spncci::SectorsSpJ>& observable_sectors,
      const std::vector<spncci::OperatorBlocks>& observable_results_matrices,
      int gex
    )
  {
    StartNewSection(out_stream,"Observables");
    out_stream << "# observable_index sector_index J_bra gex_bra J_ket gex_ket rows cols" << std::endl;
    for (int observable_index=0; observable_index<observable_results_matrices.size(); ++observable_index)
      {
        
        // retrieve sectors
        const spncci::SectorsSpJ& sectors = observable_sectors[observable_index];

        // tabulate observable on each sector
        for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
          {
            
            // retrieve sector information
            const spncci::SectorsSpJ::SectorType& sector = sectors.GetSector(sector_index);
            const HalfInt bra_J = sector.bra_subspace().J();
            const HalfInt ket_J = sector.ket_subspace().J();

            // retrieve block
            const spncci::OperatorBlock& observable_results_matrix = observable_results_matrices[observable_index][sector_index];

            // short circuit empty block
            if ((observable_results_matrix.rows()==0)||(observable_results_matrix.cols()==0))
              continue;

            out_stream
              << fmt::format(
                "{:d} {:d} {:.1f} {:d} {:.1f} {:d} {:d} {:d} ",
                observable_index,sector_index,double(bra_J),gex,double(ket_J),gex,
                observable_results_matrix.rows(),observable_results_matrix.cols()
                )
              << std::endl;
            out_stream
              << mcutils::FormatMatrix(observable_results_matrix,"13.6e")
              << std::endl;

          }
      }
  }


}  // namespace
