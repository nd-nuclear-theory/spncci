/****************************************************************
  results_output.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/results_output.h"

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
    WriteKeyValue(out_stream,"A",":d",run_parameters.A);
    WriteKeyValue(out_stream,"Nsigma0",":.1f",float(run_parameters.Nsigma0));
    WriteKeyValue(out_stream,"Nsigmamax",":d",run_parameters.Nsigmamax);
    WriteKeyValue(out_stream,"N1v",":d",run_parameters.N1v);
    WriteKeyValue(out_stream,"Nmax",":d",run_parameters.Nmax);

    StartNewSection(out_stream,"Upstream");
    WriteKeyValueList(out_stream,"nuclide",":d",run_parameters.nuclide);
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
      const spncci::SpaceSpLS spls_space
    )
  {
    StartNewSection(out_stream,"Irreps");

    StartNewSection(out_stream,"BabySpNCCI");

    StartNewSection(out_stream,"SpU3S");
    WriteKeyValue(out_stream,"subspaces",":d",spu3s_space.size());
    WriteKeyValue(out_stream,"dimension",":d",spu3s_space.Dimension());
    WriteKeyValue(out_stream,"full_dimension",":d",spu3s_space.FullDimension());

    StartNewSection(out_stream,"SpLS");

  }

  void WriteCalculationParameters(std::ostream& out_stream, double hw)
  {
    StartNewSection(out_stream,"Calculation");
    WriteKeyValue(out_stream,"hw",":.3f",hw);
  }

  void WriteEigenvalues(
      std::ostream& out_stream,
      const std::map<HalfInt,Eigen::VectorXd>& eigenvalues
    )
  {
    StartNewSection(out_stream,"Energies");
    out_stream << "# J gex i E" << std::endl;
    for (const auto& J_eigenvalues_J : eigenvalues)
      {
        // extract results for J subspace
        HalfInt J = J_eigenvalues_J.first;
        const Eigen::VectorXd& eigenvalues_J = J_eigenvalues_J.second;
        int gex = 0;  // CAVEAT: placeholder value

        // iterate over eigenvalues in J subspace
        for (int eigenstate_index=0; eigenstate_index<eigenvalues_J.size(); ++eigenstate_index)
          out_stream
            << fmt::format("{:4.1f} {:1d} {:3d} {:+9.4f}",double(J),gex,eigenstate_index,eigenvalues_J[eigenstate_index])
            << std::endl;
      }
  }


}  // namespace
