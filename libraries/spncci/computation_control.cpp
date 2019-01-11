/****************************************************************
  computation_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/computation_control.h"

// #include "SymEigsSolver.h"  // from spectra
#include "fmt/format.h"
// #include "lgi/lgi_solver.h"
// #include "mcutils/eigen.h"
#include "spncci/results_output.h"


namespace spncci
{
  void SetUpSpNCCISpaces(
      spncci::RunParameters& run_parameters,
      lgi::MultiplicityTaggedLGIVector& lgi_families,
      spncci::SpNCCISpace& spncci_space,
      spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::BabySpNCCISpace& baby_spncci_space,
      spncci::SpaceSpU3S& spu3s_space,
      spncci::SpaceSpLS& spls_space,
      spncci::SpaceSpJ& spj_space,
      std::ofstream& results_stream,
      int Nlimit,
      bool restrict_sp3r_to_u3_branching
    )
  {
    // Get LGI families
    std::string lgi_filename="lgi_families.dat";
    lgi::ReadLGISet(lgi_filename, run_parameters.Nsigma0,lgi_families);

    std::cout << "Set up SpNCCI space..." << std::endl;

    // build SpNCCI irrep branchings
    spncci::NmaxTruncator truncator(run_parameters.Nsigma0,run_parameters.Nmax);
    // spncci::NlimitTruncator truncator(run_parameters.Nsigma0,run_parameters.Nmax,Nlimit);

    // spncci::GenerateSpNCCISpace(lgi_families_truncated,truncator,spncci_space,sigma_irrep_map,restrict_sp3r_to_u3_branching);
    spncci::GenerateSpNCCISpace(lgi_families,truncator,spncci_space,sigma_irrep_map,restrict_sp3r_to_u3_branching);

    for(int i=0; i<spncci_space.size(); ++i)
      std::cout<<i<<"  "<<spncci_space[i].Str()<<spncci_space[i].gamma_max()<<std::endl;

    // diagnostics
    std::cout << fmt::format("  Irrep families {}",spncci_space.size()) << std::endl;
    std::cout << fmt::format("  TotalU3Subspaces {}",spncci::TotalU3Subspaces(spncci_space)) << std::endl;
    std::cout << fmt::format("  TotalDimensionU3S {}",spncci::TotalDimensionU3S(spncci_space)) << std::endl;

    // build baby spncci space 
    baby_spncci_space=spncci::BabySpNCCISpace(spncci_space);

    // build SpU3S gathered space
    std::cout << "Build SpU3S space..." << std::endl;
    spu3s_space=spncci::SpaceSpU3S(baby_spncci_space);
    std::cout
      << fmt::format("  subspaces {} dimension {} full_dimension {}",
                     spu3s_space.size(),spu3s_space.Dimension(),spu3s_space.FullDimension()
        )
      << std::endl;
    std::cout
      << fmt::format("  compare... TotalDimensionU3S {}",
                     TotalDimensionU3S(spncci_space)
        )
      << std::endl;
    // std::cout << spu3s_space.DebugStr(true);

	  // build SpLS branched space
	  std::cout << "Build SpLS space..." << std::endl;
	  spls_space=spncci::SpaceSpLS(spu3s_space);
	  std::cout
	    << fmt::format("  subspaces {} dimension {} full_dimension {}",
	                   spls_space.size(),spls_space.Dimension(),spls_space.FullDimension()
	      )
	    << std::endl;
	  std::cout
	    << fmt::format("  compare... TotalDimensionU3LS {}",TotalDimensionU3LS(spncci_space))
	    << std::endl;
	  // std::cout << splss_space.DebugStr(true);

	  // build SpJ branched space
	  std::cout << "Build SpJ space..." << std::endl;
	  spj_space=spncci::SpaceSpJ(run_parameters.J_values,spls_space);
	  std::cout
	    << spj_space.DebugStr(false)
	    << std::endl;
	  std::cout
	    << fmt::format("  subspaces {}",spj_space.size())
	    << std::endl;


	  // results output: basis information
	  spncci::StartNewSection(results_stream,"BASIS");
	  spncci::WriteBasisStatistics(results_stream,spncci_space,baby_spncci_space,spu3s_space,spls_space,spj_space);
	  spncci::WriteSpU3SSubspaceListing(results_stream,baby_spncci_space,run_parameters.Nsigma0);
	  spncci::WriteBabySpNCCISubspaceListing(results_stream,baby_spncci_space,run_parameters.Nsigma0);

	}

}  // namespace
