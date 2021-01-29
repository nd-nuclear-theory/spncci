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
#include "mcutils/profiling.h"

namespace spncci
{

  void InitializeSpNCCI()
  {
    // SU(3) caching
    u3::U3CoefInit();
    u3::g_u_cache_enabled = true;

    // parameters for certain calculations
    spncci::g_zero_tolerance = 1e-6;
    spncci::g_suppress_zero_sectors = true;

    // Default binary mode, unless environment variable SPNCCI_RME_MODE
    // set to "text".
    lsu3shell::g_rme_binary_format = true;
    char* spncci_rme_mode_cstr = std::getenv("SPNCCI_RME_MODE");
    if (spncci_rme_mode_cstr!=NULL)
      {
        const std::string spncci_rme_mode = std::getenv("SPNCCI_RME_MODE");
        if (spncci_rme_mode=="text")
          lsu3shell::g_rme_binary_format = false;
      }
    
    // Eigen OpenMP multithreading mode
    Eigen::initParallel();
  
    //For testing  
    // std::cout<<"u3::g_u_cache_enabled="<<u3::g_u_cache_enabled<<", should be true"<<std::endl;
    // // u3::g_u_cache_enabled = true;
    // std::cout<<"spncci::g_zero_tolerance="<<spncci::g_zero_tolerance<<". should be 1e-6"<<std::endl;
    // // spncci::g_zero_tolerance = 1e-6;
    // std::cout<<"spncci::g_suppress_zero_sectors="<<spncci::g_suppress_zero_sectors<<". Should be true"<<std::endl;
    // // spncci::g_suppress_zero_sectors = true;

    // std::cout<<"lsu3shell::g_rme_binary_format="<<lsu3shell::g_rme_binary_format<<std::endl;
    // //Constructor for RunParameters opens file "spncci.dat" 
    // //and reads in run parameters storing them in run_parameters

  }


    void SetUpSpNCCISpaces(
      spncci::RunParameters& run_parameters,
      lgi::MultiplicityTaggedLGIVector& lgi_families,
      spncci::SpNCCISpace& spncci_space,
      spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::BabySpNCCISpace& baby_spncci_space,
      std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
      bool verbose
    )
  {
    std::cout << "Set up SpNCCI space..." << std::endl;

    // Get LGI families
    std::string lgi_filename="lgi_families.dat";
    lgi::ReadLGISet(lgi_filename, run_parameters.Nsigma0,lgi_families);

    // build SpNCCI irrep branchings
    spncci::NmaxTruncator truncator(run_parameters.Nsigma0,run_parameters.Nmax);

    // spncci::GenerateSpNCCISpace(lgi_families_truncated,truncator,spncci_space,sigma_irrep_map,restrict_sp3r_to_u3_branching);
    bool restrict_sp3r_to_u3_branching=false;
    if(run_parameters.A<6)
      restrict_sp3r_to_u3_branching=true;

    // spncci::SigmaIrrepMap sigma_irrep_map;
    spncci::GenerateSpNCCISpace(lgi_families,truncator,spncci_space,sigma_irrep_map,restrict_sp3r_to_u3_branching);

    if(verbose)
      {
        for(int i=0; i<spncci_space.size(); ++i)
          std::cout<<i<<"  "<<spncci_space[i].Str()<<spncci_space[i].gamma_max()<<std::endl;

        // diagnostics
        std::cout << fmt::format("  Irrep families {}",spncci_space.size()) << std::endl;
        std::cout << fmt::format("  TotalU3Subspaces {}",spncci::TotalU3Subspaces(spncci_space)) << std::endl;
        std::cout << fmt::format("  TotalDimensionU3S {}",spncci::TotalDimensionU3S(spncci_space)) << std::endl;
      }

    // build baby spncci space 
    baby_spncci_space=spncci::BabySpNCCISpace(spncci_space);

    //Full basis by J values.  This is the space actually used in spncci
    spaces_spbasis.resize(run_parameters.J_values.size());
    for(int j=0; j<run_parameters.J_values.size(); ++j)
      {
        const HalfInt& J=run_parameters.J_values[j];
        spaces_spbasis[j]=spncci::SpaceSpBasis(baby_spncci_space,J);
      }

  }


  void SetUpSpNCCISpaces(
      spncci::RunParameters& run_parameters,
      lgi::MultiplicityTaggedLGIVector& lgi_families,
      spncci::SpNCCISpace& spncci_space,
      spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::BabySpNCCISpace& baby_spncci_space,
      std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
      spncci::KMatrixCache& k_matrix_cache, 
      spncci::KMatrixCache& kinv_matrix_cache,
      bool verbose
    )
    {
      // Set up SpNCCI Spaces 
      // spncci::SigmaIrrepMap sigma_irrep_map;
      spncci::SetUpSpNCCISpaces(run_parameters,lgi_families,spncci_space,sigma_irrep_map,baby_spncci_space,spaces_spbasis,verbose);
        ////////////////////////////////////////////////////////////////
      // precompute K matrices
      ////////////////////////////////////////////////////////////////
      std::cout << "Precompute K matrices..." << std::endl;
      bool restrict_sp3r_to_u3_branching=false;
        if(run_parameters.A<6)
          restrict_sp3r_to_u3_branching=true;
      // timing start
      mcutils::SteadyTimer timer_k_matrices;
      timer_k_matrices.Start();

      // traverse distinct sigma values in SpNCCI space, generating K
      spncci::PrecomputeKMatrices(sigma_irrep_map,k_matrix_cache,kinv_matrix_cache,restrict_sp3r_to_u3_branching);

      // timing stop
      timer_k_matrices.Stop();
      std::cout << fmt::format("(Task time: {})",timer_k_matrices.ElapsedTime()) << std::endl;

      std::cout<<"Kmatrices "<<std::endl;
      for(auto it=k_matrix_cache.begin(); it!=k_matrix_cache.end(); ++it)
        {
          std::cout<<"sigma "<<it->first.Str()<<std::endl;
          for(auto it2=it->second.begin();  it2!=it->second.end(); ++it2)
          {
            std::cout<<"  omega"<<it2->first.Str()<<std::endl;
            auto matrix=it2->second;
            std::cout<<matrix<<std::endl;
            // std::cout<<matrix.inverse()<<std::endl;
          }
        }
    }

}  // namespace
