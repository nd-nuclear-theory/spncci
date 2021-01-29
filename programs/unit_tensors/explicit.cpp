/****************************************************************
  explicit.cpp

  Tests of explicit SpNCCI basis construction in LSU3Shell basis.

  This code just tests normalization, but using clean refactored
  infrastructure.  Other deeper tests (of unit tensor matrix elements)
  were carried out in compute_unit_tensor_rmes.cpp.

  CAVEAT: Right now run parameters are hard coded in RunParameters
  constructor.

  Required data:

    Input files are generated by

      generate_lsu3shell_relative_operators.cpp

    which is invoked through scripting in

      compute_relative_tensors_lsu3shell_rmes.py

    Example: Z=3 N=3 twice_Nsigma0=22 Nmax=2 Nstep=2 N1v[=N1b]=1
   
    % python3 script/compute_relative_tensors_lsu3shell_rmes.py 3 3 22 2 2 1

    Only need .rme and .dat files.

    Not saved to repository since ~3.5 Mb...

       data/lsu3shell/lsu3shell_rme_6Li_Nmax02

    % ln -s ../../data/relative_observables/lsu3shell_rme_6Li_Nmax02/ lsu3shell_rme
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/16/17 (mac): Created.  Based on compute_unit_tensor_rmes.cpp.
  2/18/17 (mac): Implement normalization test.
  2/20/17 (mac): Branch off spncci diagonalization code.
  6/17/17 (aem): Updated to use new symplectic operator i/o
  9/27/17 (aem): Updated to use intrinic definition of lgi
  10/4/17 (aem): Adapted for all A (including A<6)
  1/28/21 (aem): Removed unnecessary dependencies on branching_u3s, branching_u3lsj and upcoupling
****************************************************************/

#include <cstdio>
#include <ctime>
#include <fstream>
#include <sys/resource.h>

#include "fmt/format.h"
#include "lgi/lgi_solver.h"
#include "mcutils/profiling.h"
#include "mcutils/eigen.h"
#include "spncci/explicit_construction.h"

////////////////////////////////////////////////////////////////
// explicit construction checks
////////////////////////////////////////////////////////////////

void 
CheckOrthonormalityExplicit(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const basis::OperatorBlocks<double>& spncci_expansions,
    double tolerance
  )
// Check orthonormality of SpNCCI basis vectors from explicit
// expansion in lsu3shell basis.
//
// Takes pairs of BabySpNCCI subspaces sharing the same U3SPN, and
// thus the same underlying lsu3shell subspace.
//
// Arguments:
//   ...
{

  for (int bra_subspace_index=0; bra_subspace_index<baby_spncci_space.size(); ++bra_subspace_index)
    for (int ket_subspace_index=bra_subspace_index; ket_subspace_index<baby_spncci_space.size(); ++ket_subspace_index)
      {
        // extract subspace info
        const spncci::BabySpNCCISubspace& bra_subspace = baby_spncci_space.GetSubspace(bra_subspace_index);
        const spncci::BabySpNCCISubspace& ket_subspace = baby_spncci_space.GetSubspace(ket_subspace_index);

        // short circuit if subspaces have different underlying lsu3shell subspaces
        if (not (bra_subspace.omegaSPN()==ket_subspace.omegaSPN()))
          continue;

        // calculate overlaps
        Eigen::MatrixXd overlap_matrix = spncci_expansions[bra_subspace_index].transpose()*spncci_expansions[ket_subspace_index];
        Eigen::MatrixXd overlap_matrix_minus_identity = overlap_matrix - Eigen::MatrixXd::Identity(overlap_matrix.rows(),overlap_matrix.cols());
        mcutils::ChopMatrix(overlap_matrix,tolerance);
        mcutils::ChopMatrix(overlap_matrix_minus_identity,tolerance);

        // check overlaps
        bool on_diagonal = (bra_subspace_index==ket_subspace_index);
        bool success = on_diagonal ? mcutils::IsZero(overlap_matrix_minus_identity)
          : mcutils::IsZero(overlap_matrix);
        
        std::cout
          << fmt::format(
              "  bra index {} labels {} ket index {} labels {}",
              bra_subspace_index,bra_subspace.LabelStr(),
              ket_subspace_index,ket_subspace.LabelStr()
            )
          << std::endl;
        std::cout << mcutils::FormatMatrix(overlap_matrix,"14.7e","  ") << std::endl;
        std::cout << fmt::format("  on_diagonal {}",on_diagonal)
                  << std::endl;
        std::cout << fmt::format("  {}",success ? "PASS" : "FAIL")
                  << std::endl;
        std::cout << mcutils::FormatMatrix(overlap_matrix,"8.5f","  ") << std::endl;
        std::cout << std::endl;
      }
}

////////////////////////////////////////////////////////////////
// run parameters
////////////////////////////////////////////////////////////////

struct RunParameters
// Structure to store input parameters for run.
//
// Data members:
//   A (int): Atomic mass.
//   ...
{

  // constructor
  RunParameters(); 

  // basis parameters
  int A;
  HalfInt Nsigma0;
  // int Nsigma0_ex_max;
  int N1v;
  int Nmax;

  // filenames
  std::string lsu3shell_rme_directory;
  std::string lsu3shell_basis_filename;
  std::string Brel_filename;
  std::string Arel_filename;
  std::string Nrel_filename;
  std::string relative_unit_tensor_filename_template;
};

RunParameters::RunParameters()
{
  // basis parameters
  
  std::array<int,2> nuclide;
  nuclide[0]=2;
  nuclide[1]=1;
  A =nuclide[0]+nuclide[1];

  bool intrinsic=true;
  Nsigma0 = lgi::Nsigma0ForNuclide(nuclide,intrinsic);
  std::cout<<"nuclide ("<<nuclide[0]<<","<<nuclide[1]<<")  Nsigma0 "<<Nsigma0<<std::endl;

  // Nsigma0_ex_max = 2;
  N1v = 1;
  Nmax = 6;

  lsu3shell_rme_directory = "lsu3shell_rme";
  lsu3shell_basis_filename = lsu3shell_rme_directory + "/" + "lsu3shell_basis.dat";
  Brel_filename = lsu3shell_rme_directory + "/" + fmt::format("Brel.rme",Nmax);
  Arel_filename = lsu3shell_rme_directory + "/" + fmt::format("Arel.rme",Nmax);
  Nrel_filename = lsu3shell_rme_directory + "/" + fmt::format("Nrel.rme",Nmax);
  relative_unit_tensor_filename_template = lsu3shell_rme_directory + "/" + "relative_unit_{:06d}.rme";

}


////////////////////////////////////////////////////////////////
// main body
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  
  // SU(3) caching
  u3::U3CoefInit();
  u3::UCoefCache u_coef_cache;
  u3::PhiCoefCache phi_coef_cache;
  u3::g_u_cache_enabled = true;

  // numerical parameter for certain calculations
  double tolerance=1e-6;

  // run parameters
  RunParameters run_parameters;

  // Eigen OpenMP multithreading mode
  Eigen::initParallel();
  Eigen::setNbThreads(0);  // disable Eigen internal multithreading

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis
  ////////////////////////////////////////////////////////////////

  std::cout << "Read lsu3shell basis..." << std::endl;

  // read lsu3shell basis (regroup into U3SPN subspaces)
  lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3ShellBasis(
      run_parameters.Nsigma0,run_parameters.lsu3shell_basis_filename,
      lsu3shell_basis_table,lsu3shell_basis_provenance,lsu3shell_space
    );

  ////////////////////////////////////////////////////////////////
  // solve for LGIs
  ////////////////////////////////////////////////////////////////

  std::cout << "Solve for LGIs..." << std::endl;

  // timing start
  mcutils::SteadyTimer timer_lgi;
  timer_lgi.Start();


    // diagnostics
    // std::cout << "Arel operator..." << std::endl;
    // std::cout << "Arel sectors" << std::endl;
    // std::cout << Aintr_sectors.DebugStr();
    // std::cout << "Arel matrices" << std::endl;
    // for (int sector_index=0; sector_index<Aintr_sectors.size(); ++sector_index)
    //   {
    //     std::cout << fmt::format("  sector {}",sector_index) << std::endl;
    //     std::cout << mcutils::FormatMatrix(Aintr_matrices[sector_index],"8.5f","  ") << std::endl;
    //   }


  u3shell::SectorsU3SPN Bintr_sectors, Aintr_sectors, Nintr_sectors;
  basis::OperatorBlocks<double> Bintr_matrices, Aintr_matrices, Nintr_matrices;
  lsu3shell::ReadLSU3ShellSymplecticOperatorRMEs(
      lsu3shell_basis_table,lsu3shell_space, 
      run_parameters.Brel_filename,Bintr_sectors,Bintr_matrices,
      run_parameters.Nrel_filename,Nintr_sectors,Nintr_matrices,
      run_parameters.A
    );

  lsu3shell::ReadLSU3ShellSymplecticRaisingOperatorRMEs(
      lsu3shell_basis_table,lsu3shell_space, 
      run_parameters.Arel_filename,Aintr_sectors,Aintr_matrices,
      run_parameters.A
    );


  const u3shell::SectorsU3SPN& Ncm_sectors = Nintr_sectors;
  basis::OperatorBlocks<double> Ncm_matrices;
  lsu3shell::GenerateLSU3ShellNcmRMEs(
      lsu3shell_space,Nintr_sectors,Nintr_matrices,
      run_parameters.A-1,
      Ncm_matrices
    );

  
  lgi::MultiplicityTaggedLGIVector lgi_families;
  basis::OperatorBlocks<double> lgi_expansions;
  std::vector<int> lsu3shell_index_lookup_table;
  
  lgi::GenerateLGIExpansion(
      lsu3shell_space, 
      Bintr_sectors,Bintr_matrices,Ncm_sectors,Ncm_matrices,
      run_parameters.Nsigma0,
      lgi_families,lgi_expansions,
      lsu3shell_index_lookup_table
    );

  // diagnostics
  std::cout << fmt::format("  LGI families {}",lgi_families.size()) << std::endl;
  lgi::WriteLGILabels(lgi_families,std::cout);

  // timing stop
  timer_lgi.Stop();
  std::cout << fmt::format("(Task time: {})",timer_lgi.ElapsedTime()) << std::endl;

  ////////////////////////////////////////////////////////////////
  // set up SpNCCI space
  ////////////////////////////////////////////////////////////////

  std::cout << "Set up SpNCCI space..." << std::endl;

  // build SpNCCI irrep branchings
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // persistent container to store branchings
  spncci::NmaxTruncator truncator(run_parameters.Nsigma0,run_parameters.Nmax);

  // If A<6, construct restricted space 
  bool restrict_sp3r_to_u3_branching=false;
  if(run_parameters.A<6)
    restrict_sp3r_to_u3_branching=true;

  spncci::GenerateSpNCCISpace(lgi_families,truncator,spncci_space,sigma_irrep_map,restrict_sp3r_to_u3_branching);

  // put SpNCCI space into standard linearized container
  spncci::BabySpNCCISpace baby_spncci_space(spncci_space);

  // diagnostics
  std::cout << fmt::format("  Irrep families {}",spncci_space.size()) << std::endl;
  std::cout << fmt::format("  TotalU3Subspaces {}",spncci::TotalU3Subspaces(spncci_space)) << std::endl;
  std::cout << fmt::format("  TotalDimensionU3 {}",spncci::TotalDimensionU3S(spncci_space)) << std::endl;


  ////////////////////////////////////////////////////////////////
  // precompute K matrices
  ////////////////////////////////////////////////////////////////

  std::cout << "Precompute K matrices..." << std::endl;

  // timing start
  mcutils::SteadyTimer timer_k_matrices;
  timer_k_matrices.Start();

  // traverse distinct sigma values in SpNCCI space, generating K
  // matrices for each
  spncci::KMatrixCache k_matrix_cache;
  spncci::KMatrixCache kinv_matrix_cache;
  
  spncci::PrecomputeKMatrices(sigma_irrep_map,k_matrix_cache, kinv_matrix_cache, restrict_sp3r_to_u3_branching);

  // diagnostics
  for (const auto& sigma_irrep_pair : sigma_irrep_map)
    {
      // extract sigma and irrep contents
      const u3::U3& sigma = sigma_irrep_pair.first;
      const sp3r::Sp3RSpace& sp_irrep = sigma_irrep_pair.second;
      for (auto& omega_matrix_pair : k_matrix_cache[sigma])
        {
          const u3::U3& omega = omega_matrix_pair.first;
          const Eigen::MatrixXd& k_matrix = omega_matrix_pair.second;
          std::cout << fmt::format("  sigma {} omega {}",sigma.Str(),omega.Str()) << std::endl;
          std::cout << k_matrix << std::endl;
        }
    }

  // timing stop
  timer_k_matrices.Stop();
  std::cout << fmt::format("(Task time: {})",timer_k_matrices.ElapsedTime()) << std::endl;

  ////////////////////////////////////////////////////////////////
  // do explicit subspace constructions
  ////////////////////////////////////////////////////////////////

  std::cout << "Explicitly construct SpNCCI basis states using Arel..." << std::endl;
  basis::OperatorBlocks<double> spncci_expansions;
  
  spncci::ConstructSpNCCIBasisExplicit(
      lsu3shell_space,spncci_space,lgi_families,
      baby_spncci_space,k_matrix_cache,kinv_matrix_cache,
      Aintr_sectors,Aintr_matrices,spncci_expansions,
      restrict_sp3r_to_u3_branching
    );

  std::cout << "Check orthonormality for all SpNCCI subspaces sharing same underlying lsu3shell subspace..." << std::endl;
  CheckOrthonormalityExplicit(baby_spncci_space,spncci_expansions,tolerance);

}
