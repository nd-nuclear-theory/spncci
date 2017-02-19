/****************************************************************
  explicit.cpp

  Tests of explicit SpNCCI basis construction in LSU3Shell basis.

  Required data:

    Input files are generated by

      generate_lsu3shell_relative_operators.cpp

    which is invoked through scripting in

      compute_relative_tensors_lsu3shell_rmes.py

    Example: Z=3 N=3 twice_Nsigma0=22 Nmax=2 Nstep=2 Nv=N1b=1
   
    % python3 script/compute_relative_tensors_lsu3shell_rmes.py 3 3 22 2 2 1

    Only need .rme and .dat files.

    Not saved to repository since ~3.5 Mb...

       data/lsu3shell/lsu3shell_rme_6Li_Nmax02
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/16/17 (mac): Created.  Based on compute_unit_tensor_rmes.cpp.
****************************************************************/

#include <cstdio>
#include <ctime>
#include <fstream>
#include <sys/resource.h>

#include "am/am.h"
#include "cppformat/format.h"
#include "lgi/lgi.h"
#include "lgi/lgi_solver.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "mcutils/eigen.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h" 
#include "spncci/unit_tensor.h"
#include "spncci/spncci_branching_u3s.h"
#include "spncci/spncci_branching_u3lsj.h"
#include "u3shell/relative_operator.h"
#include "u3shell/upcoupling.h"

namespace spncci
{

  typedef std::unordered_map<u3::U3,vcs::MatrixCache,boost::hash<u3::U3>> KMatrixCache;
  // storage for K matrices
  //
  // maps sigma -> K matrix cache for that Sp irrep (vcs::MatrixCache),
  // where then vcs::MatrixCache maps omega to K matrix
  //
  // Usage: k_matrix_cache[sigma][omega]

  void
  PrecomputeKMatrices(
      const spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::KMatrixCache& k_matrix_cache
    )
  // Precompute and cache K matrices for all symplectic irreps
  // occurring in SpNCCI space.
  //
  // May be used either in SpNCCI RME recurrence or in explicit construction of states.
  //
  // Arguments:
  //   sigma_irrep_map (input): container for distinct symplectic irreps
  //   k_matrix_cache (output): container for corresponding K matrices
  {
    for (const auto& sigma_irrep_pair : sigma_irrep_map)
      {
        // extract sigma and irrep contents
        const u3::U3& sigma = sigma_irrep_pair.first;
        const sp3r::Sp3RSpace& sp_irrep = sigma_irrep_pair.second;

        // populate K matrix cache for this irrep
        vcs::GenerateKMatrices(sp_irrep,k_matrix_cache[sigma]);

        // diagnostics
        for (auto& omega_matrix_pair : k_matrix_cache[sigma])
          {
            const u3::U3& omega = omega_matrix_pair.first;
            const Eigen::MatrixXd& k_matrix = omega_matrix_pair.second;
            std::cout << fmt::format("  sigma {} omega {}",sigma.Str(),omega.Str()) << std::endl;
            std::cout << k_matrix << std::endl;
          }
      }

    }

  void 
  ConstructSpNCCIBasisExplicit(
      const u3shell::SpaceU3SPN& lsu3shell_space,
      const basis::MatrixVector& lgi_expansions,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::KMatrixCache& k_matrix_cache,
      const u3shell::SectorsU3SPN& Arel_sectors,
      const basis::MatrixVector& Arel_matrices,
      basis::MatrixVector& spncci_expansions
    )
  // Generate expansions of SpNCCI basis states in terms of lsu3shell
  // basis states, broken up by "baby SpNCCI subspaces", i.e., U3
  // subspaces segregated by irrep family.
  //
  // Limitation: Presently only supports Nex=0 or 2.
  //
  // Arguments:
  //   ...
  {

    spncci_expansions.resize(baby_spncci_space.size());
    for (int subspace_index=0; subspace_index<baby_spncci_space.size(); ++subspace_index)
      {
        // extract subspace properties
        const BabySpNCCISubspace& subspace = baby_spncci_space.GetSubspace(subspace_index);
        int irrep_family_index = subspace.irrep_family_index();
        int Nex = int(subspace.omega().N()-subspace.sigma().N());

        // define aliases to the relevant lsu3shell subspaces
        int lgi_lsu3shell_subspace_index = lsu3shell_space.LookUpSubspaceIndex(subspace.sigmaSPN());
        int spncci_lsu3shell_subspace_index = lsu3shell_space.LookUpSubspaceIndex(subspace.omegaSPN());

        // diagnostics
        const u3shell::SubspaceU3SPN& lgi_lsu3shell_subspace = lsu3shell_space.GetSubspace(lgi_lsu3shell_subspace_index);
        const u3shell::SubspaceU3SPN& spncci_lsu3shell_subspace = lsu3shell_space.GetSubspace(spncci_lsu3shell_subspace_index);
        std::cout
          << fmt::format(
              "Constructing subspace: LGI sigmaSPN {} family index {} => omegaSPN {} (Nex {})",
              subspace.sigmaSPN().Str(),irrep_family_index,
              subspace.omegaSPN().Str(),spncci_lsu3shell_subspace.size(),Nex
            )
          << std::endl
          << fmt::format(
              "  lsu3shell subspace sigmaSPN: U3SPN {} subspace index {} dim {}",
              lgi_lsu3shell_subspace.U3SPN().Str(),
              lgi_lsu3shell_subspace_index,
              lgi_lsu3shell_subspace.size()
            )
          << std::endl
          << fmt::format(
              "  lsu3shell subspace omegaSPN: U3SPN {} subspace index {} dim {}",
             spncci_lsu3shell_subspace.U3SPN().Str(),
             spncci_lsu3shell_subspace_index,
             spncci_lsu3shell_subspace.size()
            )
          << std::endl;

        // define aliases to expansion matrices
        const Eigen::MatrixXd& lgi_expansion = lgi_expansions[irrep_family_index];
        Eigen::MatrixXd& spncci_expansion = spncci_expansions[subspace_index];

        // diagnostics
        // std::cout << fmt::format("  lgi_expansion ({},{})",lgi_expansion.rows(),lgi_expansion.cols()) << std::endl;

        // calculate expansion of baby SpNCCI subspace
        assert((Nex==0)||(Nex==2));
        if (Nex==0)
          // Nex=0 -- trivial expansion
          {
            spncci_expansion = lgi_expansion;
          }
        else if (Nex==2)
          // Nex=1 -- single laddering, free of outer multiplicity
          {
            // retrieve matrix for applicable Arel sector
            int Arel_sector_index = Arel_sectors.LookUpSectorIndex(
                spncci_lsu3shell_subspace_index,
                lgi_lsu3shell_subspace_index
              );
            assert(Arel_sector_index!=basis::kNone);
            const Eigen::MatrixXd& Arel_matrix = Arel_matrices[Arel_sector_index];

            // diagnostics
            // std::cout << fmt::format("  Arel sector index {} matrix ({},{})",Arel_sector_index,Arel_matrix.rows(),Arel_matrix.cols()) << std::endl;

            // retrieve applicable inverse K matrix
            //
            // Language issue: Subscripted access only works for nonconst map.
            // const Eigen::MatrixXd& k_matrix = k_matrix_cache[subspace.sigma()][subspace.omega()];
            const Eigen::MatrixXd& k_matrix = k_matrix_cache.at(subspace.sigma()).at(subspace.omega());

            assert((k_matrix.rows()==1)&&(k_matrix.cols()==1));
            double k_inverse = 1 / k_matrix(0,0);  // specialize to multiplicity-free branching

            // act with raising operator on LGI expansions
            int phase = ParitySign(
                u3::ConjugationGrade(subspace.sigma().SU3())
                +u3::ConjugationGrade(subspace.omega().SU3())
              );
            
            spncci_expansion = phase*k_inverse*Arel_matrix*lgi_expansion;
          }

      }
  }

  void 
  CheckOrthonormalityExplicit(
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const basis::MatrixVector& spncci_expansions
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
        mcutils::ChopMatrix(overlap_matrix);
        mcutils::ChopMatrix(overlap_matrix_minus_identity);

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
        std::cout << fmt::format("  on_diagonal {}",on_diagonal)
                  << std::endl;
        std::cout << fmt::format("  {}",success ? "PASS" : "FAIL")
                  << std::endl;
        std::cout << mcutils::FormatMatrix(overlap_matrix,"8.5f","  ") << std::endl;
        std::cout << std::endl;
      }
  }


  void
  ReadLSU3ShellSeedUnitTensorRMEs(
      const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      std::vector<u3shell::SectorsU3SPN>& lgi_unit_tensor_sectors,
      std::vector<basis::MatrixVector>& lgi_unit_tensor_lsu3shell_matrices,
      const std::string& relative_unit_tensor_filename_template
    )
  // Read lsu3shell RMEs for seed unit tensors.
  //
  // Arguments:
  //   lsu3shell_basis_table (input): lsu3shell basis data
  //   lsu3shell_space (input): lsu3shell basis
  //   lgi_unit_tensor_lables (input): labels for the unit tensors to read
  //   lgi_unit_tensor_sectors (output): U3SPN sectors (for each unit tensor)
  //   lgi_unit_tensor_lsu3shell_matrices (output): matrices of RMEs (for each unit tensor)
  //   relative_unit_tensor_filename_template (input): filename template for use with fmt::format
  {
    lgi_unit_tensor_sectors.resize(lgi_unit_tensor_labels.size());
    lgi_unit_tensor_lsu3shell_matrices.resize(lgi_unit_tensor_labels.size());
    for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
      {
        // set up aliases for current unit tensor
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];
        u3shell::SectorsU3SPN& unit_tensor_sectors = lgi_unit_tensor_sectors[unit_tensor_index];
        basis::MatrixVector& unit_tensor_lsu3shell_matrices = lgi_unit_tensor_lsu3shell_matrices[unit_tensor_index];
      
        // read rmes
        const bool spin_scalar = false;
        std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);
        unit_tensor_sectors = u3shell::SectorsU3SPN(lsu3shell_space,unit_tensor_labels,spin_scalar);
        lsu3shell::ReadLSU3ShellRMEs(
            filename,
            lsu3shell_basis_table,lsu3shell_space,
            unit_tensor_labels,unit_tensor_sectors,unit_tensor_lsu3shell_matrices
          );
      }
  }

  void
  TransformSeedUnitTensorRMEs(
      const basis::MatrixVector& lgi_expansions,
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const std::vector<u3shell::SectorsU3SPN>& lgi_unit_tensor_sectors,
      const std::vector<basis::MatrixVector>& lgi_unit_tensor_lsu3shell_matrices,
      std::vector<basis::MatrixVector>& lgi_unit_tensor_spncci_matrices
    )
  // Transform seed lsu3shell sector RMEs to obtain RMEs between LGIs.
  //
  // Arguments:
  //   lgi_expansions (input): expansions of LGIs in lsu3shell basis
  //   lgi_unit_tensor_lables (input): labels for the unit tensors to read
  //   lgi_unit_tensor_sectors (input): U3SPN sectors (for each unit tensor)
  //   lgi_unit_tensor_lsu3shell_matrices (input): matrices of RMEs (for each unit tensor)
  //   lgi_unit_tensor_spncci_matrices (output): matrices of RMEs with respect to LGIs (for each unit tensor)
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
  // Store seed RMEs in souffle (by irrep family) for recurrence.
  //
  // Arguments:
  //   lgi_unit_tensor_lables (input): labels for the unit tensors to read
  //   lgi_unit_tensor_sectors (input): U3SPN sectors (for each unit tensor)
  //   lgi_unit_tensor_spncci_matrices (input): matrices of RMEs with respect to LGIs (for each unit tensor)
  //   unit_tensor_matrices (output): container for RMEs from recurrence
  //   zero_threshold (input): floating-point threshold for zero-suppression of sector
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


}// end namespace

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
  HalfInt Nsigma_0;
  int Nsigma0_ex_max;
  int N1b;
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
  // read from command line arguments
  //
  // TODO reorder filenames 
  // if (argc<8)
  //   {
  //     std::cout << "Syntax: A twice_Nsigma0 Nsigma0_ex_max N1B Nmax <basis filename> <Nrel filename> <Brel filename> <Arel filename>" 
  //               << std::endl;
  //     std::exit(1);
  //   }
  // int A = std::stoi(argv[1]); 
  // int twice_Nsigma0= std::stoi(argv[2]);
  // int Nsigma0_ex_max=std::stoi(argv[3]);
  // int N1b=std::stoi(argv[4]);
  // int Nmax = std::stoi(argv[5]);
  // std::string lsu3shell_basis_filename = argv[6];
  // std::string Nrel_filename = argv[7];
  // std::string Brel_filename = argv[8];
  // std::string Arel_filename = argv[9];
  // HalfInt Nsigma_0=HalfInt(twice_Nsigma0,2);

  // basis parameters
  A = 6;
  int twice_Nsigma0 = 22;
  Nsigma_0=HalfInt(twice_Nsigma0,2);

  Nsigma0_ex_max = 2;
  N1b = 1;
  Nmax = 2;
  lsu3shell_rme_directory = "lsu3shell_rme";
  lsu3shell_basis_filename = lsu3shell_rme_directory + "/" + "lsu3shell_basis.dat";
  Brel_filename = lsu3shell_rme_directory + "/" + "Brel_06_Nmax02.rme";
  Arel_filename = lsu3shell_rme_directory + "/" + "Arel_06_Nmax02.rme";
  Nrel_filename = lsu3shell_rme_directory + "/" + "Nrel_06_Nmax02.rme";
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
  double zero_threshold=1e-6;

  // run parameters
  RunParameters run_parameters;

  // Eigen OpenMP multithreading mode
  Eigen::initParallel();
  Eigen::setNbThreads(0);  // disable Eigen internal multithreading

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis and operators
  ////////////////////////////////////////////////////////////////

  std::cout << "Read lsu3shell inputs..." << std::endl;

  // read lsu3shell basis (regroup into U3SPN subspaces)
  lsu3shell::LSU3BasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3Basis(
      run_parameters.Nsigma_0,run_parameters.lsu3shell_basis_filename,
      lsu3shell_basis_table,lsu3shell_basis_provenance,lsu3shell_space
    );

  // read Brel
  u3shell::OperatorLabelsU3ST Brel_labels(-2,u3::SU3(0,2),0,0,0);
  u3shell::SectorsU3SPN Brel_sectors(lsu3shell_space,Brel_labels,true);
  basis::MatrixVector Brel_matrices;
  lsu3shell::ReadLSU3ShellRMEs(
      run_parameters.Brel_filename,
      lsu3shell_basis_table,lsu3shell_space,
      Brel_labels,Brel_sectors,Brel_matrices
    );

  // read Arel
  u3shell::OperatorLabelsU3ST Arel_labels(2,u3::SU3(2,0),0,0,0);
  u3shell::SectorsU3SPN Arel_sectors(lsu3shell_space,Arel_labels,true);
  basis::MatrixVector Arel_matrices;
  lsu3shell::ReadLSU3ShellRMEs(
      run_parameters.Arel_filename,
      lsu3shell_basis_table,lsu3shell_space,
      Arel_labels,Arel_sectors,Arel_matrices
    );

  // diagnostics
  std::cout << "Arel operator..." << std::endl;
  std::cout << "Arel sectors" << std::endl;
  std::cout << Arel_sectors.DebugStr();
  std::cout << "Arel matrices" << std::endl;
  for (int sector_index=0; sector_index<Arel_sectors.size(); ++sector_index)
    {
      std::cout << fmt::format("  sector {}",sector_index) << std::endl;
      std::cout << mcutils::FormatMatrix(Arel_matrices[sector_index],"8.5f","  ") << std::endl;
    }


  // read Nrel
  u3shell::OperatorLabelsU3ST Nrel_labels(0,u3::SU3(0,0),0,0,0);
  u3shell::SectorsU3SPN Nrel_sectors(lsu3shell_space,Nrel_labels,true);
  basis::MatrixVector Nrel_matrices;
  lsu3shell::ReadLSU3ShellRMEs(
      run_parameters.Nrel_filename,
      lsu3shell_basis_table,lsu3shell_space,
      Nrel_labels,Nrel_sectors,Nrel_matrices
    );

  // read Nrel as Ncm
  u3shell::OperatorLabelsU3ST Ncm_labels(0,u3::SU3(0,0),0,0,0);
  u3shell::SectorsU3SPN Ncm_sectors(lsu3shell_space,Ncm_labels,true);
  basis::MatrixVector Ncm_matrices;
  lsu3shell::GenerateLSU3ShellNcmRMEs(
      run_parameters.A,
      run_parameters.Nrel_filename,
      lsu3shell_basis_table,lsu3shell_space,
      Ncm_matrices
    );

  ////////////////////////////////////////////////////////////////
  // solve for LGIs
  ////////////////////////////////////////////////////////////////

  std::cout << "Solve for LGIs..." << std::endl;

  lgi::MultiplicityTaggedLGIVector lgi_families;
  basis::MatrixVector lgi_expansions;
  lgi::GenerateLGIExpansion(
      lsu3shell_space, 
      Brel_sectors,Brel_matrices,Ncm_sectors,Ncm_matrices,
      run_parameters.Nsigma_0,
      lgi_families,lgi_expansions
    );

  // diagnostics
  std::cout << fmt::format("  LGI families {}",lgi_families.size()) << std::endl;
  lgi::WriteLGILabels(lgi_families,std::cout);

  ////////////////////////////////////////////////////////////////
  // set up SpNCCI space
  ////////////////////////////////////////////////////////////////

  std::cout << "Set up SpNCCI space..." << std::endl;

  // build SpNCCI irrep branchings
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // persistent container to store branchings
  spncci::NmaxTruncator truncator(run_parameters.Nsigma_0,run_parameters.Nmax);
  spncci::GenerateSpNCCISpace(lgi_families,truncator,spncci_space,sigma_irrep_map);

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

  // traverse distinct sigma values in SpNCCI space, generating K
  // matrices for each
  spncci::KMatrixCache k_matrix_cache;
  spncci::PrecomputeKMatrices(sigma_irrep_map,k_matrix_cache);

  ////////////////////////////////////////////////////////////////
  // do explicit subspace constructions
  ////////////////////////////////////////////////////////////////

  std::cout << "Explicitly construct SpNCCI basis states using Arel..." << std::endl;
  basis::MatrixVector spncci_expansions;
  spncci::ConstructSpNCCIBasisExplicit(
      lsu3shell_space,lgi_expansions,baby_spncci_space,k_matrix_cache,
      Arel_sectors,Arel_matrices,spncci_expansions
    );

  std::cout << "Check orthonormality for all SpNCCI subspaces sharing same underlying lsu3shell subspace..." << std::endl;
  spncci::CheckOrthonormalityExplicit(baby_spncci_space,spncci_expansions);

  ////////////////////////////////////////////////////////////////
  // read lsu3shell seed unit tensor rmes
  ////////////////////////////////////////////////////////////////

  std::cout << "Read seed unit tensor rmes..." << std::endl;

  // storage for seed unit tensor rmes
  //
  //   lgi_unit_tensor_labels: vector of labels for seed unit tensors
  //   lgi_unit_tensor_lsu3shell_sectors: vector of lsu3shell sectors for seed unit tensors
  //   lgi_unit_tensor_matrices: vector of matrices for these sectors
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensor_labels;
  std::vector<u3shell::SectorsU3SPN> lgi_unit_tensor_sectors;
  std::vector<basis::MatrixVector> lgi_unit_tensor_lsu3shell_matrices;

  // determine set of seed unit tensors
  //
  // i.e., those for which we calculate seed rmes among the LGIs
  //
  // Note: Should be consistant with set of tensors generated by
  // generate_lsu3shell_relative_operators.
  int Nmax_for_unit_tensors = run_parameters.Nsigma0_ex_max+2*run_parameters.N1b;  // max quanta for pair in LGI (?)
  int J0 = -1;  // all J0
  int T0 = 0;
  const bool restrict_positive_N0 = false;  // don't restrict to N0 positive
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
      Nmax_for_unit_tensors, lgi_unit_tensor_labels,
      J0,T0,restrict_positive_N0
    );

  // diagnostic
  std::cout << fmt::format("  seed unit tensors {}",lgi_unit_tensor_labels.size()) << std::endl;

  spncci::ReadLSU3ShellSeedUnitTensorRMEs(
      lsu3shell_basis_table,lsu3shell_space,
      lgi_unit_tensor_labels,
      lgi_unit_tensor_sectors,
      lgi_unit_tensor_lsu3shell_matrices,
      run_parameters.relative_unit_tensor_filename_template
    );

  ////////////////////////////////////////////////////////////////
  // transform and store seed rmes for use in SpNCCI recurrence
  ////////////////////////////////////////////////////////////////

  std::cout << "Transform and store seed unit tensor rmes..." << std::endl;

  // transform to SpNCCI LGI RMEs
  std::vector<basis::MatrixVector> lgi_unit_tensor_spncci_matrices;
  spncci::TransformSeedUnitTensorRMEs(
      lgi_expansions,
      lgi_unit_tensor_labels,
      lgi_unit_tensor_sectors,
      lgi_unit_tensor_lsu3shell_matrices,
      lgi_unit_tensor_spncci_matrices
    );

  // store unit tensor matrix elements for recurrence
  spncci::UnitTensorMatricesByIrrepFamily unit_tensor_matrices;
  spncci::StoreSeedUnitTensorRMEs(
      lgi_unit_tensor_labels,
      lgi_unit_tensor_sectors,
      lgi_unit_tensor_spncci_matrices,
      unit_tensor_matrices,
      zero_threshold
    );

  ////////////////////////////////////////////////////////////////
  // ladder unit tensor rmes to full SpNCCI basis
  ////////////////////////////////////////////////////////////////

  // TODO
}
