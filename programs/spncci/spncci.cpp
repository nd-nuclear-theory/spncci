/****************************************************************
  spncci.cpp

  Tests of explicit SpNCCI basis construction in LSU3Shell basis.

  This code just tests normalization, but using clean refactored
  infrastructure.  Other deeper tests (of unit tensor matrix elements)
  were carried out in compute_unit_tensor_rmes.cpp.

  Required data:

  * Relative operator lsu3shell rme input files are generated by

      generate_lsu3shell_relative_operators.cpp

    which is invoked through scripting in

      compute_relative_tensors_lsu3shell_rmes.py

    Example: Z=3 N=3 twice_Nsigma0=22 Nmax=2 Nstep=2 N1v[=N1b]=1

    % python3 script/compute_relative_tensors_lsu3shell_rmes.py 3 3 22 2 2 1

    Only need .rme and .dat files.

    Not saved to repository since ~3.5 Mb...

       data/lsu3shell/lsu3shell_rme_6Li_Nmax02

       * Relative Hamiltonian (and observable) upcoupled rme files are
    generated by

      generate_relative_u3st_operators

    which is invoked manually for now as

      generate_relative_u3st_operators A Nmax N1v basename

   Example:

       ../operators/generate_relative_u3st_operators 6 2 1 hamiltonian

       with hamiltonian.load containing

       20    // hw
       Tintr 1.0    // coef
       INT 1.0 4 0 0 0 relative_observables/JISP16_Nmax20_hw20.0_rel.dat // coef Jmax J0 T0 g0 interaction_filename

       ../operators/generate_relative_u3st_operators 6 2 1 Nintr

       with Nintr.load containing

       20    // hw
       Nintr 1.0    // coef

       ../operators/generate_relative_u3st_operators 6 2 1 r2intr

       with r2intr.load containing

       20    // hw
       r2intr 1.0    // coef

   % ln -s ../../data/relative_observables/




  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  2/20/17 (mac): Created (starting from explicit.cpp).
  4/9/17 (aem): Incorporated baby spncci hypersectors
  6/5/17 (mac): Read relative rather than intrinsic symplectic operators.
  6/16/17 (aem) : offload to computation and io control
  10/4/17 (aem) : Fixed basis construction and recurrence for A<6
  1/16/18 (aem) : Offloaded explicit construction and recurrence
    checks to explicit_construction.h
  1/30/18 (aem): Overhalled seed generation and recurrence
  2/5/18 (aem): Switched from using u3s sectors to u3s hypersectors
    combined with observable spaces
  2/15/18 (aem) : Removed gamma_max=0 lgi
    + Cleaned up codes and factored spncci.cpp into simpler functions
  5/2/19 (aem) : Moved ComputeManyBodyRMEs into hyperblocks_u3s
  6/21/19 (aem) : Extracted variance calculation functions into variance.{h,cpp}
  11/5/19 (aem) : Stripped out truncation tests and bundled initialization
  8/26/20 (aem) : Removed dependencies on obselete basis classes SpaceSpU3S, SpaceSpLS,SpaceSpJ
                  and update basis statistics written to results file

Notes:
branching2 currently used for branching. 
****************************************************************/

#include <cstdio>
#include <ctime>
#include <fstream>
#include <istream>
#include <iostream>
#include <sys/resource.h>
#include <omp.h>
#include <unordered_set>

#include "Spectra/SymEigsSolver.h"  // from spectra
#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "mcutils/io.h"
#include "spncci/recurrence.h"
#include "lgi/lgi_unit_tensors.h"
#include "mcutils/profiling.h"
#include "mcutils/eigen.h"
#include "spncci/computation_control.h"
#include "spncci/decomposition.h"
#include "spncci/eigenproblem.h"
#include "spncci/explicit_construction.h"
#include "spncci/io_control.h"
#include "spncci/computation_control.h"
#include "spncci/results_output.h"
#include "spncci/transform_basis.h"
#include "spncci/vcs_cache.h"
#include "spncci/hyperblocks_u3s.h"
#include "spncci/variance.h"
#include "spncci/parameters.h"
////////////////////////////////////////////////////////////////
// WIP code
//
// to extract to spncci library when ready
////////////////////////////////////////////////////////////////
namespace spncci
{
  void WriteSymmetricOperatorMatrix(
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const HalfInt& J0,
      const spncci::SpaceSpBasis& spbasis_bra, //For a given J
      const spncci::SpaceSpBasis& spbasis_ket, //For a given J
      const std::vector<spncci::LGIPair>& lgi_pairs_recurrence,
      int observable_index, int hw_index,
      const std::string& filename
    )
  {
    // open output file
    std::ofstream out_stream(filename);

    // Get J of bases
    HalfInt Jp=spbasis_bra.J();
    HalfInt J=spbasis_ket.J();

    //From list of lgi pairs from the recurrence (for which there are non-zero RMEs)
    //create a lookup table for checking if a given lgi pair in basis corresponds to a non-zero tile
    std::unordered_set<spncci::LGIPair,boost::hash<spncci::LGIPair>> lgi_pair_lookup_table;
    for(spncci::LGIPair pair : lgi_pairs_recurrence)
      lgi_pair_lookup_table.insert(pair);

    std::vector<std::vector<int>> offsets_bra;
    std::vector<std::vector<int>> offsets_ket;
    spncci::GetSpBasisOffsets(spbasis_bra,offsets_bra);
    spncci::GetSpBasisOffsets(spbasis_ket,offsets_ket);

    u3::WCoefCache w_cache;

    // Iterate through bases for bra and ket.  Each subspaces corresponds to a single irrep family (irrep subspace)
    for(int bra_subspace_index=0; bra_subspace_index<spbasis_bra.size(); ++bra_subspace_index)
      for(int ket_subspace_index=0; ket_subspace_index<spbasis_ket.size(); ++ket_subspace_index)
        {
          const spncci::SubspaceSpBasis& spbasis_subspace_bra=spbasis_bra.GetSubspace(bra_subspace_index);
          const spncci::SubspaceSpBasis& spbasis_subspace_ket=spbasis_ket.GetSubspace(ket_subspace_index);

          const int irrep_family_index_bra=spbasis_subspace_bra.irrep_family_index();
          const int irrep_family_index_ket=spbasis_subspace_ket.irrep_family_index();

          const int tile_dimension_bra=spbasis_subspace_bra.dimension();
          const int tile_dimension_ket=spbasis_subspace_ket.dimension();

          // If ket>bra, then need to construct adjoint tile and take transpose
          bool adjoint=irrep_family_index_ket>irrep_family_index_bra;

          // In some cases, no states in the irrep will contribute to a given J subspace
          // Usually only happens for low Nmax calculations of if irrep is a high Nsex irrep.
          if(tile_dimension_ket==0 || tile_dimension_bra==0)
            continue;

          spncci::LGIPair lgi_pair;
          lgi_pair=adjoint?spncci::LGIPair(irrep_family_index_ket,irrep_family_index_bra):
                            spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket);

          spncci::OperatorBlock tile;

          // If lgi pair is in lookup table,
          if(lgi_pair_lookup_table.count(lgi_pair)==0)
            {
              //Construct tile of zeros
              tile=Eigen::MatrixXd::Zero(tile_dimension_bra,tile_dimension_ket);
            }
          else
            {
              /////////////////////////////////////////////////////////////////////////////////////////////////
              //// Open files containing hyperblocks and hypersectors interating over thread number
              /////////////////////////////////////////////////////////////////////////////////////////////////
              std::vector<spncci::ObservableHypersectorLabels> list_baby_spncci_hypersectors;
              basis::OperatorHyperblocks<double> baby_spncci_observable_hyperblocks;

              spncci::GetBabySpNCCIHyperBlocks(
                observable_index,hw_index,lgi_pair,
                list_baby_spncci_hypersectors,
                baby_spncci_observable_hyperblocks
                );

              //Look up offset vectors
              const std::vector<int>& offsets_bra_subspace=offsets_bra[bra_subspace_index];
              const std::vector<int>& offsets_ket_subspace=offsets_ket[ket_subspace_index];

              // If adjoint=true, then get tile for LGI pair (ket_lgi,bra_lgi) and transpose,
              // otherwise, just get tile corresponding to LGI pair
              if(adjoint)
                {
                  spncci::OperatorBlock temp_tile;
                  spncci::GetOperatorTile(
                    baby_spncci_space,observable_space,spbasis_subspace_ket,spbasis_subspace_bra,
                    offsets_ket_subspace,offsets_bra_subspace,J0,J,Jp,hw_index,observable_index,
                    lgi_pair,w_cache,list_baby_spncci_hypersectors,baby_spncci_observable_hyperblocks,
                    temp_tile
                  );

                  tile=temp_tile.transpose();
                }
              else
                {
                  spncci::GetOperatorTile(
                    baby_spncci_space,observable_space,spbasis_subspace_bra,spbasis_subspace_ket,
                    offsets_bra_subspace,offsets_ket_subspace,J0,Jp,J,hw_index,observable_index,
                    lgi_pair,w_cache,list_baby_spncci_hypersectors,baby_spncci_observable_hyperblocks,
                    tile
                  );
                }
            }

          // write tile to file
          mcutils::WriteBinary<int32_t>(out_stream, 2*sizeof(int32_t));
          mcutils::WriteBinary<int32_t>(out_stream, tile_dimension_bra);
          mcutils::WriteBinary<int32_t>(out_stream, tile_dimension_ket);
          mcutils::WriteBinary<int32_t>(out_stream, 2*sizeof(int32_t));
          int32_t num_matrix_elements = tile_dimension_bra*tile_dimension_ket;
          mcutils::WriteBinary<int32_t>(out_stream, num_matrix_elements*sizeof(double));
          mcutils::WriteBinary<double>(out_stream, tile.data(), num_matrix_elements);
          mcutils::WriteBinary<int32_t>(out_stream, num_matrix_elements*sizeof(double));
        }
    out_stream.close();
  }

}// end namespace

////////////////////////////////////////////////////////////////
// main body
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  std::cout<<"entering spncci"<<std::endl;
  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  //Initializes extern variables and calls Eigen::initParallel() and u3::U3CoefInit()
  spncci::InitializeSpNCCI();

  std::cout << "Reading control file..." << std::endl;
  spncci::RunParameters run_parameters;
  std::cout<<"Nmax="<<run_parameters.Nmax<<std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////
  //// set up SpNCCI spaces
  ///////////////////////////////////////////////////////////////////////////////////////
  lgi::MultiplicityTaggedLGIVector lgi_families;
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map; //Container for spncci space irreps.  Must stay in scope for spncci_space 
  spncci::BabySpNCCISpace baby_spncci_space;

  std::vector<spncci::SpaceSpBasis> spaces_spbasis;
  spncci::KMatrixCache k_matrix_cache, kinv_matrix_cache;

  // 
  bool verbose_basis=true;
  spncci::SetUpSpNCCISpaces(
      run_parameters,lgi_families,spncci_space,sigma_irrep_map,baby_spncci_space,
      spaces_spbasis,k_matrix_cache, kinv_matrix_cache,verbose_basis
    );

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Results file: Writing parameters and basis statistics
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  // open output files
  std::ofstream results_stream("spncci.res");

  // results output: code information
  spncci::StartNewSection(results_stream,"CODE");
  spncci::WriteCodeInformation(results_stream,run_parameters);

  // results output: run parameters
  spncci::StartNewSection(results_stream,"PARAMETERS");
  spncci::WriteRunParameters(results_stream,run_parameters);

  // results output: basis information
  spncci::StartNewSection(results_stream,"BASIS");
  // spncci::WriteBasisStatistics(results_stream,spncci_space,baby_spncci_space,spu3s_space,spls_space,spj_space);
  spncci::WriteBasisStatistics(results_stream,spncci_space,baby_spncci_space,spaces_spbasis);
  // spncci::WriteSpU3SSubspaceListing(results_stream,baby_spncci_space,run_parameters.Nsigma0);
  spncci::WriteBabySpNCCISubspaceListing(results_stream,baby_spncci_space,run_parameters.Nsigma0);


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Enumerate unit tensor space
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  int J0_for_unit_tensors = -1;  // all J0
  int T0_for_unit_tensors = -1;  // all T0
  const bool restrict_positive_N0 = false;  // don't restrict to N0 positive

  // get full set of possible unit tensor labels up to Nmax, N1v truncation
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
      run_parameters.Nmax, run_parameters.N1v,
      unit_tensor_labels,J0_for_unit_tensors,T0_for_unit_tensors,
      restrict_positive_N0
    );

  // generate unit tensor subspaces
  u3shell::RelativeUnitTensorSpaceU3S
    unit_tensor_space(run_parameters.Nmax,run_parameters.N1v,unit_tensor_labels);

  // std::cout<<"unit tensor space "<<unit_tensor_space.size()<<std::endl;
  // std::cout<<unit_tensor_space.Str()<<std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //  Read in observables
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout << "Reading observables..." << std::endl;

  // Initialize containers for rmes and their symmetries
  // Stored by hw, then by observable
  std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>> observables_relative_rmes(run_parameters.hw_values.size());
  std::vector<std::vector<u3shell::IndexedOperatorLabelsU3S>> observable_symmetries_u3s(run_parameters.num_observables);

  spncci::ReadRelativeObservables(
      run_parameters.Nmax, run_parameters.N1v, run_parameters.hw_values,
      run_parameters.observable_directory,run_parameters.observable_filenames,
      unit_tensor_space, observables_relative_rmes, observable_symmetries_u3s
    );

  // Create observable spaces for each observable including Hamiltonian
  std::cout<<"create observable space"<<std::endl;
  std::vector<u3shell::ObservableSpaceU3S> observable_spaces(run_parameters.num_observables);
  for(int ob_num=0; ob_num<run_parameters.num_observables; ++ob_num)
      observable_spaces[ob_num]=u3shell::ObservableSpaceU3S(observable_symmetries_u3s[ob_num]);

  ////////////////////////////////////////////////////////////////
  // terminate counting only run
  ////////////////////////////////////////////////////////////////
  if (run_parameters.count_only)
    {
      // termination
      results_stream.close();

      std::cout << "End of counting-only run" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"setting up lgi unit tensor blocks"<<std::endl;
  // Get list of unit tensor labels between lgi's
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
      run_parameters.Nsigmamax, run_parameters.N1v,
      lgi_unit_tensor_labels,J0_for_unit_tensors,T0_for_unit_tensors,
      restrict_positive_N0
    );

  //Get look-up table for lgi index in full space.  Used for looking up seed filenames
  // which are index by full space index
  std::cout<<"reading lgi table "<<std::endl;
  std::vector<int> lgi_full_space_index_lookup;
  lgi::ReadLGILookUpTable(lgi_full_space_index_lookup,lgi_families.size());
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"Starting recurrence and contraction"<<std::endl;
  // Get list of lgi pairs with non-zero matrix elements between them.
  // Restricted to ket<=bra.
  //
  std::vector<spncci::LGIPair> lgi_pairs;
  spncci::GetLGIPairsForRecurrence(lgi_full_space_index_lookup,spncci_space,run_parameters.Nmax,lgi_pairs);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // If transforming LGI basis
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  spncci::OperatorBlocks lgi_transformations;
  if(run_parameters.transform_lgi)
    {
      std::cout<<"reading in lgi transformations"<<std::endl;
      std::string lgi_transformations_filename="lgi_transformations.dat";
      spncci::ReadTransformationMatrices(lgi_transformations_filename,lgi_transformations);
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"begin parallel region"<<std::endl;

  //Declaring shared variables
  int num_files;

  mcutils::SteadyTimer timer_recurrence;
  timer_recurrence.Start();

  std::ofstream timing_output;
  timing_output.open("recurrence_timing.dat");

  #pragma omp parallel shared(num_files)
    {
      // Parallelization is currently set up so that each thread needs at least on lgi pair
      #pragma omp single
        {
          int num_threads=omp_get_num_threads();
          if(num_threads>lgi_pairs.size())
            {
              std::cout<<"Too many threads.  Only "<<lgi_pairs.size()<<" needed."<<std::endl;
              assert(num_threads<=lgi_pairs.size());
            }

          num_files=num_threads;
        }

      //private coefficient caches--avoids locks and barriers
      u3::UCoefCache u_coef_cache;
      u3::PhiCoefCache phi_coef_cache;

      #pragma omp for schedule(dynamic) nowait
      // For each LGI pair, compute SU(3)xSU(2) reduced many-body matrix elements of unit tensors.
      // Then contract unit tensors with relative matrix elements of observables
      // Write observable hypersectors and hyperblocks to separate file for each LGI pair, each
      // observable and each hw.
      //
      // Note: Only observable hypersectors with irrep_family_bra>=irrep_family_ket written to files
      // If diagonal sector, only upper triangle stored.  --Is this still correct?
      for(int i=0; i<lgi_pairs.size(); ++i)
        {
          const spncci::LGIPair& lgi_pair=lgi_pairs[i];
          mcutils::SteadyTimer timer_pair;
          int bra,ket;
          std::tie(bra,ket)=lgi_pair;

          timer_pair.Start();
          #pragma omp critical
          	timing_output<<fmt::format("Computing for lgi pair {:3d} : ({:3d},{:3d})",i,bra,ket)<<std::endl;

          spncci::ComputeManyBodyRMEs(
              run_parameters,lgi_families,lgi_full_space_index_lookup,
              spncci_space,baby_spncci_space,unit_tensor_space, observable_spaces,
              observables_relative_rmes,k_matrix_cache,kinv_matrix_cache,
              lgi_transformations,u_coef_cache,phi_coef_cache,lgi_pair
            );
          timer_pair.Stop();

          std::string out_str=fmt::format("Recurrence for pair {:3d} : ({:3d},{:3d}) took {}",i, bra,ket,timer_pair.ElapsedTime());
          #pragma omp critical
            timing_output<<out_str<<std::endl;

        }// end lgi_pair

      //After recurrence completed, dealocate coefficient caches
      // std::cout<<"ucoef cache size: "<<u_coef_cache.size()<<std::endl;
      u_coef_cache.clear();
      phi_coef_cache.clear();

    } //end parallel region
     
  timing_output.close();
  timer_recurrence.Stop();
  std::cout<<"Recurrence: "<<timer_recurrence.ElapsedTime()<<std::endl;
  ////////////////////////////////////////////////////////////////
  // calculation mesh master loop
  ////////////////////////////////////////////////////////////////

  // timing start
  mcutils::SteadyTimer timer_mesh;
  timer_mesh.Start();

  // W coefficient cache -- needed for observable branching
  u3::WCoefCache w_cache;

  // for each hw value, solve eigenproblem and get expectation values
  std::cout << "Calculation mesh master loop..." << std::endl;
  for(int hw_index=0; hw_index<run_parameters.hw_values.size(); ++hw_index)
    {
      // retrieve mesh parameters
      double hw = run_parameters.hw_values[hw_index];

      // results output: log start of individual mesh calculation
      spncci::StartNewSection(results_stream,"RESULTS");
      spncci::WriteCalculationParameters(results_stream,hw);

      ////////////////////////////////////////////////////////////////
      // eigenproblem
      ////////////////////////////////////////////////////////////////
      std::cout<<"Solve eigenproblem..."<<std::endl;

      std::vector<spncci::Vector> eigenvalues(run_parameters.J_values.size());  // eigenvalues by J subspace
      std::vector<spncci::Matrix> eigenvectors(run_parameters.J_values.size());  // eigenvectors by J subspace

      // Construct and diagonalize Hamiltonian, do decompositions
      {
        const int observable_index = 0;  // for Hamiltonian
        for(int subspace_index=0; subspace_index<run_parameters.J_values.size(); ++subspace_index)
          {
            // for eigenproblem
            const HalfInt& J=run_parameters.J_values[subspace_index];
            spncci::Vector& eigenvalues_J = eigenvalues[subspace_index];
            spncci::Matrix& eigenvectors_J = eigenvectors[subspace_index];
            HalfInt J00 = run_parameters.observable_J0_values[observable_index];

            const spncci::SpaceSpBasis& spbasis_bra=spaces_spbasis[subspace_index];
            const spncci::SpaceSpBasis& spbasis_ket=spaces_spbasis[subspace_index];

            const u3shell::ObservableSpaceU3S& observable_space=observable_spaces[observable_index];

            
            if(false)
              {
                // Write basis information to file
                //
                // Number of states in full basis
                // For each subspace (irrep family): Number of states in given subspace
                std::string basis_filename = fmt::format("basis-J{:03.1f}.dat", float(J));
                std::cout << "  Writing basis information to file " << basis_filename << std::endl;
                std::ostringstream basis_str_stream;
                std::size_t true_dimension = 0;
                for (std::size_t subspace_index=0; subspace_index < spbasis_bra.size(); ++subspace_index)
                {
                  if (spbasis_bra.GetSubspace(subspace_index).full_dimension() <= 0)
                    continue;
                  ++true_dimension;
                  basis_str_stream << fmt::format(
                      "{:d}", spbasis_bra.GetSubspace(subspace_index).full_dimension()
                    ) << std::endl;
                }
                std::ofstream basis_stream(basis_filename);
                basis_stream << fmt::format("{:d}", true_dimension) << std::endl;
                basis_stream << basis_str_stream.str();
                basis_stream.close();
              }


            //Start timer for Hamiltonian construction
            mcutils::SteadyTimer timer_hamiltonian;
            timer_hamiltonian.Start();

            std::cout<<"  Constructing Hamiltonian matrix"<<std::endl;

            spncci::OperatorBlock hamiltonian_matrix;
            spncci::ConstructSymmetricOperatorMatrix(
              baby_spncci_space,observable_space,
              J00,spbasis_bra, spbasis_ket, lgi_pairs,
              observable_index,hw_index,hamiltonian_matrix
            );


            if(false)
              { 
                std::string filename = fmt::format(
                "hamiltonian-hw{:04.1f}-J{:03.1f}",
                run_parameters.hw_values.at(hw_index), float(J)
                );
                std::cout<<"filename is "<<filename<<std::endl;
                
                spncci::WriteSymmetricOperatorMatrix(
                  baby_spncci_space,observable_space,
                  J00,spbasis_bra, spbasis_ket, lgi_pairs,
                  observable_index,hw_index,filename
                );
              }
            
            // assert(0);
            timer_hamiltonian.Stop();
            std::cout<<fmt::format("    time: {}",timer_hamiltonian.ElapsedTime())<<std::endl;

            std::cout << fmt::format("  Diagonalizing: J={}",J) << std::endl;

            mcutils::SteadyTimer timer_eigensolver;
            timer_eigensolver.Start();
            std::cout<<"num eigenvalues :"<<run_parameters.num_eigenvalues<<std::endl;
            std::cout<<"num convergence :"<<run_parameters.eigensolver_num_convergence<<std::endl;
            std::cout<<"basis dimension :"<<hamiltonian_matrix.rows()<<std::endl;
            spncci::SolveHamiltonian(
                hamiltonian_matrix,
                run_parameters.num_eigenvalues,
                run_parameters.eigensolver_num_convergence,  // whatever exactly this is...
                run_parameters.eigensolver_max_iterations,
                run_parameters.eigensolver_tolerance,
                eigenvalues_J,eigenvectors_J
              );
            timer_eigensolver.Stop();
            std::cout<<fmt::format("   time: {}",timer_eigensolver.ElapsedTime())<<std::endl;

            //Testing
            //Eigen matrix cast to float is not a good idea
            int binary_float_precision=8;
            std::string eigv_filename=fmt::format("eigenvector_{:02d}_{:02.1f}.dat",TwiceValue(J),hw);
            std::cout<<"write eigenvector"<<std::endl;
            spncci::WriteEigenvectors(eigenvectors_J,J,eigv_filename,binary_float_precision);
            //////////////////////////////////////////////////////////////////
          }

        // Write eigenvalues to results_stream
        std::cout<<"Writing eigenvalues"<<std::endl;
        spncci::WriteEigenvalues(results_stream,run_parameters.J_values,eigenvalues,run_parameters.gex);

        // Generate decompositions and write to results_stream
        std::cout<<"Generating decompositions"<<std::endl;
        // NOTE: Function now reads eigenvectors from file before doing decomposition.  
        // TODO: Extract to postprocessing code 
        spncci::GenerateDecompositions(baby_spncci_space,spaces_spbasis,run_parameters, hw,results_stream);

      }// End Hamiltonian section

      {
      //////////////////////////////////////////////////////////////
      // calculate observable RMEs
      ////////////////////////////////////////////////////////////////
      // TODO: Set flag to only do observable 0 (Hamiltonian) if debugging 

      std::cout << "Calculate observable results..." << std::endl;
      mcutils::SteadyTimer timer_observables;
      timer_observables.Start();

      // observable_results_matrices:
      //   - vector over observable_index
      //   - vector over sector_index
      //   - matrix over (bra_eigenstate_index,ket_eigenstate_index)
      std::vector<spncci::OperatorBlocks> observable_results_matrices;
      observable_results_matrices.resize(run_parameters.num_observables);

      //std::cout<<"Get allowed (Jp,J) for each observable"<<std::endl;
      std::vector<std::vector<std::pair<int,int>>> observable_sectors(run_parameters.num_observables);
      for (int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
        {
          auto& sectors=observable_sectors[observable_index];
          HalfInt J0 = run_parameters.observable_J0_values[observable_index];
          for(int i=0; i<run_parameters.J_values.size(); ++i)
            for(int j=0; j<run_parameters.J_values.size(); ++j)
              {
                HalfInt Jp=run_parameters.J_values[i];
                HalfInt J=run_parameters.J_values[j];
                if(am::AllowedTriangle(J,J0,Jp))
                  sectors.emplace_back(i,j);
              }
        }

      // std::cout<<"calculate observable expectation value for each allowed (Jp,J)"<<std::endl;
      for (int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
        {
          std::cout<<"Starting loop for observable "<<observable_index<<std::endl;
          HalfInt J0 = run_parameters.observable_J0_values[observable_index];
          auto&sectors=observable_sectors[observable_index];

          observable_results_matrices[observable_index].resize(sectors.size());

          const u3shell::ObservableSpaceU3S& observable_space=observable_spaces[observable_index];

          for(int sector_index=0; sector_index<sectors.size(); ++sector_index)
            {
              int bra_index,ket_index;
              std::tie(bra_index,ket_index)=sectors[sector_index];

              const spncci::SpaceSpBasis& spbasis_bra=spaces_spbasis[bra_index];
              const spncci::SpaceSpBasis& spbasis_ket=spaces_spbasis[ket_index];

              const HalfInt bra_J=run_parameters.J_values[bra_index];
              const HalfInt ket_J=run_parameters.J_values[ket_index];

              spncci::OperatorBlock observable_block;
              spncci::ConstructSymmetricOperatorMatrix(
                  baby_spncci_space,observable_space,
                  J0,spbasis_bra,spbasis_ket,lgi_pairs,
                  observable_index, hw_index,
                  observable_block
                );

              
              Eigen::MatrixXd& observable_results_matrix = observable_results_matrices[observable_index][sector_index];

              observable_results_matrix = eigenvectors[bra_index].transpose()
                * observable_block
                * eigenvectors[ket_index];


              std::cout
                << fmt::format("Observable {} bra_J {} ket_J {}",observable_index,bra_J,ket_J)
                << std::endl;
            }
        }

      // end timing
      timer_observables.Stop();
      std::cout << fmt::format("  (Observables: {})",timer_observables.ElapsedTime()) << std::endl;

      // results output: observables

      spncci::WriteObservables(
        results_stream,run_parameters.J_values,observable_sectors,
        observable_results_matrices,run_parameters.gex
      );

      }//End observable section 

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }

// timing stop
timer_mesh.Stop();
std::cout << fmt::format("(Mesh master loop: {})",timer_mesh.ElapsedTime()) << std::endl;

results_stream.close();

// assert(0);
}
