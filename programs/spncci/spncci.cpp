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
  University of Notre Dame

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
****************************************************************/

#include <cstdio>
#include <ctime>
#include <fstream>
#include <istream>
#include <iostream>
#include <sys/resource.h>
#include <omp.h>  

#include "SymEigsSolver.h"  // from spectra
#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "mcutils/io.h"
#include "spncci/recurrence.h"
#include "lgi/lgi_unit_tensors.h"
#include "mcutils/profiling.h"
#include "mcutils/eigen.h"
#include "spncci/branching.h"
#include "spncci/branching2.h"
#include "spncci/branching_u3s.h"
#include "spncci/branching_u3lsj.h"
#include "spncci/computation_control.h"
#include "spncci/decomposition.h"
#include "spncci/eigenproblem.h"
#include "spncci/explicit_construction.h"
#include "spncci/io_control.h"
#include "spncci/computation_control.h"
#include "spncci/parameters.h"
#include "spncci/results_output.h"
#include "spncci/transform_basis.h"

////////////////////////////////////////////////////////////////
// WIP code
//
// to extract to spncci library when ready
//////////////////////////////////////////////////////////////// 
namespace spncci
{}// end namespace

////////////////////////////////////////////////////////////////
// main body
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  std::cout<<"entering spncci"<<std::endl;
  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////
  bool check_unit_tensors=false;

  // SU(3) caching
  u3::U3CoefInit();
  u3::g_u_cache_enabled = true;

  // parameters for certain calculations
  spncci::g_zero_tolerance = 1e-6;
  spncci::g_suppress_zero_sectors = true;


  // Default binary mode, unless environment variable SPNCCI_RME_MODE
  // set to "text".
  //
  // This is meant as an ad hoc interface until text mode i/o is abolished.
  lsu3shell::g_rme_binary_format = true;
  char* spncci_rme_mode_cstr = std::getenv("SPNCCI_RME_MODE");
  if (spncci_rme_mode_cstr!=NULL)
    {
      const std::string spncci_rme_mode = std::getenv("SPNCCI_RME_MODE");
      if (spncci_rme_mode=="text")
        lsu3shell::g_rme_binary_format = false;
    }

  // run parameters
  std::cout << "Reading control file..." << std::endl;
  spncci::RunParameters run_parameters;

  // Eigen OpenMP multithreading mode
  Eigen::initParallel();
  // Eigen::setNbThreads(0);

  // open output files
  std::ofstream results_stream("spncci.res");

  // results output: code information
  spncci::StartNewSection(results_stream,"CODE");
  spncci::WriteCodeInformation(results_stream,run_parameters);

  // results output: run parameters
  spncci::StartNewSection(results_stream,"PARAMETERS");
  spncci::WriteRunParameters(results_stream,run_parameters);

  std::cout<<"Nmax="<<run_parameters.Nmax<<std::endl;

  // /////////////////////////////////////////////////////////////////////////////////////
  // // set up SpNCCI space
  // ////////////////////////////////////////////////////////////////
  bool restrict_sp3r_to_u3_branching=false;
    if(run_parameters.A<6)
      restrict_sp3r_to_u3_branching=true;

  lgi::MultiplicityTaggedLGIVector lgi_families;
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;
  spncci::BabySpNCCISpace baby_spncci_space;
  spncci::SpaceSpU3S spu3s_space;
  spncci::SpaceSpLS spls_space;
  spncci::SpaceSpJ spj_space;
  
  //Read in lgi families and generate spaces at different branching levels 
  //Nlimit allows for different irreps to be truncated to different Nmax
  int Nlimit=run_parameters.Nmax;
  spncci::SetUpSpNCCISpaces(
      run_parameters,lgi_families,spncci_space,sigma_irrep_map,
      baby_spncci_space,spu3s_space,spls_space,spj_space,
      results_stream,Nlimit,restrict_sp3r_to_u3_branching
    );

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

  // for(auto tensor :unit_tensor_labels)
  //   std::cout<<tensor.Str()<<std::endl;

  // generate unit tensor subspaces 
  u3shell::RelativeUnitTensorSpaceU3S 
    unit_tensor_space(run_parameters.Nmax,run_parameters.N1v,unit_tensor_labels);

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
  // set up indexing for branching
  ////////////////////////////////////////////////////////////////
  std::cout << "Set up basis indexing for branching..." << std::endl;
  //////////////////////////////////////////////////////////////////
  // NEW BRANCHING 
  //TODO: MAKE Vector with indices corresponding to run_parameters.J_values
  // vector be comes space with subspaces SpaceSpBasis in sectorsJ etcs. 
  std::vector<spncci::SpaceSpBasis> spaces_spbasis(run_parameters.J_values.size());
  for(int j=0; j<run_parameters.J_values.size(); ++j)
    {
      const HalfInt& J=run_parameters.J_values[j];
      spaces_spbasis[j]=spncci::SpaceSpBasis(baby_spncci_space,J);
    }
  //////////////////////////////////////////////////////////////////


//Will eventually remove.  For now just taking out of scope.
{ 
  ////////////////////////////////////////////////////////////////
  // Enumerate U3S sectors for observables 
  ////////////////////////////////////////////////////////////////
  std::cout << "Enumerating u3s sectors..." << std::endl;

  // enumerate u3S space from baby spncci for each observable 
  spncci::SpaceU3S space_u3s(baby_spncci_space);


  // Generate vector of hypersectors for each observable
  std::vector<spncci::ObservableHypersectorsU3S> 
    observable_hypersectors_by_observable(run_parameters.num_observables);
  for(int ob_num=0; ob_num<run_parameters.num_observables; ++ob_num)
    observable_hypersectors_by_observable[ob_num]=spncci::ObservableHypersectorsU3S(space_u3s,observable_spaces[ob_num]);

  // Write observable u3s hypersector information to results file
  spncci::WriteU3SHypersectorSectorInformation(
      results_stream,space_u3s,run_parameters.num_observables, 
      observable_hypersectors_by_observable
    );

  // set up basis indexing for branching
  std::map<HalfInt,spncci::SpaceLS> spaces_lsj;  // map: J -> space
  std::map<HalfInt,spncci::SpaceSpLS> spaces_splsj;
  for (const HalfInt J : run_parameters.J_values)
    {

      std::cout << fmt::format("Build LS space for J={}...",J.Str()) << std::endl;
      spaces_lsj[J] = spncci::SpaceLS(space_u3s,J);
      std::cout
        << fmt::format(
            "  subspaces {} dimension {}",
            J.Str(),
            spaces_lsj[J].size(),spaces_lsj[J].Dimension()
          ) << std::endl;

      // comparison tests with new basis branching construction
      std::cout << fmt::format("Build SpLS space for J={}...",J.Str()) << std::endl;
      spaces_splsj[J]=spncci::SpaceSpLS(spu3s_space,J);
      const auto& spls_space=spaces_splsj.at(J);
      // spncci::SpaceSpLS spls_space(spu3s_space,J);
      std::cout
        << fmt::format("  subspaces {} dimension {} full_dimension {}",
                       spls_space.size(),spls_space.Dimension(),spls_space.FullDimension()
          )
        << std::endl;
      std::cout
        << fmt::format("  compare... TotalDimensionU3LSJConstrained {}",TotalDimensionU3LSJConstrained(spncci_space,J))
        << std::endl;

      //////////////////////////////////////////////////////////////////
    }

  // determine J sectors for each observable
  std::vector<spncci::SectorsSpJ> observable_sectors;
  observable_sectors.resize(run_parameters.num_observables);

  for (int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
    {
      const int J0=run_parameters.observable_J0_values[observable_index];
      observable_sectors[observable_index] = spncci::SectorsSpJ(spj_space,J0);
    }


}


  ////////////////////////////////////////////////////////////////
  // terminate counting only run
  ////////////////////////////////////////////////////////////////
  // We now have to do all termination manually.  But, when the
  // control code is properly refactored, we can just have a single
  // termination, and the rest of the run can be in an "if
  // (!count_only)"...

  if (run_parameters.count_only)
    {

      // termination
      results_stream.close();

      std::cout << "End of counting-only run" << std::endl;
      std::exit(EXIT_SUCCESS);
    }






  ////////////////////////////////////////////////////////////////
  // precompute K matrices
  ////////////////////////////////////////////////////////////////
  std::cout << "Precompute K matrices..." << std::endl;

  // timing start
  mcutils::SteadyTimer timer_k_matrices;
  timer_k_matrices.Start();

  // traverse distinct sigma values in SpNCCI space, generating K
  // matrices for each
  // spncci::KMatrixCache k_matrix_cache;
  spncci::KMatrixCache k_matrix_cache, kinv_matrix_cache;
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

  ///////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"setting up lgi unit tensor blocks"<<std::endl;
  // Get list of unit tensor labels between lgi's 
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
      run_parameters.Nsigmamax, run_parameters.N1v,
      lgi_unit_tensor_labels,J0_for_unit_tensors,T0_for_unit_tensors,
      restrict_positive_N0
    );

  // //FOR TESTING
  // // explicit construction of spncci basis
  // basis::OperatorBlocks<double> spncci_expansions;
  // if(check_unit_tensors)
  //   spncci::ExplicitBasisConstruction(
  //     run_parameters,spncci_space,baby_spncci_space,
  //     k_matrix_cache, kinv_matrix_cache,
  //     restrict_sp3r_to_u3_branching,spncci_expansions
  //     );

  //Get look-up table for lgi index in full space.  Used for looking up seed filenames 
  // which are index by full space index 
  std::cout<<"reading lgi table "<<std::endl;
  std::vector<int> lgi_full_space_index_lookup;
  lgi::ReadLGILookUpTable(lgi_full_space_index_lookup,lgi_families.size());
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
  std::cout<<"Starting recurrence and contraction"<<std::endl;
    
  std::vector<spncci::LGIPair> lgi_pairs;
  spncci::GetLGIPairsForRecurrence(lgi_families,spncci_space,sigma_irrep_map,lgi_pairs);


  spncci::ObservableHypersectorsByLGIPairTable 
    observable_hypersectors_mesh(run_parameters.num_observables);

  // // by observable, by hw, by lgi pair
  // spncci::ObservableHyperblocksByLGIPairTable observable_hyperblocks_mesh(run_parameters.num_observables);

  // // Presize table
  // for(int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
  //   observable_hyperblocks_mesh[observable_index].resize(run_parameters.hw_values.size());

  //TODO: If doing change of basis for irrep families, read in transformation matrices 
  spncci::OperatorBlocks lgi_transformations;
    
  if(run_parameters.transform_lgi)
    {
      std::cout<<"reading in lgi transformations"<<std::endl;
      std::string lgi_transformations_filename="lgi_transformations.dat";
      spncci::ReadTransformationMatrices(lgi_transformations_filename,lgi_transformations);
      // for(int i=0; i<lgi_families.size(); ++i)
      //   {
      //     std::cout<<"---------------------------------------"<<std::endl;
      //     std::cout<<"irrep family "<<i<<std::endl;
      //     int j=lgi_full_space_index_lookup[i];  
      //     std::cout<<"full space index "<<j<<std::endl;        
      //     std::cout<<lgi_transformations[j]<<std::endl<<std::endl;
      //     std::cout<<"---------------------------------------"<<std::endl<<std::endl;
      //   }
    }


  // assert(0);


  std::cout<<"begin parallel region"<<std::endl;
  int num_files;
  // std::vector<std::vector<int>> num_lgi_pairs_per_thread(run_parameters.num_observables);
  std::vector<int> num_lgi_pairs_per_thread;
  #pragma omp parallel shared(observable_hypersectors_mesh,num_files,num_lgi_pairs_per_thread)
  // observable_hyperblocks_mesh,
    {
      
      #pragma omp single
      {
        int num_threads=omp_get_num_threads();

        if(num_threads>lgi_pairs.size())
          {
            std::cout<<"Too many threads.  Only "<<lgi_pairs.size()<<" needed."<<std::endl;
              assert(num_threads<=lgi_pairs.size());
          }

        num_lgi_pairs_per_thread.resize(num_threads);
        num_files=num_threads;
        
        // for(int i=0; i<run_parameters.num_observables; ++i)
        //   num_lgi_pairs_per_thread[i].resize(num_threads);
      }

      
      
      // Private containers 
      //
      //coefficient caches
      u3::UCoefCache u_coef_cache;
      u3::PhiCoefCache phi_coef_cache;

      mcutils::SteadyTimer timer_recurrence;
      timer_recurrence.Start();

      #pragma omp for schedule(dynamic) nowait
      // for(int i=0; i<12; ++i)
      for(int i=0; i<lgi_pairs.size(); ++i)
        {

          // Brought hyperbocks inside for loop to ensure deallocation 
          // by observable, by hw, by lgi pair
          spncci::ObservableHyperblocksTable observable_hyperblocks_table(run_parameters.num_observables);
          // Formerly:
          //    spncci::ObservableHyperblocksByLGIPairTable observable_hyperblocks_by_lgi_table(run_parameters.num_observables);
          // Presize table
          for(int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
            observable_hyperblocks_table[observable_index].resize(run_parameters.hw_values.size());

          const spncci::LGIPair& lgi_pair=lgi_pairs[i];
          int irrep_family_index_bra,irrep_family_index_ket;
          std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
          // std::cout<<irrep_family_index_bra<<"  "<<irrep_family_index_ket<<std::endl;
          // (irrep1,irrep2)
          // spncci::LGIPair lgi_pair(irrep_family_index_bra,irrep_family_index_ket);
          
          basis::OperatorHyperblocks<double> unit_tensor_hyperblocks;
          spncci::BabySpNCCIHypersectors baby_spncci_hypersectors;
          
          // Generate Unit tensor blocks if lgi pair seed files found.
          // If files not found, function returns false.
          bool files_found
            =GenerateUnitTensorHyperblocks(
                lgi_pair,run_parameters.Nmax, run_parameters.N1v,
                lgi_families,lgi_full_space_index_lookup,
                spncci_space,baby_spncci_space,unit_tensor_space,
                k_matrix_cache,kinv_matrix_cache,lgi_transformations,
                run_parameters.transform_lgi,u_coef_cache,phi_coef_cache,
                baby_spncci_hypersectors,unit_tensor_hyperblocks
              );

          if(not files_found)
            continue;

          // Check if hypersectors are diagonal in irrep family. 
          bool is_diagonal=irrep_family_index_ket==irrep_family_index_bra;

          // Initialize hyperblocks for (irrep2,irrep1)
          spncci::LGIPair lgi_pair2(irrep_family_index_ket,irrep_family_index_bra);
          basis::OperatorHyperblocks<double> unit_tensor_hyperblocks2;
          spncci::BabySpNCCIHypersectors baby_spncci_hypersectors2;
          
          if(not is_diagonal)
            {  
              // std::cout<<"conjugate pair"<<std::endl;
              bool files_found2
              =GenerateUnitTensorHyperblocks(
                  lgi_pair2,run_parameters.Nmax, run_parameters.N1v,
                  lgi_families,lgi_full_space_index_lookup,
                  spncci_space,baby_spncci_space,unit_tensor_space,
                  k_matrix_cache,kinv_matrix_cache,lgi_transformations,
                  run_parameters.transform_lgi,u_coef_cache,phi_coef_cache,
                  baby_spncci_hypersectors2,unit_tensor_hyperblocks2
                );

              // If we've made it this far (passed files_found) then files for (irrep2,irrep1) should exist
              assert(files_found2);
            }
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // Contract and regroup
          //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          // std::cout<<"contract "<<std::endl;
          spncci::ObservableHypersectorsTable observable_hypersectors_table(run_parameters.num_observables);
          spncci::ContractBabySpNCCIHypersectors(
            lgi_pair,run_parameters.num_observables, run_parameters.hw_values.size(),
            baby_spncci_space,observable_spaces,unit_tensor_space,
            baby_spncci_hypersectors,baby_spncci_hypersectors2,
            unit_tensor_hyperblocks,unit_tensor_hyperblocks2,
            observables_relative_rmes,observable_hypersectors_table,
            // observable_hypersectors_by_lgi_table,
            observable_hyperblocks_table
          );
          

          // spncci::WriteBabySpncciObservableRMEs(lgi_pair,observable_hyperblocks_table);
          spncci::WriteBabySpncciObservableRMEs(
            lgi_pair,observable_hypersectors_table,
            observable_hyperblocks_table
          );

          num_lgi_pairs_per_thread[omp_get_thread_num()]++;
          
          int lgi1, lgi2;
          std::tie(lgi1,lgi2)=lgi_pair;
          std::cout<<fmt::format("finished lgi pair {}  {}",lgi1,lgi2)<<std::endl;

        }// end lgi_pair

        timer_recurrence.Stop();

      //After lgi_pair for loop complete, each thread dealocates
      //  phi and U coefficient caches 
      u_coef_cache.clear();
      phi_coef_cache.clear();
    
    } //end parallel region
 



  ////////////////////////////////////////////////////////////////
  // calculation mesh master loop
  ////////////////////////////////////////////////////////////////

  std::cout << "Calculation mesh master loop..." << std::endl;

  // timing start
  mcutils::SteadyTimer timer_mesh;
  timer_mesh.Start();

  // Set num threads to one
  // omp_set_num_threads(1); 


  // W coefficient cache -- needed for observable branching
  u3::WCoefCache w_cache;

  // for each hw value, solve eigen problem and get expectation values 
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
            //////////////////////////////////////////////////////////////////
            // NEW BRANCHING
            //////////////////////////////////////////////////////////////////
            HalfInt J00 = run_parameters.observable_J0_values[observable_index];
            
            const spncci::SpaceSpBasis& spbasis_bra=spaces_spbasis[subspace_index];
            const spncci::SpaceSpBasis& spbasis_ket=spaces_spbasis[subspace_index];
            // std::cout<<" spbasis "<<J<<std::endl;
            // std::cout<<spbasis_bra.DebugStr(true)<<std::endl;

            const u3shell::ObservableSpaceU3S& observable_space=observable_spaces[observable_index];
            // std::cout<<"constructing "<<std::endl;
            spncci::OperatorBlock hamiltonian_matrix;
            spncci::ConstructOperatorMatrix(
              baby_spncci_space,
              observable_space,
              J00,
              // w_cache,
              spbasis_bra, 
              spbasis_ket, //For a given J
              num_lgi_pairs_per_thread,
              observable_index, hw_index,
              hamiltonian_matrix
            );

            spncci::WriteMatrixToFile(hamiltonian_matrix, hw);
            // std::cout<<hamiltonian_matrix<<std::endl;
            // long int num_nonzero_rmes=0;
            // for(int i=0; i<hamiltonian_matrix.rows(); ++i)
            //   for(int j=0; j<=i; ++j)
            //     {
            //       if(fabs(hamiltonian_matrix(i,j))>10e-4)
            //         num_nonzero_rmes++;
            //     }
            // std::cout<<"number of non-zero rmes "<<num_nonzero_rmes<<std::endl;

            std::cout << fmt::format("  Diagonalizing: J={}",J) << std::endl;
            spncci::SolveHamiltonian(
                hamiltonian_matrix,
                run_parameters.num_eigenvalues,
                run_parameters.eigensolver_num_convergence,  // whatever exactly this is...
                run_parameters.eigensolver_max_iterations,
                run_parameters.eigensolver_tolerance,
                eigenvalues_J,eigenvectors_J
              );

            //////////////////////////////////////////////////////////////////
          }

        // results output: eigenvalues
        spncci::WriteEigenvalues(results_stream,run_parameters.J_values,eigenvalues,run_parameters.gex);
        // spncci::WriteEigenvalues(results_stream,spj_space,eigenvalues,run_parameters.gex);

        ////////////////////////////////////////////////////////////////
        // do decompositions
        ////////////////////////////////////////////////////////////////

        std::cout << "Calculate eigenstate decompositions..." << std::endl;
        mcutils::SteadyTimer timer_decompositions;
        timer_decompositions.Start();

        // decomposition matrices:
        //   - vector over J 
        //   - matrix over (basis_subspace_index,eigenstate_index)
        //
        // That is, decompositions are stored as column vectors, within a
        // matrix, much like the eigenstates themselves.
        std::vector<spncci::Matrix> Nex_decompositions;
        std::vector<spncci::Matrix> baby_spncci_decompositions;
        Nex_decompositions.resize(run_parameters.J_values.size());
        baby_spncci_decompositions.resize(run_parameters.J_values.size());

        // calculate decompositions
        spncci::CalculateNexDecompositions(
          spaces_spbasis,eigenvectors,Nex_decompositions,
          run_parameters.Nsigma0,run_parameters.Nmax
        );
  
        spncci::CalculateBabySpNCCIDecompositions(
          spaces_spbasis,eigenvectors,baby_spncci_decompositions,
          baby_spncci_space.size()
        );

        timer_decompositions.Stop();
        std::cout << fmt::format("  (Decompositions: {})",timer_decompositions.ElapsedTime()) << std::endl;

        // // results output: decompositions
        spncci::WriteDecompositions(
          results_stream,"Nex",".6f",spaces_spbasis,
          Nex_decompositions,run_parameters.gex
        );

        spncci::WriteDecompositions(
          results_stream,"BabySpNCCI",".4e",spaces_spbasis,
          baby_spncci_decompositions,run_parameters.gex
        );

        /////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        // Writing irrep family blocks to files for use in lgi basis transformation

        if(true) //TEMP While doing higher Nmax runs
        // if(not run_parameters.transform_lgi)
        {
          //TODO: Remove restriction to 3 and make input
          int num_eigenvalues=std::min(run_parameters.num_eigenvalues,3);
          std::cout<<"basis transformation "<<std::endl;
          int num_irrep_families=lgi_families.size();
          std::vector<std::vector<spncci::OperatorBlocks>> irrep_family_blocks;

          spncci::RegroupIntoIrrepFamilies(
            spaces_spbasis,num_irrep_families,num_eigenvalues,
            eigenvectors,irrep_family_blocks
          );

          std::string test_filename=fmt::format("irrep_family_blocks_{}",hw);
          spncci::WriteIrrepFamilyBlocks(
            run_parameters.J_values,  num_irrep_families,num_eigenvalues,
            lgi_full_space_index_lookup,irrep_family_blocks,test_filename
          );

        }
      }// End Hamiltonian section

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      {
      //////////////////////////////////////////////////////////////
      // calculate observable RMEs
      ////////////////////////////////////////////////////////////////

      std::cout << "Calculate observable results..." << std::endl;
      mcutils::SteadyTimer timer_observables;
      timer_observables.Start();

      // observable_results_matrices:
      //   - vector over observable_index
      //   - vector over sector_index
      //   - matrix over (bra_eigenstate_index,ket_eigenstate_index)
      std::vector<spncci::OperatorBlocks> observable_results_matrices;
      observable_results_matrices.resize(run_parameters.num_observables);

      //Get sectors 
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


      // std::cout<<"calculate observable results"<<std::endl;
      for (int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
        {
          HalfInt J0 = run_parameters.observable_J0_values[observable_index];
          auto&sectors=observable_sectors[observable_index];

          observable_results_matrices[observable_index].resize(sectors.size());

          // std::cout<<"Get corresponding observable space"<<std::endl;
          const u3shell::ObservableSpaceU3S& observable_space=observable_spaces[observable_index];

          for(int sector_index=0; sector_index<sectors.size(); ++sector_index)
            {
              int bra_index,ket_index;
              std::tie(bra_index,ket_index)=sectors[sector_index];

              const spncci::SpaceSpBasis& spbasis_bra=spaces_spbasis[bra_index];
              const spncci::SpaceSpBasis& spbasis_ket=spaces_spbasis[ket_index];

              const HalfInt bra_J=run_parameters.J_values[bra_index];
              const HalfInt ket_J=run_parameters.J_values[ket_index];

              // std::cout<<"constructing "<<std::endl;
              spncci::OperatorBlock observable_block;
              spncci::ConstructOperatorMatrix(
                baby_spncci_space,observable_space,J0,
                // w_cache,
                spbasis_bra, spbasis_ket,
                num_lgi_pairs_per_thread,observable_index, hw_index,observable_block
              );


              // std::cout<<"calculate observable results"<<std::endl;
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
      
      }//observable 

    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  }

// timing stop
timer_mesh.Stop();
std::cout << fmt::format("(Mesh master loop: {})",timer_mesh.ElapsedTime()) << std::endl;

results_stream.close();


}
