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
       INT 1.0 4 0 0 0 relative_observables/JISP16_Nmax20_hw20.0_rel.dat      // coef Jmax J0 T0 g0 interaction_filename

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
****************************************************************/

#include <cstdio>
#include <ctime>
#include <fstream>
#include <sys/resource.h>

#include "SymEigsSolver.h"  // from spectra
#include "cppformat/format.h"

#include "lgi/lgi_solver.h"
// to vett as moved into computation_control 
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"
#include "spncci/branching.h"
#include "spncci/computation_control.h"
#include "spncci/parameters.h"
// #include "spncci/results_output.h"


////////////////////////////////////////////////////////////////
// WIP code
//
// to extract to spncci library when ready
////////////////////////////////////////////////////////////////
  
namespace spncci
{
void PrintHypersectors(
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
  const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  )
{
  for(int hypersector_index=0; hypersector_index<baby_spncci_hypersectors.size(); ++hypersector_index)
  {
    const auto& hypersector=baby_spncci_hypersectors.GetHypersector(hypersector_index);
    
    int unit_tensor_subspace_index, ket_subspace_index,bra_subspace_index, rho0;
    std::tie(bra_subspace_index, ket_subspace_index,unit_tensor_subspace_index,rho0)=hypersector.Key();

    const auto& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
    const auto& bra_subspace=baby_spncci_space.GetSubspace(bra_subspace_index);
    const auto& ket_subspace=baby_spncci_space.GetSubspace(ket_subspace_index);

    std::cout<<"hypersector "<<hypersector_index<<" "<< bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()
    <<"  "<<unit_tensor_subspace.LabelStr()<<rho0<<std::endl;
    for(int i=0; i<unit_tensor_subspace.size(); ++i)
    {
      int T0,Sp,Tp,S,T;
      std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(i);
      std::cout<<fmt::format("{}  {} {}  {} {}",T0,Sp,Tp,S,T)<<std::endl;
      std::cout<<unit_tensor_hyperblocks[hypersector_index][i]<<std::endl<<std::endl;
    }
  }

}

void PrintU3SSector(
  const std::vector<double>& hw_values,
  const std::vector<std::vector<spncci::SectorLabelsU3S>>& observables_sectors_u3s,
  std::vector<std::vector<basis::OperatorBlocks<double>>>& observables_blocks_u3s, //can't be constant because of chop function
  const spncci::SpaceU3S& space_u3s,
  int num_observables
  )
// Prints out U3SSectors and blocks 
{
  for(int observable_index=0; observable_index<num_observables; ++observable_index)
    for(int h=0; h<hw_values.size(); ++h)
      {
        std::cout<<"observable "<<observable_index<<" hw "<<hw_values[h]<<std::endl;
        const std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];
        basis::OperatorBlocks<double>& blocks_u3s=observables_blocks_u3s[h][observable_index];
        for(int i=0; i<blocks_u3s.size(); ++i)
          {
            auto& block=blocks_u3s[i];
            const auto& sector=sectors_u3s[i];            
            const auto& bra=space_u3s.GetSubspace(sector.bra_index());
            const auto& ket=space_u3s.GetSubspace(sector.ket_index());
            const auto& op=sector.operator_labels();
            if(not mcutils::IsZero(block))
            {
              std::cout<<"block number "<<i<<std::endl;
              std::cout<<bra.Str()<<"  "<<ket.Str()<<"  "<<op.Str()<<"  "<<sector.kappa0()<<"  "<<sector.L0()<<"  "<<sector.rho0()<<std::endl;
              mcutils::ChopMatrix(block, 1e-6);
              std::cout<<block<<std::endl<<std::endl;
            }
          }
      }
}

  void
  WriteEigenValues(
    const std::vector<HalfInt>& J_values, double hw, 
    int Nmax, int Nsigma0_ex_max,
    std::map<HalfInt,Eigen::VectorXd>& eigenvalues,
    std::vector<std::string>& observable_filenames,
    std::vector<int>& scalar_observable_indices,
    std::vector<std::map<HalfInt,Eigen::VectorXd>>& scalar_observable_expectations,
    std::vector<int>& nonscalar_observable_indices,
    std::vector<std::map<spncci::JPair,Eigen::MatrixXd>>& nonscalar_observable_expectations
  )
  // for observables with J0=0, line them up with energy eigenvalue and read off diagonal matrix elements 
  // for observables with J0!=0, then have their own section --probably do this in the code as well, i.e.,
  {
    std::string filename=fmt::format("eigenvalues_Nmax{:02d}_Nsigma_ex{:02d}.dat",Nmax,Nsigma0_ex_max);
    std::cout<<"writing to file"<<std::endl;
    std::fstream fs;
    const int width=3;
    const int precision=16;
    fs.open (filename, std::fstream::out | std::fstream::app);
    fs << std::setprecision(precision);

    fs << "OUPTPUT from spncci Version 1"<<std::endl<<std::endl;;
    fs << "Scalar observables:";
    for(int i=0; i<scalar_observable_indices.size(); ++i)
      fs <<"  "<<observable_filenames[scalar_observable_indices[i]];
    fs << std::endl;

    fs <<"Nonscalar observables:";
    for(int i=0; i<nonscalar_observable_indices.size(); ++i)
      fs <<"  "<<observable_filenames[nonscalar_observable_indices[i]];

    fs << std::endl<<fmt::format("hw {:2.1f}", hw)<<std::endl;

    for(HalfInt J : J_values)
      {
        Eigen::VectorXd& eigenvalues_J=eigenvalues[J];
        // Eigen::VectorXd& observables=observable_expectations[J];

        for(int i=0; i<eigenvalues_J.size(); ++i)
          {
            double eigenvalue=eigenvalues_J(i);
            std::cout<<fmt::format("{:2d}   {}   {:8.5f}",i, J,eigenvalue);
            fs << fmt::format("{:2d}   {}   {:8.5f}",i, J,eigenvalue);
            
            for(int j=0; j<scalar_observable_indices.size(); ++j)  
            {          
              std::cout<<fmt::format("   {:8.5f}",scalar_observable_expectations[j][J](i))
              <<std::endl<<std::endl;
              fs <<fmt::format("   {:8.5f}",scalar_observable_expectations[j][J](i))
              <<std::endl<<std::endl;
            }
          }
      }
    for(HalfInt J : J_values)
      {
        Eigen::VectorXd& eigenvalues_J_initial=eigenvalues[J];
        for(int i=0; i<eigenvalues_J_initial.size(); ++i)
          for(HalfInt Jp :J_values)
          {
            Eigen::VectorXd& eigenvalues_J_final=eigenvalues[Jp];
            spncci::JPair J_pair(Jp,J);
            for(int ip=0; ip<eigenvalues_J_final.size(); ++ip)
              {
                double eigenvalue_initial=eigenvalues_J_initial(i);
                double eigenvalue_final=eigenvalues_J_final(ip);
                fs << fmt::format("{:2d}   {}   {:2d}   {}   {:8.5f}   {:8.5f}",
                      i,J,ip,Jp,eigenvalue_initial,eigenvalue_final);

                for(int j=0; j<nonscalar_observable_indices.size(); ++j)
                  {
                    // std::cout<<"(ip, i) ("<<ip<<","<<i<<")"<<std::endl;
                    Eigen::MatrixXd& obserable_matrix=nonscalar_observable_expectations[j][J_pair];
                    // std::cout<<obserable_matrix<<std::endl;
                    double observable=obserable_matrix(ip,i);
                    fs<<fmt::format("  {:8.5f}",observable);
                  }
                fs<<std::endl;
              }
          }
      }
    fs<<std::endl<<std::endl;
    fs.close();
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
  
  // SU(3) caching
  u3::U3CoefInit();
  u3::UCoefCache u_coef_cache;
  u3::PhiCoefCache phi_coef_cache;
  u3::g_u_cache_enabled = true;

  // parameters for certain calculations
  spncci::g_zero_tolerance = 1e-6;
  spncci::g_suppress_zero_sectors = true;

  // rme input mode
  //
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
  spncci::RunParameters run_parameters(argc,argv);

  // Eigen OpenMP multithreading mode
  Eigen::initParallel();
  // Eigen::setNbThreads(0);

  // open results output file
   //std::ofstream results_stream("spncci.res");
   //spncci::WriteResultsFileHeader(results_stream);
   //results_stream.close();

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis
  ////////////////////////////////////////////////////////////////

  std::cout << "Read lsu3shell basis..." << std::endl;
  // read lsu3shell basis (regroup into U3SPN subspaces)
  lsu3shell::LSU3BasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3Basis(
      run_parameters.Nsigma_0,run_parameters.lsu3shell_basis_filename,
      lsu3shell_basis_table,lsu3shell_basis_provenance,lsu3shell_space
    );

  ////////////////////////////////////////////////////////////////
  // solve for LGIs
  ////////////////////////////////////////////////////////////////
  std::cout << "Solve for LGIs..." << std::endl;

  // timing start
  Timer timer_lgi;
  timer_lgi.Start();

  lgi::MultiplicityTaggedLGIVector lgi_families;
  basis::MatrixVector lgi_expansions;
  
  spncci::GetLGIExpansion(
      lsu3shell_space,lsu3shell_basis_table,
      run_parameters.Brel_filename,run_parameters.Nrel_filename,
      run_parameters.A, run_parameters.Nsigma_0,
      lgi_families, lgi_expansions
    );

  // diagnostics
  std::cout << fmt::format("  LGI families {}",lgi_families.size()) << std::endl;
  if (false)
    lgi::WriteLGILabels(lgi_families,std::cout);

  // timing stop
  timer_lgi.Stop();
  std::cout << fmt::format("(Task time: {})",timer_lgi.ElapsedTime()) << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////
  // set up SpNCCI space
  ////////////////////////////////////////////////////////////////

  std::cout << "Set up SpNCCI space..." << std::endl;

  // build SpNCCI irrep branchings
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // persistent container to store branchings
  spncci::NmaxTruncator truncator(run_parameters.Nsigma_0,run_parameters.Nmax);
  spncci::GenerateSpNCCISpace(lgi_families,truncator,spncci_space,sigma_irrep_map);

  for(int i=0; i<spncci_space.size(); ++i)
    std::cout<<i<<"  "<<spncci_space[i].Str()<<spncci_space[i].gamma_max()<<std::endl;

  // diagnostics
  std::cout << fmt::format("  Irrep families {}",spncci_space.size()) << std::endl;
  std::cout << fmt::format("  TotalU3Subspaces {}",spncci::TotalU3Subspaces(spncci_space)) << std::endl;
  std::cout << fmt::format("  TotalDimensionU3 {}",spncci::TotalDimensionU3S(spncci_space)) << std::endl;

  // build baby spncci space 
  spncci::BabySpNCCISpace baby_spncci_space(spncci_space);

  // build SpU3S gathered space
  std::cout << "Build SpU3S space..." << std::endl;
  spncci::SpaceSpU3S spu3s_space(baby_spncci_space);
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
  spncci::SpaceSpLS spls_space(spu3s_space);
  std::cout
    << fmt::format("  subspaces {} dimension {} full_dimension {}",
                   spls_space.size(),spls_space.Dimension(),spls_space.FullDimension()
      )
    << std::endl;
  std::cout
    << fmt::format("  compare... TotalDimensionU3LS {}",TotalDimensionU3LS(spncci_space))
    << std::endl;
  // std::cout << splss_space.DebugStr(true);


  ////////////////////////////////////////////////////////////////
  // precompute K matrices
  ////////////////////////////////////////////////////////////////
  std::cout << "Precompute K matrices..." << std::endl;

  // timing start
  Timer timer_k_matrices;
  timer_k_matrices.Start();

  // traverse distinct sigma values in SpNCCI space, generating K
  // matrices for each
  spncci::KMatrixCache k_matrix_cache;
  bool intrinsic = true;
  spncci::PrecomputeKMatrices(sigma_irrep_map,k_matrix_cache,intrinsic);

  // timing stop
  timer_k_matrices.Stop();
  std::cout << fmt::format("(Task time: {})",timer_k_matrices.ElapsedTime()) << std::endl;

  // std::cout<<"Kmatrices "<<std::endl;
  // for(auto it=k_matrix_cache.begin(); it!=k_matrix_cache.end(); ++it)
  //   {
  //     std::cout<<"sigma "<<it->first.Str()<<std::endl;
  //     for(auto it2=it->second.begin();  it2!=it->second.end(); ++it2)
  //     {
  //       std::cout<<"  omega"<<it2->first.Str()<<std::endl;
  //       std::cout<<it2->second<<std::endl;
  //     }
  //   }

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

  // ///////////////////////////////////////////////////////////////////////////////////////////////////
  // //  For testing, get lsu3shell expansion of full spncci basis
  // ///////////////////////////////////////////////////////////////////////////////////////////////////        
  /////////////////////////////////////////////////////////////////////////////////////
  // // For explicit construction 
  // u3shell::SectorsU3SPN Aintr_sectors;
  // basis::MatrixVector Aintr_matrices;
  // spncci::ReadLSU3ShellSymplecticRaisingOperatorRMEs(
  //     lsu3shell_basis_table,lsu3shell_space, 
  //     run_parameters.Arel_filename,Aintr_sectors,Aintr_matrices,
  //     run_parameters.A
  //   );

  // if(run_parameters.Nmax==run_parameters.Nsigma0_ex_max)
  // {
  // basis::MatrixVector spncci_expansions;
  // spncci::ConstructSpNCCIBasisExplicit(
  //     lsu3shell_space,spncci_space,lgi_expansions,baby_spncci_space,
  //     k_matrix_cache,Aintr_sectors,Aintr_matrices,spncci_expansions
  //   );
  // }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //  Read in observables  
  ///////////////////////////////////////////////////////////////////////////////////////////////////        
  // Initialize containers for rmes and their symmetries 
  // Stored by hw, then by observable
  std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>> observables_relative_rmes(run_parameters.hw_values.size());
  std::vector<std::vector<u3shell::IndexedOperatorLabelsU3S>> observable_symmetries_u3s(run_parameters.num_observables); 

  spncci::ReadRelativeObservables(
    run_parameters.Nmax, run_parameters.N1v, run_parameters.hw_values,
    run_parameters.observable_directory,run_parameters.observable_filenames, 
    unit_tensor_space, observables_relative_rmes, observable_symmetries_u3s
  );

  // // set up U3S sectors for each of the observables 
  // // Was the function ConstructObservablesU3S
  ///////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"setting up u3 sectors "<<std::endl;
  // enumerate u3S space from baby spncci for each observable 
  spncci::SpaceU3S space_u3s(baby_spncci_space);

  // vector of sectors for each observable
  std::vector<std::vector<spncci::SectorLabelsU3S>> observables_sectors_u3s(run_parameters.num_observables);
  
  // vector of blocks for u3 sectors for each hbar omega,for each observable
  std::vector<std::vector<basis::OperatorBlocks<double>>> observables_blocks_u3s(run_parameters.hw_values.size());

  // for each observable, enumerate sectors 
  for(int observable_index=0; observable_index< run_parameters.num_observables; ++observable_index) 
    {
      std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];
      spncci::GetSectorsU3S(space_u3s,observable_symmetries_u3s[observable_index],sectors_u3s);
    }

  // For each hbar omega, zero initialize block for each observable
  // based on basis::SetOperatorToZero in operator.h
  for(int h=0; h<run_parameters.hw_values.size(); ++h)
    {
      std::vector<basis::OperatorBlocks<double>>& observables_blocks=observables_blocks_u3s[h];
      observables_blocks.resize(run_parameters.num_observables);

      for(int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
        {
          basis::OperatorBlocks<double>& blocks=observables_blocks[observable_index];
          std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];
          blocks.resize(sectors_u3s.size());
          for(int sector_index=0; sector_index<sectors_u3s.size(); ++sector_index)
            {
              int rows=space_u3s.GetSubspace(sectors_u3s[sector_index].bra_index()).full_dimension();
              int cols=space_u3s.GetSubspace(sectors_u3s[sector_index].ket_index()).full_dimension();
              blocks[sector_index]=basis::OperatorBlock<double>::Zero(rows,cols);
            }
        }

    } 
  ///////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"seting up lgi unit tensor blocks"<<std::endl;
  
  // map of {lgi pair : list of hypersector indices organized by Nsum}
  // Read in lsu3shell unit tensors
  // transform block for each unit tensor to spncci
  // identify unit tensors with non-zero rmes's between each lgi pair 
  // Generate unit tensor labels for recurrence for each lgi pair
  // put seed blocks into hypersector blocks for each lgi pair 

  // Get list of unit tensor labels between lgi's 
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
    run_parameters.Nsigma0_ex_max, run_parameters.N1v,
    lgi_unit_tensor_labels,J0_for_unit_tensors,T0_for_unit_tensors,
    restrict_positive_N0
    );

  //////////////////////////////////////////////////////////////////////////////////////////
  // for each unit tensor, read in unit tensor lsu3shell rmes and transform to spncci basis
  //////////////////////////////////////////////////////////////////////////////////////////
  // timing start
  Timer timer_read_seeds;
  timer_read_seeds.Start();

  std::cout << "Get seed unit tensor rmes..." << std::endl;
  // diagnostic
  std::cout << fmt::format("  seed unit tensors {}",lgi_unit_tensor_labels.size()) << std::endl;
  
    // Container for lgi unit tensor blocks 
  std::map< std::pair<int,int>, std::map<std::pair<int,int>, basis::OperatorBlocks<double>>> lgi_unit_tensor_blocks;
  
  // Get unit tensor seeds obtained from lsu3shell rmes transformed to spncci basis.
  spncci::GetUnitTensorSeedBlocks(
    lgi_unit_tensor_labels,unit_tensor_space,
    run_parameters.relative_unit_tensor_filename_template,
    lsu3shell_space, lsu3shell_basis_table,
    lgi_expansions, baby_spncci_space,
    lgi_unit_tensor_blocks
    );

  timer_read_seeds.Stop();
  std::cout << fmt::format("(Task time: {})",timer_read_seeds.ElapsedTime()) << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //create a map of vectors of unit tensor subspace indices keyed by spncci irrep pairs 
  std::map<std::pair<int,int>,std::set<int>>lgi_unit_tensor_subset;
  for(auto it=lgi_unit_tensor_blocks.begin(); it!=lgi_unit_tensor_blocks.end(); ++it)
    for(auto it2=it->second.begin(); it2!=it->second.end(); ++it2)
      lgi_unit_tensor_subset[it->first].insert(it2->first.first);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // for each lgi pair lgi1, lgi2 compute all unit tensor hypersectors for which Nnp<=Nn and 
  // the conjugate hypersectors for Nnp>Nn, i.e., compute <lgi1 Nnp=0 | |lgi2 Nn=2> and 
  // <lgi2 Nn=0| |lgi1 Nnp=2> etc. 

  // Need to add seeds for both lgi pair and conjugate lgi pair
  for(auto it=lgi_unit_tensor_blocks.begin(); it!=lgi_unit_tensor_blocks.end(); ++it)
    {
      int irrep_family_index_bra,irrep_family_index_ket;
      std::tie(irrep_family_index_bra,irrep_family_index_ket)=it->first;
      
      if(irrep_family_index_bra>irrep_family_index_ket)
        continue;      

      // Generate hypersectors for recurrence 
      std::vector<std::vector<int>> unit_tensor_hypersector_subsets;
      spncci::BabySpNCCIHypersectors baby_spncci_hypersectors;
      spncci::GenerateRecurrenceHypersectors(
        unit_tensor_space,baby_spncci_space,lgi_unit_tensor_subset,
        run_parameters.Nmax, irrep_family_index_bra, irrep_family_index_ket,
        unit_tensor_hypersector_subsets,
        baby_spncci_hypersectors
      );

      // std::cout<<"number of hypersectors "<<baby_spncci_hypersectors.size()<<std::endl;

      // // std::cout<<"checking hypersector subsets"<<std::endl;
      // for(int N=0; N<unit_tensor_hypersector_subsets.size(); N++)
      //   for(int hypersector_index : unit_tensor_hypersector_subsets[N])
      //     {
      //       // std::cout<<"N="<<N<<std::endl;
      //       const auto& hypersector=baby_spncci_hypersectors.GetHypersector(hypersector_index);
      //       int unit_tensor_subspace_index, ket_subspace_index,bra_subspace_index, rho0;
      //       std::tie(bra_subspace_index, ket_subspace_index,unit_tensor_subspace_index,rho0)=hypersector.Key();
    
      //       const auto& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
      //       const auto& bra_subspace=baby_spncci_space.GetSubspace(bra_subspace_index);
      //       const auto& ket_subspace=baby_spncci_space.GetSubspace(ket_subspace_index);

      //       // std::cout<<"hypersector "<<hypersector_index<<" "<< bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()
      //       // <<"  "<<unit_tensor_subspace.LabelStr()<<rho0<<std::endl;

      //     }

      // zero initialize hypersectors 
      basis::OperatorHyperblocks<double> unit_tensor_hyperblocks;
      basis::SetHyperoperatorToZero(baby_spncci_hypersectors,unit_tensor_hyperblocks);


      // For testing with explicit construction
      // TODO: extract into separate loop
      // if(run_parameters.Nmax==run_parameters.Nsigma0_ex_max)
      // {
      //   basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_explicit;
      // basis::SetHyperoperatorToZero(baby_spncci_hypersectors,unit_tensor_hyperblocks_explicit);
      // }
      ///////////////////////////////////////////////////////////////////////////////////////////////
      
      // Populate hypersectors with seeds
      spncci::PopulateHypersectorsWithSeeds(
        irrep_family_index_bra, irrep_family_index_ket,
        unit_tensor_hypersector_subsets[0],
        baby_spncci_space,baby_spncci_hypersectors, 
        lgi_unit_tensor_blocks,unit_tensor_hyperblocks
        );


      // Recurse over unit tensor hypersectors 
      // std::cout<<"entering the recurrence for "<<irrep_family_index_bra<<" "<<irrep_family_index_ket<<std::endl;
      spncci::ComputeUnitTensorHyperblocks(
        run_parameters.Nmax,run_parameters.N1v,u_coef_cache,phi_coef_cache,k_matrix_cache,
        spncci_space,baby_spncci_space,unit_tensor_space,
        baby_spncci_hypersectors, unit_tensor_hypersector_subsets,
        unit_tensor_hyperblocks
        );

      // std::cout<<"checking hypersectors"<<std::endl;

      // spncci::PrintHypersectors(
      //   baby_spncci_space,unit_tensor_space, 
      //   baby_spncci_hypersectors,unit_tensor_hyperblocks
      //   );

      // spncci::CheckUnitTensorRecurrence(
      //   irrep_family_index_bra, irrep_family_index_ket,
      //   unit_tensor_space,lgi_unit_tensor_labels,
      //   run_parameters.relative_unit_tensor_filename_template,
      //   lsu3shell_space, lsu3shell_basis_table,
      //   spncci_space,baby_spncci_space,spncci_expansions,
      //   baby_spncci_hypersectors,unit_tensor_hyperblocks
      // );
      
      // std::cout<<"contracting over observables "<<std::endl;
      for(int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
        for(int h=0; h<run_parameters.hw_values.size(); ++h)
          {
            // std::cout<<"observable "<<observable_index<<" hw "<<run_parameters.hw_values[h]<<std::endl;
            const u3shell::RelativeRMEsU3SSubspaces& relative_observable=observables_relative_rmes[h][observable_index];
            const std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];
            basis::OperatorBlocks<double>& blocks_u3s=observables_blocks_u3s[h][observable_index];
      
            ContractAndRegroupU3S(
                unit_tensor_space,baby_spncci_space,
                space_u3s,relative_observable,
                baby_spncci_hypersectors,unit_tensor_hyperblocks,
                sectors_u3s,blocks_u3s
              );      
          }
    }// end lgi_pair


    // spncci::PrintU3SSector(
    //   run_parameters.hw_values,
    //   observables_sectors_u3s,observables_blocks_u3s,  
    //   space_u3s, run_parameters.num_observables
    // );


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // At this point observable rmes should be fully computed and unit tensor cache, Ucoef cache and Kmatrix cache deleted 
  // Delete Kmatrix
  // Delete Unit tensor Cache
  // Delete Ucoef Cache 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // loop over hw values, branch matrix sectors and compute eigenvalues
  
  ////////////////////////////////////////////////////////////////
  // branch observables
  ////////////////////////////////////////////////////////////////
  std::cout << "Set up basis indexing for branching..." << std::endl;

  // set up basis indexing for branching
  std::map<HalfInt,spncci::SpaceLS> spaces_lsj;  // map: J -> space
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
      spncci::SpaceSpLS spls_space(spu3s_space,J);
      std::cout
        << fmt::format("  subspaces {} dimension {} full_dimension {}",
                       spls_space.size(),spls_space.Dimension(),spls_space.FullDimension()
          )
        << std::endl;
      std::cout
        << fmt::format("  compare... TotalDimensionU3LSJConstrained {}",TotalDimensionU3LSJConstrained(spncci_space,J))
        << std::endl;
    }



  std::cout << "Construct branched observable matrices..." << std::endl;

  // timing start
  Timer timer_branching;
  timer_branching.Start();

  // for each hw value, solve eigen problem and get expectation values 
  for(int h=0; h<run_parameters.hw_values.size(); ++h)
  {

    ////////////////////////////////////////////////////////////////
    // Formerly computational_control ConstructBranchedObservables
    ////////////////////////////////////////////////////////////////
    // populate fully-branched many-body matrices for observables

    // map: observable -> J ->  matrix
    std::vector<std::map<std::pair<HalfInt, HalfInt>,Eigen::MatrixXd>> observable_matrices;  
    observable_matrices.resize(run_parameters.num_observables);
    spncci::ConstructBranchedObservables(
      space_u3s,
      observables_sectors_u3s,
      observables_blocks_u3s[h], 
      spaces_lsj,
      run_parameters.num_observables,
      run_parameters.J_values,
      run_parameters.observable_Jvalues, 
      observable_matrices);

    // timing stop
    timer_branching.Stop();
    std::cout << fmt::format("(Task time: {})",timer_branching.ElapsedTime()) << std::endl;

    // diagnostics: branched matrices
    if (false)
      {
        for (int observable_index=0; observable_index<run_parameters.num_observables; ++observable_index)
          for (const HalfInt J : run_parameters.J_values)
            for (const HalfInt Jp : run_parameters.J_values)
              {
                const HalfInt J0=run_parameters.observable_Jvalues[observable_index];
                std::cout<<"Jvalues "<<Jp<<" "<<J0<<" "<<J<<std::endl;
                if(not am::AllowedTriangle(J,J0,Jp))
                  continue;

                const Eigen::MatrixXd& observable_matrix = observable_matrices[observable_index].at(std::make_pair(Jp,J));

                const HalfInt bra_J = Jp;
                const HalfInt ket_J = J;
                std::cout
                  << fmt::format("Observable {} bra_J {} ket_J {} J0 {}",observable_index,bra_J,ket_J, J0)
                  << std::endl;
                std::cout
                  << mcutils::FormatMatrix(observable_matrix,"8.5f")
                  << std::endl
                  << std::endl;

            }
      }
    
    std::map<HalfInt,Eigen::VectorXd> eigenvalues;  // map: J -> eigenvalues
    std::map<HalfInt,spncci::MatrixType> eigenvectors;  // map: J -> eigenvectors
    std::cout<<"solving Hamiltonian"<<std::endl;

    for (const HalfInt J : run_parameters.J_values)
    {
      // auto& matrix=observable_matrices[0].at(std::make_pair(J,J));
      // for(int i=0; i<matrix.rows(); ++i)
      //   for(int j=0; j<matrix.cols(); ++j)
      //     if(fabs(matrix(i,j)-matrix(j,i))>10e-6)
      //       std::cout<<i<<"  "<<j<<"  "<<matrix(i,j)<<"  "<<matrix(j,i)<<std::endl;

      spncci::SolveHamiltonian(observable_matrices[0].at(std::make_pair(J,J)),J,
        run_parameters.num_eigenvalues,
        run_parameters.eigensolver_num_convergence,  // whatever exactly this is...
        run_parameters.eigensolver_max_iterations,
        run_parameters.eigensolver_tolerance,
        eigenvalues, eigenvectors
        );
    }

  }



}
