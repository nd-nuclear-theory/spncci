/****************************************************************
  computation_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/computation_control.h"

#include "SymEigsSolver.h"  // from spectra
#include "cppformat/format.h"
#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"

namespace spncci
{

void GetLGIExpansion(
    const u3shell::SpaceU3SPN& lsu3shell_space, 
    const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
    const std::string& Brel_filename,
    const std::string& Nrel_filename,
    int A, HalfInt Nsigma_0,
    lgi::MultiplicityTaggedLGIVector& lgi_families,
    basis::MatrixVector& lgi_expansions
  )
  {
  u3shell::SectorsU3SPN Bintr_sectors, Nintr_sectors;
  basis::MatrixVector Bintr_matrices, Nintr_matrices;
  spncci::ReadLSU3ShellSymplecticOperatorRMEs(
      lsu3shell_basis_table,lsu3shell_space, 
      Brel_filename,Bintr_sectors,Bintr_matrices,
      Nrel_filename,Nintr_sectors,Nintr_matrices,
      A
    );

  const u3shell::SectorsU3SPN& Ncm_sectors = Nintr_sectors;
  basis::MatrixVector Ncm_matrices;
  lsu3shell::GenerateLSU3ShellNcmRMEs(
      lsu3shell_space,Nintr_sectors,Nintr_matrices,
      A,Ncm_matrices
    );

  lgi::GenerateLGIExpansion(
      lsu3shell_space, 
      Bintr_sectors,Bintr_matrices,Ncm_sectors,Ncm_matrices,
      Nsigma_0,lgi_families,lgi_expansions
    );
}



void GetUnitTensorSeedBlocks(
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const std::string& relative_unit_tensor_filename_template,
  const u3shell::SpaceU3SPN& lsu3shell_space, 
  const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
  const basis::MatrixVector& lgi_expansions,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  std::map< std::pair<int,int>, std::map<std::pair<int,int>, basis::OperatorBlocks<double>>>& lgi_unit_tensor_blocks
  )
{
  //If parallellize, iterate first over unit_tensor_index in parallel and then write to new map in serial?, maybe
  // Alternatively, we zero initalize lgi_unit_tensor_blocks
  // std::unordered_map<u3shell::SectorsU3SPN,basis::MatrixVector, boost::hash<u3shell::SectorsU3SPN>> unit_tensor_spncci_blocks;
  for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
    {
      // Get labels of corresponding unit tensor 
      const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];

      // File name containing lsu3shell rmes of unit tensor 
      std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);

      // generate unit tensor sectors in lsu3shell basis 
      const bool spin_scalar = false; // All spin values allowed
      u3shell::SectorsU3SPN unit_tensor_sectors(lsu3shell_space,unit_tensor_labels,spin_scalar);
      
      // read in lsu3shell rms of unit tensor 
      basis::MatrixVector unit_tensor_lsu3shell_matrices;
      lsu3shell::ReadLSU3ShellRMEs(
          filename,
          lsu3shell_basis_table,lsu3shell_space,
          unit_tensor_labels,unit_tensor_sectors,unit_tensor_lsu3shell_matrices
        );

      // transform seed rmes to SpNCCI basis (among LGIs)
      // basis::MatrixVector unit_tensor_spncci_matrices=unit_tensor_spncci_blocks[unit_tensor_sectors]; // temporary container if omp?
      basis::MatrixVector unit_tensor_spncci_matrices; // temporary
      lgi::TransformOperatorToSpBasis(
          unit_tensor_sectors,lgi_expansions,
          unit_tensor_lsu3shell_matrices,unit_tensor_spncci_matrices
        );

      // Extract unit tensor labels 
      u3::SU3 x0; 
      HalfInt S0,T0,Sp,Tp,S,T;
      int etap,eta;
      std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=unit_tensor_labels.FlatKey();

      // Look up unit tensor subspace
      u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
      int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
      auto& subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

      // Look up index of unit tensor in subspace
      int unit_tensor_state_index=subspace.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0), int(Sp),int(Tp),int(S),int(T)));


      //unit_tensor_spncci_matrices
      //baby_spncci_space
      //unit_tensor_sectors
      //lgi_unit_tensor_blocks

      // Turn into function, iterate over map, then over sector index


      // transfer spncci unit tensor rmes from temporary container to lgi hyperblocks container for use in spncci recurrence
      for(int sector_index=0; sector_index<unit_tensor_sectors.size(); ++sector_index)
          {
            // Check that the sector has non-zero rmes as defined by the zero_tolerance
            // If not, then don't need to copy over rmes.
            if(mcutils::IsZero(unit_tensor_spncci_matrices[sector_index],spncci::g_zero_tolerance))
              continue;
            
            // extract U3SPN sector information
            const typename u3shell::SectorsU3SPN::SectorType& sector = unit_tensor_sectors.GetSector(sector_index);
            const int bra_subspace_index = sector.bra_subspace_index();
            const int ket_subspace_index = sector.ket_subspace_index();
            const u3shell::SubspaceU3SPN bra_subspace=sector.bra_subspace();
            const u3shell::SubspaceU3SPN ket_subspace=sector.ket_subspace();
            const u3::U3& bra_sigma = bra_subspace.U3();
            const u3::U3& ket_sigma = ket_subspace.U3();
            const int rho0 = sector.multiplicity_index();
            // const int rho0 = unit_tensor_sectors.GetSector(sector_index).multiplicity_index();

            // Get baby spncci index 
            spncci::BabySpNCCISubspaceLabels 
              baby_spncci_bra(bra_sigma,bra_subspace.Sp(),bra_subspace.Sn(), bra_subspace.S(),bra_sigma);
            spncci::BabySpNCCISubspaceLabels 
              baby_spncci_ket(ket_sigma,ket_subspace.Sp(),ket_subspace.Sn(), ket_subspace.S(), ket_sigma);

            int baby_spncci_index_bra=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_bra);
            int baby_spncci_index_ket=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_ket);

            int irrep_family_index_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra).irrep_family_index();
            int irrep_family_index_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket).irrep_family_index();
            
            std::pair<int,int> irrep_family_pair(irrep_family_index_bra,irrep_family_index_ket);
            auto& subspace_blocks=lgi_unit_tensor_blocks[irrep_family_pair];
            
            // If multiplicity tagged unit tensor subspace not already initialized in container,
            //  need to resize container for given unit tensor subspace.
            std::pair<int,int> multiplicity_tagged_unit_tensor_subspace(unit_tensor_subspace_index,rho0);
            if(not subspace_blocks.count(multiplicity_tagged_unit_tensor_subspace))
              subspace_blocks[multiplicity_tagged_unit_tensor_subspace].resize(unit_tensor_space.GetSubspace(unit_tensor_subspace_index).size());

            // Transfer unit tensor rme blocks from temporary container to hypersector container
            subspace_blocks[multiplicity_tagged_unit_tensor_subspace][unit_tensor_state_index]=unit_tensor_spncci_matrices[sector_index];
           }
    }
}


void
  GenerateRecurrenceUnitTensors(
      int Nmax,
      const std::set<int>& lgi_operator_subset,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      std::map<spncci::NnPair,std::set<int>>& operator_subsets
    )
  {

    u3::SU3 x0; 
    HalfInt S0;
    int etap,eta;

    // Initialize operator_subsets with indices of unit tensor subspace whose rme's were computed between lgi
    std::set<int>& operator_subset=operator_subsets[spncci::NnPair(0,0)];
    for(int unit_tensor_subspace_index : lgi_operator_subset)
        operator_subset.insert(unit_tensor_subspace_index);

    // Recursively generate unit tensor family labels which will appear in recurrence and get corresponding
    // unit tensor subspace indices for each <Nnp,Nn> pair 
    for(int Nsum=0; Nsum<=2*Nmax; Nsum+=2)
      for(int Nn=0; Nn<=std::min(Nsum,Nmax); Nn+=2)
        {
          int Nnp=Nsum-Nn;
          // Nnp must be non-negative and less than truncation
          if((Nnp<0)||(Nnp>Nmax))
            continue;

          // For each source unit tensor subspace (coming from Nnp,Nn sector), 
          // generate related unit tensors in Nnp,Nn+2 sector and Nnp+2,Nn+2 sector
          spncci::NnPair NnpNn2(Nnp,Nn+2);
          spncci::NnPair NnpNn2_conjugate(Nn+2,Nnp);

          const std::set<int>& source_operator_subspace_indices=operator_subsets[spncci::NnPair(Nnp,Nn)];
          for(int source_operator_subspace_index : source_operator_subspace_indices)
            {
              // Extract source unit tensor labels
              std::tie(x0,S0,etap,eta)=unit_tensor_space.GetSubspace(source_operator_subspace_index).labels();
              
              // if etap-2 is negative, then not a valid unit tensor subspace
              if((etap-2)>=0)
                {
                  // Get list of possible x0p values based on requirement that x0p in (etap-2,0)x(0,eta)
                  MultiplicityTagged<u3::SU3>::vector x0p_set1=KroneckerProduct(u3::SU3(etap-2,0), u3::SU3(0,eta));
                  // For each possible x0p
                  for(auto& x0p_tagged : x0p_set1)
                    {
                      u3::SU3 x0p(x0p_tagged.irrep);

                      // look up unit tensor index 
                      u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0p,S0,etap-2,eta);
                      int operator_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
                      
                      // If unit tensor subspace exists, add to operator_subsets for given Nnp,Nn sector
                      if(operator_subspace_index!=-1)
                        operator_subsets[NnpNn2].insert(operator_subspace_index);

                      // Look up conjugate unit tensor subspace
                      u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0p),S0,eta,etap-2);
                      int operator_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
                      
                      // If conjugate unit tensor subspace exists, add to NnpNn2 conjugate list
                      if(operator_subspace_index_conj!=-1)
                        operator_subsets[NnpNn2_conjugate].insert(operator_subspace_index_conj);
                    }
                } 

              // Get list of possible x0p values based on x0p in (etap,0)x(0,eta+2)
              MultiplicityTagged<u3::SU3>::vector x0p_set2=KroneckerProduct(u3::SU3(etap,0), u3::SU3(0,eta+2)); 
              for(auto& x0p_tagged : x0p_set2)
                {
                  u3::SU3 x0p(x0p_tagged.irrep);
                  
                  // look up unit tensor subspace index 
                  u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0p,S0,etap,eta+2);
                  int operator_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
                  
                  // check if unit tensor subspace exits 
                  if(operator_subspace_index!=-1)
                    operator_subsets[NnpNn2].insert(operator_subspace_index);

                  // look up conjugate unit tensor subspace 
                  u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0p),S0,eta+2,etap);
                  int operator_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
                  
                  // check if conjugate unit tensor subspace exits 
                  if(operator_subspace_index_conj!=-1)
                    operator_subsets[NnpNn2_conjugate].insert(operator_subspace_index_conj);
                }

              // Add all source unit tensor subspace to Nnp+2,Nn+2 target sector.  
              // All source unit tensor subspace will appear in Nnp+2,Nn+2 sector 
              operator_subsets[spncci::NnPair(Nnp+2,Nn+2)].insert(source_operator_subspace_index);

              // Add all conjugates of source unit tensors to Nn+2,Nnp+2 sector
              u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
              int operator_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
              
              // Check the conjugate unit tensor subspace exists
              if(operator_subspace_index_conj!=-1)
                operator_subsets[spncci::NnPair(Nn+2,Nnp+2)].insert(operator_subspace_index_conj);
            }
        }
  }


void GenerateRecurrenceHypersectors(
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const std::map<std::pair<int,int>,std::set<int>>& lgi_unit_tensor_subset,
    int Nmax, int irrep_family_index_bra, int irrep_family_index_ket,
    std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
    spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors
  )
{
  //////////////////////////////////////////////////////////////////////////////////////////
  // get seeds and their labels for lgi pair
  std::pair<int,int> lgi_pair(irrep_family_index_bra,irrep_family_index_ket);    

  // Generate set of recurrence unit tensors for given starting seed subset 
  //
  // get starting set
  const std::set<int>& lgi_unit_tensors=lgi_unit_tensor_subset.at(lgi_pair);
  // Gereate set
  std::map<spncci::NnPair,std::set<int>> unit_tensor_subsets;
  spncci::GenerateRecurrenceUnitTensors(
    Nmax,lgi_unit_tensors,
    unit_tensor_space,unit_tensor_subsets
    );

  // Repeat for conjugate lgi pair 
  //
  // Define conjugate pair 
  std::pair<int,int> lgi_pair_conjugate(irrep_family_index_ket,irrep_family_index_bra);
  
  // get set of unit tensors for conjugate pair
  const std::set<int>& lgi_unit_tensors_conjugate=lgi_unit_tensor_subset.at(lgi_pair_conjugate);
  // Generate recurrence for conjugate, accumulating in unit_tensor_subset.
  spncci::GenerateRecurrenceUnitTensors(
    Nmax,lgi_unit_tensors_conjugate,
    unit_tensor_space,unit_tensor_subsets
    );

  // generate baby spncci hypersectors for given irrep family from unit tensor subspace
  // Resize 
  unit_tensor_hypersector_subsets.resize(2*Nmax);
  baby_spncci_hypersectors=spncci::BabySpNCCIHypersectors(
      baby_spncci_space, unit_tensor_space, 
      unit_tensor_subsets, unit_tensor_hypersector_subsets,
      irrep_family_index_bra, irrep_family_index_ket
    );
}



void PopulateHypersectorsWithSeeds(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const std::vector<int>& unit_tensor_hypersector_subset,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    const std::map< std::pair<int,int>, std::map<std::pair<int,int>, basis::OperatorBlocks<double>>>& 
      lgi_unit_tensor_blocks,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  )
{
  // get seeds for given lgi pair
  std::pair<int,int> lgi_pair(irrep_family_index_bra,irrep_family_index_ket);
  auto& seed_blocks=lgi_unit_tensor_blocks.at(lgi_pair);

  // get set of unit tensors for conjugate pair
  // conjugate lgi
  std::pair<int,int> lgi_pair_conjugate(irrep_family_index_ket,irrep_family_index_bra);
  const auto& seed_blocks_conjugate=lgi_unit_tensor_blocks.at(lgi_pair_conjugate);

  // Populate hypersectors with seeds
  for(int hypersector_index : unit_tensor_hypersector_subset)
    {
      const auto& hypersector=baby_spncci_hypersectors.GetHypersector(hypersector_index);
      int unit_tensor_subspace_index=hypersector.operator_subspace_index();
      int rho0=hypersector.multiplicity_index();
      std::pair<int,int> seed_unit_tensor_key(unit_tensor_subspace_index,rho0);
      
      // Check if hypersector is conjugate
      const auto& bra_subspace=baby_spncci_space.GetSubspace(hypersector.bra_subspace_index());
      const auto& ket_subspace=baby_spncci_space.GetSubspace(hypersector.ket_subspace_index());

      // If conjugate then get blocks from seed_blocks_conjugate
      // otherwise, get seeds from seed blocks.
      bool is_conjugate=(bra_subspace.irrep_family_index()>ket_subspace.irrep_family_index());
      const basis::OperatorBlocks<double>& seeds=is_conjugate?
      seed_blocks_conjugate.at(std::pair<int,int>(unit_tensor_subspace_index,rho0)):
      seed_blocks.at(std::pair<int,int>(unit_tensor_subspace_index,rho0));
      
      for(int i=0; i<seeds.size(); ++i)
        {
          if(seeds[i].rows()==0)
            continue;
          
          unit_tensor_hyperblocks[hypersector_index][i]=seeds[i];
        }
    }
}


  void CheckUnitTensorRecurrence(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
    const std::string& relative_unit_tensor_filename_template,
    const u3shell::SpaceU3SPN& lsu3shell_space, 
    const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
    const spncci::SpNCCISpace& spncci_space,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const basis::MatrixVector& spncci_expansions,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  )
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Checking unit tensors
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Zero initialize hypersectors 
    basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_explicit;
    basis::SetHyperoperatorToZero(baby_spncci_hypersectors,unit_tensor_hyperblocks_explicit);

    const u3::U3& sigmap=spncci_space[irrep_family_index_bra].sigma();
    const u3::U3& sigma =spncci_space[irrep_family_index_ket].sigma();

    // for each unit tensor
    for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
      {
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor = lgi_unit_tensor_labels[unit_tensor_index];

        // get unit tensor labels 
        u3::SU3 x0; 
        HalfInt S0,T0,Sp,Tp,S,T;
        int etap,eta;
        std::tie(x0,S0,T0,etap,Sp,Tp,eta,S,T)=unit_tensor.FlatKey();

        // Look up unit tensor subspace
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int unit_tensor_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        auto& subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        // Get unit tensor index in subspace 
        int unit_tensor_state_index
          =subspace.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0), int(Sp),int(Tp),int(S),int(T)));

        // Construct lsu3shell sectors for unit tensor
        const bool spin_scalar = false;
        u3shell::SectorsU3SPN unit_tensor_sectors;
        unit_tensor_sectors = u3shell::SectorsU3SPN(lsu3shell_space,unit_tensor,spin_scalar);
        
        // read in lsu3shell rms for unit tensor 
        basis::MatrixVector unit_tensor_lsu3shell_blocks;
        std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);
        lsu3shell::ReadLSU3ShellRMEs(
            filename,
            lsu3shell_basis_table,lsu3shell_space,
            unit_tensor,unit_tensor_sectors,unit_tensor_lsu3shell_blocks
          );

        // Compute unit tensor hyperblocks from lsu3shell rmes using explicit basis construction
        spncci::ComputeUnitTensorSectorsExplicit(
          sigmap, sigma, unit_tensor,unit_tensor_space,
          lsu3shell_space,unit_tensor_sectors,unit_tensor_lsu3shell_blocks,
          baby_spncci_space, spncci_expansions,baby_spncci_hypersectors,
          unit_tensor_hyperblocks_explicit
        );

        // Compute conjugate unit tensor hyperblocks if not diagonal sectors 
        if(not (sigmap==sigma))
        {
          spncci::ComputeUnitTensorSectorsExplicit(
            sigma, sigmap, unit_tensor,unit_tensor_space,
            lsu3shell_space,unit_tensor_sectors,unit_tensor_lsu3shell_blocks,
            baby_spncci_space, spncci_expansions,baby_spncci_hypersectors,
            unit_tensor_hyperblocks_explicit
          );
        }
      } // end unit tensor index 

    // Comparing recurrence to explicit hyperblocks and printing error message if difference between
    // hyperblocks exceeds tolerance of 1e-4   
    bool errors=false;
    for(int i=0; i<unit_tensor_hyperblocks.size(); ++i)
      for(int j=0; j<unit_tensor_hyperblocks[i].size(); ++j)
        {
          auto& hypersector=baby_spncci_hypersectors.GetHypersector(i);
          int bra, ket, tensor, rho0;
          std::tie(bra,ket,tensor,rho0)=hypersector.Key();
          auto& bra_subspace=baby_spncci_space.GetSubspace(bra);
          auto& ket_subspace=baby_spncci_space.GetSubspace(ket);
          auto& tensor_subspace=unit_tensor_space.GetSubspace(tensor);
          const Eigen::MatrixXd matrix1=unit_tensor_hyperblocks[i][j];
          const Eigen::MatrixXd matrix2=unit_tensor_hyperblocks_explicit[i][j];

          if(not mcutils::IsZero(matrix1-matrix2, 1e-4))
            {
              errors=true;
              std::cout<<"hyperblock "<<i<<" sub-block "<<j<<" is not correct"<<std::endl;
              std::cout<<bra_subspace.LabelStr()<<"  "<<ket_subspace.LabelStr()<<"  "<<tensor_subspace.LabelStr()<<"  "
              << rho0<<std::endl;
              std::cout<<"the matrix should be "<<bra_subspace.size()<<" x "<<ket_subspace.size()<<std::endl;
              std::cout<<"gammma_max: "<<ket_subspace.gamma_max()<<" upsilon_max "<<ket_subspace.upsilon_max()<<std::endl;
              std::cout<<"matrix1"<<std::endl<<matrix1<<std::endl<<"matrix2"
              <<std::endl<<matrix2<<std::endl;
            }
        }
    // If no error found, print no errors.
    assert(not errors);
    if(not errors)
      std::cout<<"no errors"<<std::endl;

  }




  void ConstructBranchedObservables(
    u3::WCoefCache& w_cache,
    const spncci::SpaceU3S& space_u3s,
    const std::vector<std::vector<spncci::SectorLabelsU3S>>& observable_sectors_u3s,
    const std::vector<basis::MatrixVector>& observable_matrices_u3s,
    std::map<HalfInt,spncci::SpaceLS>& spaces_lsj,
    int num_observables,
    const std::vector<HalfInt>& J_values,
    const std::vector<int>& observable_Jvalues,
    std::vector<std::map<spncci::JPair,spncci::MatrixType>>& observable_matrices
    )
  {
    // populate fully-branched many-body matrices for observables
    // map: observable -> J ->  matrix
    // std::vector<std::map<HalfInt,Eigen::MatrixXd>> observable_matrices;  
    observable_matrices.resize(num_observables);
    for (int observable_index=0; observable_index<num_observables; ++observable_index)
      {
        int J0=observable_Jvalues[observable_index];
        for (const HalfInt bra_J : J_values)
          for (const HalfInt ket_J : J_values)  
            {
              if(not am::AllowedTriangle(bra_J,J0,ket_J))
                continue;
              spncci::JPair JpJ(bra_J,ket_J);
              // set up aliases (for current observable and J space)
              const std::vector<spncci::SectorLabelsU3S>& sectors_u3s = observable_sectors_u3s[observable_index];
              const basis::MatrixVector& matrices_u3s = observable_matrices_u3s[observable_index];
              
              // determine allowed LS sectors
              const spncci::SpaceLS& bra_space_lsj = spaces_lsj[bra_J];
              const spncci::SpaceLS& ket_space_lsj = spaces_lsj[ket_J];

              Eigen::MatrixXd& observable_matrix = observable_matrices[observable_index][JpJ];

              // determine set of (L0,S0) labels for this observable (triangular with J0)
              std::vector<spncci::OperatorLabelsLS> operator_labels_ls;
              // Note: to update when J0 varies by observable
              spncci::GenerateOperatorLabelsLS(J0,operator_labels_ls);

              std::vector<spncci::SectorLabelsLS> sectors_lsj;
              spncci::GetSectorsLS(bra_space_lsj,ket_space_lsj,operator_labels_ls,sectors_lsj);

              // branch LS sectors to LSJ
              basis::MatrixVector matrices_lsj;  
              spncci::ContractAndRegroupLSJ(
                  bra_J,J0,ket_J,
                  w_cache,
                  space_u3s,sectors_u3s,matrices_u3s,
                  bra_space_lsj,ket_space_lsj,sectors_lsj,matrices_lsj
                );

              // collect LSJ sectors into J matrix
              //
              // Note: Interface needs to be generalized to handle J_bra != J_ket.
              ConstructOperatorMatrix(
                  bra_space_lsj,ket_space_lsj,sectors_lsj,matrices_lsj,
                  observable_matrix
                );
            }
      }
  }


void 
  SolveHamiltonian(
      const spncci::MatrixType& hamiltonian_matrix,
      const HalfInt& J,
      int num_eigenvalues,
      int eigensolver_num_convergence,  // whatever exactly this is...
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      std::map<HalfInt,Eigen::VectorXd>& eigenvalues,  // map: J -> eigenvalues
      std::map<HalfInt,spncci::MatrixType>& eigenvectors  // map: J -> eigenvectors
    )
  {    

    // set up aliases
    // spncci::MatrixType& hamiltonian_matrix = observable_matrices[0][J];
    std::cout << fmt::format("  Diagonalizing: J={}",J) << std::endl;

    // define eigensolver and compute
    typedef spncci::MatrixType MatrixType;  // allow for possible future switch to more compact single-precision matrix
    typedef double FloatType;
    Spectra::DenseSymMatProd<FloatType> matvec(hamiltonian_matrix);
    Spectra::SymEigsSolver<FloatType,Spectra::SMALLEST_ALGE,Spectra::DenseSymMatProd<FloatType>>
      eigensolver(
          &matvec,
          num_eigenvalues,
          eigensolver_num_convergence
        );
    eigensolver.init();
    int converged_eigenvectors = eigensolver.compute(
        eigensolver_max_iterations,
        eigensolver_tolerance,
        Spectra::SMALLEST_ALGE  // sorting rule
      );
    int eigensolver_status = eigensolver.info();
    std::cout
      << fmt::format("  Eigensolver reports: eigenvectors {} status {}",converged_eigenvectors,eigensolver_status)
      << std::endl;
    assert(converged_eigenvectors=eigensolver.eigenvalues().size());  // should this always be true?
    assert(converged_eigenvectors=eigensolver.eigenvectors().cols());  // should this always be true?

    // save eigenresults
    eigenvalues[J] = eigensolver.eigenvalues();
    eigenvectors[J] = eigensolver.eigenvectors();
    std::cout << fmt::format("  Eigenvalues (J={}):",J) << std::endl
              << mcutils::FormatMatrix(eigenvalues[J],"8.5f","    ")
              << std::endl;

    // check eigenvector norms
    Eigen::VectorXd eigenvector_norms(eigenvectors[J].cols());
    for (int eigenvector_index=0; eigenvector_index<converged_eigenvectors; ++eigenvector_index)
      {
        eigenvector_norms(eigenvector_index) = eigenvectors[J].col(eigenvector_index).norm();
        const double norm_tolerance=1e-8;
        assert(fabs(eigenvector_norms(eigenvector_index)-1)<norm_tolerance);
      }
      if (true)
        {
          std::cout << fmt::format("  Norms (J={}):",J) << std::endl
                    << mcutils::FormatMatrix(eigenvector_norms,"8.5f","    ")
                    << std::endl;
        }

      // normalize eigenvectors -- redundant with Spectra eigensolver
      for (int eigenvector_index=0; eigenvector_index<converged_eigenvectors; ++eigenvector_index)
        eigenvectors[J].col(eigenvector_index).normalize();

      // diagnostics
      if (false)
        {
          std::cout << fmt::format("  Eigenvectors -- norm (J={}):",J) << std::endl
                    << mcutils::FormatMatrix(eigenvectors[J],"8.5f","    ")
                    << std::endl;
        }
  }


}  // namespace
