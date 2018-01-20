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

void GetLGIUnitTensorSubspaceIndices(
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
  std::set<int>& lgi_operator_subset
  )
  {
    // for each unit tensor, extract labels and identify subspace index for both the tensor and its conjugate
    for(auto unit_tensor : lgi_unit_tensors)
      {
        // Extract unit tensor labels 
        u3::SU3 x0; 
        HalfInt S0;
        int etap,eta;
        std::tie(x0,S0,std::ignore,etap,std::ignore,std::ignore,eta,std::ignore,std::ignore)=unit_tensor.FlatKey();

        // look up unit tensor index 
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
        int operator_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
        lgi_operator_subset.insert(operator_subspace_index);
        
        // Add conjugate tensor for Nn=0 sectors 
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);        
        int operator_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
        lgi_operator_subset.insert(operator_subspace_index_conj);


      }
  }

void GetUnitTensorSubspaceIndices(
  const u3::SU3& x0, HalfInt S0, int etap, int eta,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  std::set<int>& operator_subset
  )
  // for a given unit tensor subspace, look up its subspace index for both the tensor and its conjugate
  // and add the unit tensor subspace to the accumulating set of operator subspaces appear in the recurrence
  {
    // look up unit tensor index 
    u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels(x0,S0,etap,eta);
    int operator_subspace_index=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels);
    
    // If unit tensor subspace exists, add to operator_subsets for given Nnp,Nn sector
    if(operator_subspace_index!=-1)
      operator_subset.insert(operator_subspace_index);

    // Look up conjugate unit tensor subspace
    u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
    int operator_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);
    
    // If conjugate unit tensor subspace exists, add to NnpNn2 conjugate list
    if(operator_subspace_index_conj!=-1)
      operator_subset.insert(operator_subspace_index_conj);
  }

void GenerateRecurrenceUnitTensors(
    int Nrel_max,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensors,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    std::vector<int>& operator_subset
  )
  // TODO: use this function to replace below
  // Generate vector of operator subspaces which will appear in the recurrence from
  // unit tensors with non-zero rmes between the given lgi pair
  //
  // Note: Nrel_max=Nmax+2N1v
  {
    // Unit tensors subspaces are identified recursively starting from those between the lgi
    //
    // initialize subspace container 
    std::vector<std::set<int>> operator_subspace_indices_by_Nsum(Nrel_max/2+1);

    // Get lgi unit tensor subspaces 

    auto& lgi_operator_subset=operator_subspace_indices_by_Nsum[0];
    GetLGIUnitTensorSubspaceIndices(unit_tensor_space,lgi_unit_tensors,lgi_operator_subset);

    // Iteratively obtain unit tensor subspaces from lgi unit tensors for recurrence 
    // In each interaction, new unit tensors must satisfy:
    //    x0 x (2,0) -> x0p
    // and either 
    //    (etap-2,0)x(eta,0) -> x0p
    // or
    //    (etap,0)x(eta+2,0) -> x0p
    //
    for(int m=1; m<operator_subspace_indices_by_Nsum.size(); ++m)
      for(int source_operator_subspace_index : operator_subspace_indices_by_Nsum[m-1])
        {
          // Extract source unit tensor labels
          // initialize labels 
          u3::SU3 x0; 
          HalfInt S0;
          int etap,eta;
          std::tie(x0,S0,etap,eta)=unit_tensor_space.GetSubspace(source_operator_subspace_index).labels();
          
          // case 1
          
          if((etap-2)>=0)
            {
              // Get list of possible x0p values from etap and eta
              MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(u3::SU3(etap-2,0), u3::SU3(0,eta));
              // For each possible x0p
              for(auto& x0p_tagged : x0p_set)
                {
                  u3::SU3 x0p(x0p_tagged.irrep);

                  // If x0 x (2,0) -> x0p constraint satisfied, add unit tensors list
                  if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)>0)
                      GetUnitTensorSubspaceIndices(
                        x0p,S0,etap-2,eta,unit_tensor_space,
                        operator_subspace_indices_by_Nsum[m]
                      );               
                }
            }

          // case 2
            
          else if ((eta+2)<=Nrel_max)
            {
             // Get list of possible x0p values from etap and eta
              MultiplicityTagged<u3::SU3>::vector x0p_set=KroneckerProduct(u3::SU3(etap,0), u3::SU3(0,eta+2));
              // For each possible x0p
              for(auto& x0p_tagged : x0p_set)
                {
                  u3::SU3 x0p(x0p_tagged.irrep);

                  // If x0 x (2,0) -> x0p constraint satisfied, add unit tensors list
                  if(u3::OuterMultiplicity(x0p,u3::SU3(2,0),x0)>0)
                      GetUnitTensorSubspaceIndices(
                        x0p,S0,etap,eta+2,unit_tensor_space,
                        operator_subspace_indices_by_Nsum[m]
                      );               
                }
            }
        }

    // Flattening operator subsets into one vector
    for(auto& subset : operator_subspace_indices_by_Nsum)
      for(auto index : subset)
        operator_subset.push_back(index);
  }




  void GetUnitTensorSeedBlocks(
      const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const std::string& relative_unit_tensor_filename_template,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
      const basis::MatrixVector& lgi_expansions,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      std::map<std::pair<int,int>,std::map<std::pair<int,int>,basis::OperatorBlocks<double>>>& lgi_unit_tensor_blocks
    )
  {
    //If parallellize, iterate first over unit_tensor_index in parallel and then write to new map in serial?, maybe
    // Alternatively, we zero initalize lgi_unit_tensor_blocks
    // std::unordered_map<u3shell::SectorsU3SPN,basis::MatrixVector, boost::hash<u3shell::SectorsU3SPN>> 
    //    unit_tensor_spncci_blocks;
    for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
      {
        // Get labels of corresponding unit tensor 
        const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];

        // File name containing lsu3shell rmes of unit tensor 
        std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);

        // generate unit tensor sectors in lsu3shell basis
        const bool spin_scalar = false; // All spin values allowed
        u3shell::SectorsU3SPN unit_tensor_sectors(lsu3shell_space,unit_tensor_labels,spin_scalar);
      
        // read in lsu3shell rmes of unit tensor 
        basis::MatrixVector unit_tensor_lsu3shell_matrices;
        lsu3shell::ReadLSU3ShellRMEs(
            filename,
            lsu3shell_basis_table,lsu3shell_space,
            unit_tensor_labels,unit_tensor_sectors,unit_tensor_lsu3shell_matrices
          );

        // transform seed rmes to SpNCCI basis (among LGIs)
        basis::MatrixVector unit_tensor_spncci_matrices; // temporary container
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
        int unit_tensor_state_index
              =subspace.LookUpStateIndex(std::tuple<int,int,int,int,int>(int(T0),int(Sp),int(Tp),int(S),int(T)));

        // transfer spncci unit tensor rmes from temporary container to lgi hyperblocks for use in spncci recurrence
        for(int sector_index=0; sector_index<unit_tensor_sectors.size(); ++sector_index)
          {
            // Check that the sector has non-zero rmes as defined by the zero_tolerance
            if(mcutils::IsZero(unit_tensor_spncci_matrices[sector_index],spncci::g_zero_tolerance))
              continue;
            
            // extract U3SPN sector information from unit tensor sectors definded in lsu3shell basis
            const typename u3shell::SectorsU3SPN::SectorType& sector = unit_tensor_sectors.GetSector(sector_index);
            const int bra_subspace_index = sector.bra_subspace_index();
            const int ket_subspace_index = sector.ket_subspace_index();
            const u3shell::SubspaceU3SPN bra_subspace=sector.bra_subspace();
            const u3shell::SubspaceU3SPN ket_subspace=sector.ket_subspace();
            const u3::U3& bra_sigma = bra_subspace.U3();
            const u3::U3& ket_sigma = ket_subspace.U3();
            const int rho0 = sector.multiplicity_index();

            // Look up baby spncci index from labels
            spncci::BabySpNCCISubspaceLabels 
              baby_spncci_bra(bra_sigma,bra_subspace.Sp(),bra_subspace.Sn(), bra_subspace.S(),bra_sigma);
            spncci::BabySpNCCISubspaceLabels 
              baby_spncci_ket(ket_sigma,ket_subspace.Sp(),ket_subspace.Sn(), ket_subspace.S(), ket_sigma);

            int baby_spncci_index_bra=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_bra);
            int baby_spncci_index_ket=baby_spncci_space.LookUpSubspaceIndex(baby_spncci_ket);

            // if baby_spncci_index_bra or baby_spncci_index_ket not found, then subspace belongs to an irrep which
            // is not part of the truncated spncci space, so sector is ignored
            if(baby_spncci_index_bra==-1 || baby_spncci_index_ket==-1)
              continue;

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
    // unit_tensor_hypersector_subsets.resize(2*Nmax);
    unit_tensor_hypersector_subsets.resize(Nmax+1);
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
      
        // Check if there exists non-zero seeds 
        bool seeds_nonzero=
          is_conjugate?
          seed_blocks_conjugate.count(std::pair<int,int>(unit_tensor_subspace_index,rho0)):
          seed_blocks.count(std::pair<int,int>(unit_tensor_subspace_index,rho0));

        // if no non-zero seeds, then continue
        if(not seeds_nonzero)
          continue;

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





}  // namespace
