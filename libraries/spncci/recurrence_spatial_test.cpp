/****************************************************************
  recurrence_seeds_test.cpp

  Anna E. McCoy
  INT

  SPDX-License-Identifier: MIT

 9/24/21 (aem): Created.
****************************************************************/
#include "spncci/recurrence_seeds.h"

#include <cppitertools/itertools.hpp>
#include <string>
#include <numeric>
#include <iostream>
#include <fstream>

#include "am/halfint_fmt.h"
#include "fmt/format.h"
#include "lgi/lgi.h"
#include "lgi/lgi_unit_tensors.h"
// #include "mcutils/profiling.h"
#include "mcutils/eigen.h"
#include "spncci/recurrence_indexing.h"
#include "spncci/recurrence_spatial.h"
#include "spncci/recurrence.h"
#include "spncci/spncci_basis.h"
#include "u3shell/relative_operator.h"



///////////////////////////////////////////////////////////////////////////////////////
// Mapping from output of recurrence to unit tensor hypersectors 
//
// Within Recurrence matrix (for fixed sigma,sigma',parity_bar), rows are indexed by 
// -> omega, omega'
//    -> upsilon, upsilon'
//      -> x0
//        -> rho0
//          -> Nbar,Nbar'
// while columns are 
// -> S,S'
//    -> Sp,Sn,Sp',Sn'
//      -> gamma,gamma'
//        -> S0,T0,Sbar,Sbarp,Tbar,Tbarp
//
// Hyperblocks are for fixed lgi pair, i.e., fixed sigma' Sp' Sn' S' and sigma Sp Sn S
// Each hyperblock corresponds to omega',omega, rho0
// Within each block, rmes are index by (gamma',upsilon')x(gamma,upsilon)   
//
basis::OperatorHyperblocks<double> MapRecurrenceBlockToHypersectors(
  const spncci::spatial::RecurrenceSpace& spatial_recurrence_space,
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>& spin_recurrence_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& unit_tensor_hypersectors,
  // const spncci::spatial::RecurrenceU3Sectors& recurrence_u3_sectors,
  const basis::OperatorBlocks<double> seed_blocks
  // const spncci::recurrence::SpatialRecurrenceMatrix& spatial_recurrence_matrix
  )
{
  basis::OperatorHyperblocks<double> unit_tensor_hyperblocks;
  basis::SetHyperoperatorToZero(unit_tensor_hypersectors,unit_tensor_hyperblocks);
  for(int hypersector_index=0; hypersector_index<unit_tensor_hypersectors.size(); ++hypersector_index)
    {
      const auto& hypersector=unit_tensor_hypersectors.GetHypersector(hypersector_index);
      const auto& [baby_spncci_index_bra,baby_spncci_index_ket,unit_tensor_subspace_index,rho0] 
              = hypersector.Key();
      const auto& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
      const auto& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);
      const auto& [sigma_bra,Sp_bra,Sn_bra,S_bra,omega_bra] = baby_spncci_subspace_bra.labels();
      const auto& [sigma_ket,Sp_ket,Sn_ket,S_ket,omega_ket] = baby_spncci_subspace_ket.labels();
      const int upsilon_max_bra=baby_spncci_subspace_bra.upsilon_max();
      const int upsilon_max_ket=baby_spncci_subspace_ket.upsilon_max();
      const int gamma_max_bra=baby_spncci_subspace_bra.gamma_max();
      const int gamma_max_ket=baby_spncci_subspace_ket.gamma_max();

      const auto& unit_tensor_subspace = unit_tensor_space.GetSubspace(unit_tensor_subspace_index);
      const auto& [x0,S0,Nbarp,Nbar] = unit_tensor_subspace.labels();
      // std::cout<<unit_tensor_subspace.LabelStr()<<std::endl;
      basis::OperatorBlocks<double>& hyperbocks = unit_tensor_hyperblocks[hypersector_index];

      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      // spatial::RecurrenceSpace() []
      // ->spatial::RecurrenceSp3RSpace() [sigma,sigma',parity_bar]
      //   ->spatial::RecurrenceNnsumSpace() [Nsum]
      //     ->spatial::RecurrenceU3Space() [omega,omega']->(upsilon x upsilon')
      //       ->spatial::RecurrenceOperatorSubspace() [x0] -> rho0
      //         ->spatial::RecurrenceOperatorState() [Nbar,Nbar']
      /////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Even though sigma,sigma' is fixed, we need to get the appropriate RecurrenceSp3RSpace for given parity par
      int sp3r_recurrence_index = spatial_recurrence_space.LookUpSubspaceIndex({sigma_ket,sigma_bra,Nbar%2});
      const auto& sp3r_recurrence_space = spatial_recurrence_space.GetSubspace(sp3r_recurrence_index);
      const auto& seed_block = seed_blocks[sp3r_recurrence_index];

      // std::cout<<"seed block "<<std::endl;
      // std::cout<<seed_block<<std::endl<<std::endl;

      // // Nnsum
      unsigned int Nnsum(omega_ket.N() - sigma_ket.N() + omega_bra.N() - sigma_bra.N());
      const auto& Nnsum_space = sp3r_recurrence_space.LookUpSubspace(Nnsum);

      // omega, omega'
      int u3_recurrence_space_index = Nnsum_space.LookUpSubspaceIndex({omega_ket,omega_bra});
      const auto& u3_recurrence_space = Nnsum_space.GetSubspace(u3_recurrence_space_index); 
      // std::cout<<"hi"<<std::endl;
      // TODO: reinstate. Removing temporarily since only checking seeds 
      // [omega,omega']x[sigma,sigma',parity_bar] block
      // const basis::OperatorBlock<double>& recurrence_omega_block
      //   = spatial_recurrence_matrix.GetRecurrenceBlock(Nnsum)[u3_recurrence_space_index];

      // Calculate unit tensor rmes from recurrence matrix and seeds
      // Product matrix is spncci::spatial::RecurrenceU3Space x spncci::spin::RecurrenceLGISpace
      
      // TODO: reinstate calculation once recurrence implemented 
      // basis::OperatorBlock<double> unit_tensor_omega_block = recurrence_omega_block*seed_block;
      basis::OperatorBlock<double> unit_tensor_omega_block = seed_block;

      ////////////////////////////////////////////////////////////////////////////////      
      // Reorganize unit_tensor_omega_block into hyperblocks
      ////////////////////////////////////////////////////////////////////////////////
      // Get index and offset information for spatial operators
      int spatial_operator_index = u3_recurrence_space.LookUpSubspaceIndex(x0);
      // std::cout<<spatial_operator_index<<std::endl;
      int spatial_operator_subspace_offset = u3_recurrence_space.GetSubspaceOffset(spatial_operator_index,rho0);
      const auto& spatial_operator_subspace = u3_recurrence_space.GetSubspace(spatial_operator_index);
      int spatial_operator_state_index = spatial_operator_subspace.LookUpStateIndex({Nbar,Nbarp});
      // std::cout<<"got everything"<<std::endl;
      // Get index and offset information for spins and spin operators 
      ////////////////////////////////////////////////////////////////////////////////
      // spin::RecurrenceSpace() []
      //   spin::RecurrenceLGISpace() [sigma,sigma',exchange_symm_bar]
      //   ->spin::RecurrenceSpinSpace() [S,S']
      //     -> spin::RecurrenceSpinSubspace() [Sp,Sn,Sp',Sn']/[T,T']->(gamma,gamma')
      //       ->spin::RecurrenceOperatorState() [S0,T0,Sbar,Sbar',Tbar,Tbar']
      ////////////////////////////////////////////////////////////////////////////////
      const auto& recurrence_lgi_space = spin_recurrence_space.LookUpSubspace({sigma_ket,sigma_bra, (Nbar+1)%2});
      // S, S'
      int spin_space_index = recurrence_lgi_space.LookUpSubspaceIndex({S_ket,S_bra});
      int spin_space_offset=recurrence_lgi_space.GetSubspaceOffset(spin_space_index);
      const auto& spin_space = recurrence_lgi_space.GetSubspace(spin_space_index);
      // std::cout<<"here"<<std::endl;
      // Sp, Sn, Sp', Sn
      int spin_subspace_index = spin_space.LookUpSubspaceIndex({{Sp_ket,Sn_ket},{Sp_bra,Sn_bra}});
      // std::cout<<"spin subspace index "<<spin_subspace_index<<std::endl;
      if (spin_subspace_index ==-1)
        continue;

      const auto& spin_subspace = spin_space.GetSubspace(spin_subspace_index);
      // std::cout<<"got all the spins"<<std::endl;
      for(int gamma_ket=1; gamma_ket<=gamma_max_ket; ++gamma_ket)
        for(int gamma_bra=1; gamma_bra<=gamma_max_bra; ++gamma_bra)
          {
            // std::cout<<"gammas "<<gamma_ket<<" "<<gamma_bra<<std::endl;
            int spin_subspace_offset = spin_space.GetSubspaceOffset(spin_subspace_index,gamma_ket,gamma_bra);

            ///////////////////////////////////
            for(int upsilon_ket=1; upsilon_ket<=upsilon_max_ket; ++upsilon_ket)
              for(int upsilon_bra=1; upsilon_bra<=upsilon_max_bra; ++upsilon_bra)
                {
                  // std::cout<<"upsilons "<<upsilon_ket<<" "<<upsilon_bra<<std::endl;
                  int u3_recurrence_space_offset = Nnsum_space.GetSubspaceOffset(u3_recurrence_space_index,upsilon_ket,upsilon_bra);
                  
                  // Want tile for fixed omega',omega, upsilon',upsilon, x0, rho0 and fixed S,S', SpSn, Sp'Sn', gamma,gamma'
                  const auto& operator_tile = unit_tensor_omega_block.block(
                      spatial_operator_subspace_offset, spin_space_offset+spin_subspace_offset,
                      spatial_operator_subspace.dimension(), spin_subspace.dimension()
                      );
          
                  // std::cout<<fmt::format(" operator tile {}  {}  {}  {}",spatial_operator_subspace_offset, spin_space_offset+spin_subspace_offset,
                  //     spatial_operator_subspace.dimension(), spin_subspace.dimension())<<std::endl;
                  // // std::cout<<operator_tile<<std::endl<<std::endl;;
                  for(int unit_tensor_state_index=0; unit_tensor_state_index<unit_tensor_subspace.size(); ++unit_tensor_state_index)
                    {
                      const auto&[T0,Sbarp,Tbarp,Sbar,Tbar]=unit_tensor_subspace.GetState(unit_tensor_state_index).labels();
                      // std::cout<<fmt::format("U[{}{}] {} {} {} {}",S0,T0, Sbarp,Tbarp,Sbar,Tbar)<<std::endl;
                      auto& hyperblock=hyperbocks[unit_tensor_state_index];

                      //S0,T0,Sbar, Sbar', Tbar, Tbar'
                      int spin_operator_state_index = spin_subspace.LookUpStateIndex(
                            {int(S0),int(T0),int(Sbar),int(Sbarp),int(Tbar),int(Tbarp)}
                          );
                      if (spin_operator_state_index == -1)
                        continue;

                      int row = (gamma_bra-1)*upsilon_max_bra+(upsilon_bra-1);
                      int col = (gamma_ket-1)*upsilon_max_ket+(upsilon_ket-1);


                      // std::cout<<"row and col "<<row<<" "<<col<<std::endl;
                      // std::cout<<spatial_operator_state_index<<" x "<<spin_operator_state_index<<std::endl;
                      hyperblock(row,col)=operator_tile(spatial_operator_state_index,spin_operator_state_index);
                      // std::cout<<hyperblock(row,col)<<std::endl<<std::endl;
                    } 
              }
        }
    }
  return unit_tensor_hyperblocks;
}



int main(int argc, char** argv)
{
  if (argc<4)
  {
    std::cout<<"Usage: recurrence seeds_test <Z> <N> <Nmax>"<<std::endl;
  }

  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nmax=std::stoi(argv[3]);
  // int N1v=std::stoi(argv[4]);

  const lgi::NuclideType nuclide({Z,N});
  int N1v = spncci::ValenceShellForNuclide(nuclide);

  const std::string lgi_filename="lgi_families.dat";
  
  HalfInt Nsigma0 = lgi::Nsigma0ForNuclide(nuclide,true);
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  lgi::ReadLGISet(lgi_filename,Nsigma0,lgi_vector);

  std::vector<int> lgi_full_space_index_lookup;
  lgi::ReadLGILookUpTable(lgi_full_space_index_lookup, lgi_vector.size());
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Set up recurrence spin and spatial indexing 
  // Read seeds from file.  Seeds for each sigma,sigma',parity_bar are stored in a single block
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"setting up recurrence indexing "<<std::endl;
  spncci::spin::Space<lgi::LGI> spin_space(lgi_vector, Nmax);
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST> spin_recurrence_space(spin_space, spin_space);

  auto it =iter::imap([](MultiplicityTagged<lgi::LGI> l) { return l.irrep.U3(); }, lgi_vector)| iter::unique_everseen;
  const spncci::spatial::Space spatial_space(std::vector<u3::U3>(it.begin(), it.end()), Nsigma0, Nmax);
  const spncci::spatial::RecurrenceSpace spatial_recurrence_space(spatial_space,spatial_space,N1v,Nsigma0);
  
  std::cout<<"Read seeds from file "<<std::endl;
  basis::OperatorBlocks<double> 
  seed_blocks=spncci::recurrence::GetRecurrenceSeedsFromFile(
        lgi_vector,lgi_full_space_index_lookup,
        spatial_recurrence_space,spin_recurrence_space
      );

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Enumerate unit tensor space
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"set up unit tensor space"<<std::endl;
  int J0_for_unit_tensors = -1;  // all J0
  int T0_for_unit_tensors = -1;  // all T0
  const bool restrict_positive_N0 = false;  // don't restrict to N0 positive
  // get full set of possible unit tensor labels up to Nmax, N1v truncation
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(
    Nmax, N1v,unit_tensor_labels,J0_for_unit_tensors,
    T0_for_unit_tensors,restrict_positive_N0
  );
  // generate unit tensor space
  u3shell::RelativeUnitTensorSpaceU3S unit_tensor_space(Nmax,N1v,unit_tensor_labels);
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Set up baby spncci hypersectors
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<<"set up baby spncci space"<<std::endl;
  spncci::NmaxTruncator truncator(Nsigma0,Nmax);  
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;
  bool restrict_sp3r_to_u3_branching=false;
  spncci::GenerateSpNCCISpace(lgi_vector,truncator,spncci_space,sigma_irrep_map,restrict_sp3r_to_u3_branching);
  spncci::BabySpNCCISpace baby_spncci_space=spncci::BabySpNCCISpace(spncci_space);
  
  std::vector<spncci::LGIPair> lgi_pair_vector;
  for(int irrep_family_index_ket=0; irrep_family_index_ket<lgi_vector.size(); irrep_family_index_ket++)
    for(int irrep_family_index_bra=0; irrep_family_index_bra<lgi_vector.size(); irrep_family_index_bra++)
      {
        lgi_pair_vector.emplace_back(irrep_family_index_bra,irrep_family_index_ket);
      }

  spncci::OperatorBlocks lgi_transformations; // Initalized but never populated or used
  bool transform_lgi_families = false;
  //Populated in DoRecurrenceInitialization
  std::map<spncci::NnPair,std::set<int>> unit_tensor_subspace_subsets; 
  spncci::BabySpNCCIHypersectors baby_spncci_hypersector_seeds;
  spncci::BabySpNCCIHypersectors baby_spncci_hypersector_seeds_conj;
  basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_seeds;
  basis::OperatorHyperblocks<double> unit_tensor_hyperblocks_seeds_conj;

  std::cout<<"for each lgi pair -- compare"<<std::endl;
  for(const spncci::LGIPair& lgi_pair : lgi_pair_vector)
    { 
      const auto& [irrep_family_index_bra,irrep_family_index_ket] = lgi_pair;
      // if ((irrep_family_index_bra!=4) || (irrep_family_index_ket!=5))
      //   continue;
      // std::cout<<fmt::format("{} {}",
      //   lgi_vector[irrep_family_index_bra].irrep.Str(),
      //   lgi_vector[irrep_family_index_ket].irrep.Str())<<std::endl;
      // std::cout<<irrep_family_index_bra<<" "<<irrep_family_index_ket<<std::endl;

      // std::cout<<"do recurrence initalization "<<std::endl;
      spncci::DoRecurrenceInitialization(Nmax,N1v,lgi_pair,lgi_vector,lgi_full_space_index_lookup,
        baby_spncci_space,unit_tensor_space,lgi_transformations,transform_lgi_families,
        unit_tensor_subspace_subsets,
        baby_spncci_hypersector_seeds,baby_spncci_hypersector_seeds_conj,
        unit_tensor_hyperblocks_seeds,unit_tensor_hyperblocks_seeds_conj
      );

      // std::cout<<"map to hyperblocks"<<std::endl;
      basis::OperatorHyperblocks<double> 
      unit_tensor_hyperblocks_comparison
        =MapRecurrenceBlockToHypersectors(
          spatial_recurrence_space,spin_recurrence_space,
          baby_spncci_space,unit_tensor_space,
          baby_spncci_hypersector_seeds,seed_blocks
        );

      // std::cout<<"do comparison"<<std::endl;
      for(int i=0; i<unit_tensor_hyperblocks_comparison.size(); ++i)
        {
          const auto& blocks1 = unit_tensor_hyperblocks_seeds[i];
          const auto& blocks2 = unit_tensor_hyperblocks_comparison[i];
          for(int j=0; j<blocks1.size(); ++j)
            {
              const auto& block1 = blocks1[j];
              const auto& block2 = blocks2[j];

              if(not mcutils::IsZero(block1-block2))
                std::cout<<block1<<std::endl<<std::endl<<block2
              <<std::endl<<"------------"<<std::endl;
            }
        }


    }
}
