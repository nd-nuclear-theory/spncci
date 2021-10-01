/****************************************************************
  recurrence_seeds.cpp

  Anna E. McCoy[1] and Patrick J. Fasano[2,3]
  [1] Institute for Nuclear Theory
  [2] University of Notre Dame
  [3] Lawrence Berkeley National Laboratory

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/recurrence_seeds.h"

#include <algorithm>
#include <fstream>
#include <iostream>

#include "cppitertools/itertools.hpp"
#include "spncci/io_control.h"
// #include "cppitertools/enumerate.hpp"

#include "basis/basis.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"


#include "spncci/spncci_basis.h"
#include "u3shell/unit_tensor_space_u3s.h"

namespace spncci
{
namespace recurrence
{
using LGILookupTable = std::unordered_map<u3::U3,std::vector<std::pair<int,lgi::LGI>>,boost::hash<u3::U3>>;

LGILookupTable LGIGroups(const lgi::MultiplicityTaggedLGIVector& lgi_vector)
// Regroup mutiplicity tagged lgi by sigma
// Retruns a map keyed by sigma with value (irrep_family_index,lgi)
{
  std::unordered_map<u3::U3,std::vector<std::pair<int,lgi::LGI>>,boost::hash<u3::U3>> lgi_groups;
  
  auto extract_sigma=[](const auto& enum_lgi) {return enum_lgi.element.irrep.U3();};
  auto compare_sigma=[&extract_sigma](const auto& enum_lgi1,const auto& enum_lgi2) {return extract_sigma(enum_lgi1)<extract_sigma(enum_lgi2);};
  auto regrouped_lgi_vector=iter::groupby(iter::sorted(iter::enumerate(lgi_vector),compare_sigma),extract_sigma);
  
  for(auto&& [sigma,group] : regrouped_lgi_vector)
    {
      for (const auto& [index,tagged_lgi] : group)
      {
        lgi_groups[sigma].emplace_back(index,tagged_lgi.irrep);
      }
    }
  return lgi_groups;
}


///////////////////////////////////////////////////////////////////////////////////////
// Returns seed matrix for pair of [sigma_ket,sigma_bra] and parity_bar
// Matrix rows are 
//    -> x0, rho0, [Nbar, Nbarp]
//
// Matrix cols are
//    -> [S_ket,S_bra], [Sp_ket,Sn_ket,]
//
//
///////////////////////////////////////////////////////////////////////////////////////
basis::OperatorBlock<double>
ReadSeedBlockFromFile(
  const LGILookupTable& lgi_groups,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::spatial::RecurrenceU3Space&  lgi_spatial_space,
  const spncci::spin::RecurrenceLGISpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>& lgi_spin_space,
  basis::OperatorBlock<double>& seed_block
)
{
  // Sort lgi's by U(3)
  const auto&[sigma_ket,sigma_bra] = lgi_spatial_space.labels();
  int exchange_symm_bar = lgi_spin_space.exchange_symm_bar();
  int parity_bar = (exchange_symm_bar+1)%2;
  seed_block
    =basis::OperatorBlock<double>::Zero(lgi_spatial_space.dimension(),lgi_spin_space.dimension());

  for(const auto&[irrep_family_index_bra,lgi_bra] : lgi_groups.at(sigma_bra))
    for(const auto&[irrep_family_index_ket,lgi_ket] : lgi_groups.at(sigma_ket))
      {        
        const HalfInt& S_bra = lgi_bra.S(); 
        const HalfInt& S_ket = lgi_ket.S();
        const auto& [Sp_bra,Sn_bra] = lgi_bra.upstream_labels();
        const auto& [Sp_ket,Sn_ket] = lgi_ket.upstream_labels();
        int spin_space_index = lgi_spin_space.LookUpSubspaceIndex({S_ket,S_bra});
        
        if(spin_space_index==-1)
          continue;
        
        const auto& spin_space=lgi_spin_space.GetSubspace(spin_space_index);
        int spin_subspace_index = spin_space.LookUpSubspaceIndex({{Sp_ket,Sn_ket},{Sp_bra,Sn_bra}});
        if(spin_subspace_index==-1)
          continue;

        const int degeneracy_ket=spin_space.GetKetSubspaceDegeneracy(spin_subspace_index);
        const int degeneracy_bra=spin_space.GetBraSubspaceDegeneracy(spin_subspace_index);
        
        const auto& spin_subspace = spin_space.GetSubspace(spin_subspace_index);

        ///////////////////////////////////////////////////////////////////////////
        // Read in list of unit tensors between lgi pair and conjugates from files
        // Corresponding rho0 values stored separately for later hypersector lookup
        // but could just be recalculated.
        ///////////////////////////////////////////////////////////////////////////

        // Get index corresponding to lgi in the full space.
        // Index may differ from lgi index in basis if space has been truncated
        int index1=lgi_full_space_index_lookup[irrep_family_index_bra];
        int index2=lgi_full_space_index_lookup[irrep_family_index_ket];

        // Read in operators 
        std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensors;
        std::vector<int> rho0_values;

        std::string lgi_unit_tensor_filename=fmt::format("seeds/operators_{:06d}_{:06d}.dat",index1,index2);
        lgi::ReadUnitTensorLabels(lgi_unit_tensor_filename,lgi_unit_tensors,rho0_values);

        // Reads in unit tensor seed blocks and stores them in a vector of blocks. Order
        // corresponds to order of (unit_tensor,rho0) pairs in corresponding operator file. 
        basis::OperatorBlocks<double> unit_tensor_seed_blocks;
        std::string seed_filename=fmt::format("seeds/seeds_{:06d}_{:06d}.rmes",index1,index2);
        lgi::ReadBlocks(seed_filename, lgi_unit_tensors.size(), unit_tensor_seed_blocks);

        for(int t=0; t<lgi_unit_tensors.size(); ++t)
          {
            //Extract unit tensor labels 
            const auto& [x0,S0,T0,Nbarp,Sbarp,Tbarp,Nbar,Sbar,Tbar] = lgi_unit_tensors[t].FlatKey();
            if (Nbar%2!=parity_bar)
              continue;

            //Select source block 
            const auto& block = unit_tensor_seed_blocks[t];
            int rho0 = rho0_values[t];

            // Look up spin operator.  Note that cast to int necessary because unit tensor spin labels are 
            // currently HalfInt type even though they are always integer values. 
            // TODO: Patrick, is there a better way to do this casting?                
            int spin_operator_index 
                = spin_subspace.LookUpStateIndex({int(S0),int(T0),int(Sbar),int(Sbarp),int(Tbar),int(Tbarp)});
            //Get spatial operator indices 
            int spatial_operator_index = lgi_spatial_space.LookUpSubspaceIndex(x0);
            const auto& spatial_operator_subspace
                = lgi_spatial_space.GetSubspace(spatial_operator_index);
            
            int spatial_state_index = spatial_operator_subspace.LookUpStateIndex({Nbar,Nbarp});
            
            //Calculate target row index 
            int target_index_row 
              = lgi_spatial_space.GetSubspaceOffset(spatial_operator_index,rho0)
                + spatial_state_index;

            int g=1;

            for(int gamma_ket=1; gamma_ket<=degeneracy_ket; ++gamma_ket)
              for(int gamma_bra=1; gamma_bra<=degeneracy_bra; ++gamma_bra)
            // for (int g_ket=0; g_ket<degeneracy_ket; ++g_ket)
            //   for (int g_bra=0; g_bra<degeneracy_bra; ++g_bra)
                {
                  //Calcuate target column index
                  // Lookup for Subspace offset require degeneracy index g which is 
                  // g=(gamma_ket-1)*gammap_max + (gammap-1)
                  // Note: degeneracy is 1 based

                  int target_index_col 
                    = lgi_spin_space.GetSubspaceOffset(spin_space_index)
                      +spin_space.GetSubspaceOffset(spin_subspace_index,gamma_ket,gamma_bra)
                      + spin_operator_index;

                  seed_block(target_index_row,target_index_col)
                    = block(gamma_bra-1,gamma_ket-1);

                  // Increment g
                  g++;
                }
          }
      }    
  return seed_block;
}// end function

///////////////////////////////////////////////////////////////////////////////////////
basis::OperatorBlocks<double>
GetRecurrenceSeedsFromFile(
  const lgi::MultiplicityTaggedLGIVector& lgi_vector,
  const std::vector<int>& lgi_full_space_index_lookup,
  const spncci::spatial::RecurrenceSpace& spatial_recurrence_space,
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>& spin_recurrence_space
)
{
  // Sort lgi's by U(3)
  const auto& lgi_groups=LGIGroups(lgi_vector);
  
  basis::OperatorBlocks<double> seed_blocks(spatial_recurrence_space.size());
  for(int sp3r_space_index=0; sp3r_space_index<spatial_recurrence_space.size(); ++sp3r_space_index)
    {      
      const auto& recurrence_sp3r_space = spatial_recurrence_space.GetSubspace(sp3r_space_index);
      const auto& lgi_labels = recurrence_sp3r_space.labels();
      const auto&[sigma_ket,sigma_bra,parity_bar] = lgi_labels;

      // We only care about the RecurrenceNnsumSpace for Nsum=0
      // Within The RecurrenceNnsumSpace, there is only one lgi pair
      // and thus only one RecurrenceU3Space 
      const auto& lgi_spatial_space = recurrence_sp3r_space.GetSubspace(0).GetSubspace(0);
            
      // To enforce antisymmetry
      uint8_t exchange_symmetry = (parity_bar+1)%2;

      // Get corresponding spin space 
      const auto& lgi_spin_space=spin_recurrence_space.LookUpSubspace({sigma_ket,sigma_bra,exchange_symmetry});
      
      ReadSeedBlockFromFile(
        lgi_groups,lgi_full_space_index_lookup,
        lgi_spatial_space,lgi_spin_space,
        seed_blocks[sp3r_space_index]
      );

     }

  return seed_blocks;
}// end function

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
// Indexing is 
//

// Indexing within hypersectors is:
//    
//
void MapRecurrenceBlockToHypersectors(
  const spncci::spatial::RecurrenceSpace& spatial_recurrence_space,
  const spncci::spin::RecurrenceSpace<lgi::LGI, spncci::spin::UnitTensorLabelsST>& spin_recurrence_space,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const BabySpNCCIHypersectors& unit_tensor_hypersectors,
  const basis::OperatorBlocks<double> seed_blocks
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
      basis::OperatorBlocks<double>& hyperbocks = unit_tensor_hyperblocks[hypersector_index];

      // Select sigma,sigma' space and corresponding seed block
      int sp3r_recurrence_index = spatial_recurrence_space.LookUpSubspaceIndex({sigma_ket,sigma_bra,Nbar%2});
      const auto& sp3r_recurrence_space = spatial_recurrence_space.GetSubspace(sp3r_recurrence_index);
      const auto& seed_block = seed_blocks[sp3r_recurrence_index];

      // Get indexing for Nnsum space
      int Nnsum(omega_ket.N() - sigma_ket.N() + omega_bra.N() - sigma_bra.N());
      int Nnsum_space_index=sp3r_recurrence_space.LookUpSubspaceIndex(Nnsum);
      int Nnsum_offset = sp3r_recurrence_space.GetSubspaceOffset(Nnsum_space_index);

      // Get relevant omega, omega' space 
      const auto& u3_recurrence_space 
            = sp3r_recurrence_space.GetSubspace(Nnsum_space_index).LookUpSubspace({omega_ket,omega_bra});

      // For each unit tensor state
      // for(int unit_tensor_state_index=0; unit_tensor_state_index<unit_tensor_subspace.size(); ++unit_tensor_state_index)
      //   {
      //     const auto&[T0,Sbarp,Tbarp,Sbar,Tbar]=unit_tensor_subspace.GetState(unit_tensor_state_index).labels();
      //     auto& hyperblock=hyperbocks[unit_tensor_state_index];

      //     int spin_operator_index = spin_subspace.LookUpStateIndex(
      //           {int(S0),int(T0),int(Sbar),int(Sbarp),int(Tbar),int(Tbarp)}
      //         );
            
      //     //Get spatial operator indices 
      //     int spatial_operator_index = u3_recurrence_space.LookUpSubspaceIndex(x0);
      //     const auto& spatial_operator_subspace
      //         = u3_recurrence_space.GetSubspace(spatial_operator_index);
          
      //     int spatial_state_index 
      //         = spatial_operator_subspace.LookUpStateIndex({Nbar,Nbarp});

      //   }
    }
}

}  // namespace recurrence
}  // namespace spncci
