/****************************************************************
  branching_u3s.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching_u3s.h"

#include <fstream>
#include <iostream>

#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "spncci/unit_tensor.h"

namespace spncci
{


  std::string SubspaceU3S::Str() const
  {
    return fmt::format("{}",labels_.Str());
  }


  SubspaceU3S::SubspaceU3S(const u3::U3S& omegaS,const BabySpNCCISpace& baby_spncci_space)
    {
      labels_=omegaS;
      int sector_index=0;
      int state_index=0;
      for(int baby_spncci_index=0; baby_spncci_index<baby_spncci_space.size(); ++baby_spncci_index)
        {
          const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_index);
          if(omegaS==baby_spncci_subspace.omegaS())
            {
              assert(omegaS==baby_spncci_subspace.omegaS());
              int state_dim=baby_spncci_subspace.size();
              
              PushStateLabels(StateLabelsType(baby_spncci_index));
              sector_index_lookup_[state_index]=sector_index;
              // state index
              ++state_index;
              // starting position of state in sector row or column. 
              sector_index+=state_dim;
            }
        }
      sector_size_=sector_index;

    }

  SpaceU3S::SpaceU3S(const BabySpNCCISpace& baby_spncci_space)
  {
    for(int baby_spncci_index=0; baby_spncci_index<baby_spncci_space.size(); ++baby_spncci_index)
      {
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_index);
        u3::U3S omegaS(baby_spncci_subspace.omega(),baby_spncci_subspace.S());
        if(ContainsSubspace(omegaS))
          continue;
        SubspaceU3S subspace(omegaS,baby_spncci_space);
        PushSubspace(subspace);
      }
  }

  std::string SectorLabelsU3S::Str() const
  {
    return fmt::format("( {} {}  {}{} {} : {} {}  {}", 
      bra_index(),ket_index(), N0(), x0().Str(),S0(),kappa0(),L0(),rho0());
  }

  void GetSectorsU3S(
    const spncci::SpaceU3S& space, 
    const std::vector<u3shell::IndexedOperatorLabelsU3S>& relative_tensor_labels,
    std::vector<spncci::SectorLabelsU3S>& sector_vector
    )
  // TODO: replace with a normal sectors object
  {
    int u3s_sector_vector_index=0;
    for(auto tensor_labels:relative_tensor_labels)
      for(int j=0; j<space.size(); ++j)
        for(int i=0; i<space.size(); ++i)
          {
            int kappa0,L0; 
            u3shell::OperatorLabelsU3S op_labels;
            std::tie(op_labels, kappa0,L0)=tensor_labels;
            assert(kappa0!=0);
            u3::U3S omegapSp(space.GetSubspace(i).GetSubspaceLabels());
            u3::U3S omegaS(space.GetSubspace(j).GetSubspaceLabels());
            int rho0_max=u3::OuterMultiplicity(omegaS.SU3(), op_labels.x0(),omegapSp.SU3());
            // Check if allowed U(1) coupling
            if(omegaS.U3().N()+op_labels.N0()!=omegapSp.U3().N())
              continue;
            // Check if allowed SU(2) coupling 
            if(not am::AllowedTriangle(omegaS.S(), op_labels.S0(), omegapSp.S()))
              continue;
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                spncci::SectorLabelsU3S sector(i,j,op_labels, kappa0,L0,rho0);
                // should not alread be a sector in map
                // assert(not u3s_sectors.count(sector));
                // if(not u3s_sectors.count(sector))
                //   {
                //     u3s_sectors[sector]=u3s_sector_vector_index;
                //     ++u3s_sector_vector_index;
                //   }
                sector_vector.push_back(sector);
              }
          }
    // sector_vector.resize(u3s_sectors.size());
    // for(auto it=u3s_sectors.begin(); it!=u3s_sectors.end(); ++it)
    //   sector_vector[it->second]=it->first;
  }

  void 
  ContractAndRegroupU3S(
      int Nmax, int N1b,
      const std::vector<spncci::SectorLabelsU3S>& sector_labels_vector,
      const u3shell::RelativeRMEsU3ST& interaction_rme_cache,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      spncci::SpaceU3S& target_space,
      spncci::UnitTensorMatricesByIrrepFamily& unit_tensor_sector_cache,
      basis::MatrixVector& matrix_vector
    )
  {
    // Initial sectors to zero matrices
    matrix_vector.resize(sector_labels_vector.size());
    for(int s=0; s<matrix_vector.size(); ++s)
    {
      const spncci::SectorLabelsU3S& sector=sector_labels_vector[s];
      const spncci::SubspaceU3S& ket_subspace=target_space.GetSubspace(sector.ket_index());
      const spncci::SubspaceU3S& bra_subspace=target_space.GetSubspace(sector.bra_index());
      int sector_dim_bra=bra_subspace.sector_dim();
      int sector_dim_ket=ket_subspace.sector_dim();
      matrix_vector[s]=Eigen::MatrixXd::Zero(sector_dim_bra,sector_dim_ket);
    }

    // iterate over interaction get unit tensor,kappa0,L0
    // iterate over U3S sectors to get target sectors
    // get corresponding unit tensor sector
    // Contract interaction with unit tensors rmes and accumulate in U3S sectors.   
    for(auto it=interaction_rme_cache.begin(); it!=interaction_rme_cache.end(); ++it)
      {
        // Extract labels 
        int kappa0,L0;
        u3shell::RelativeUnitTensorLabelsU3ST tensor_u3st;
        std::tie(tensor_u3st,kappa0,L0)=it->first;
        double interaction_rme=it->second;

        //Check that unit tensor has rme between states in Nmax truncated basis
        int rp=tensor_u3st.bra().eta();
        int r=tensor_u3st.ket().eta();
        if((r>Nmax+2*N1b)||(rp>Nmax+2*N1b))
          continue;

        // Iterate over U3 sectors to get target sectors
        #pragma omp parallel for schedule(runtime)
        for(int s=0; s<sector_labels_vector.size(); ++s)
        {
          const spncci::SectorLabelsU3S& sector=sector_labels_vector[s];
          
          // Checking if sector is target sector
          bool allowed=sector.operator_labels()==u3shell::OperatorLabelsU3S(tensor_u3st.operator_labels());
          allowed&=sector.kappa0()==kappa0;
          allowed&=(sector.L0()==L0);
          if(not allowed)
              continue;
          
          //get subspace labels
          const spncci::SubspaceU3S& ket_subspace=target_space.GetSubspace(sector.ket_index());
          const spncci::SubspaceU3S& bra_subspace=target_space.GetSubspace(sector.bra_index());
          const u3::U3& omegap=bra_subspace.GetSubspaceLabels().U3();
          const u3::U3& omega=ket_subspace.GetSubspaceLabels().U3();
          int rho0=sector.rho0();

          // Iterating through U3 subspace "states" which correspond to baby spncci subspaces.
          for(int i=0; i<bra_subspace.size(); ++i)
            for(int j=0; j<ket_subspace.size(); ++j)
              {
                // (indexp,index)-> position of upper left corner of subsector in sector
                int indexp=bra_subspace.sector_index(i);
                int index=ket_subspace.sector_index(j);

                // extracting baby spncci information
                int baby_spncci_index_bra,baby_spncci_index_ket;
                std::tie(baby_spncci_index_bra)=bra_subspace.GetStateLabels(i);
                std::tie(baby_spncci_index_ket)=ket_subspace.GetStateLabels(j);
                
                const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
                  =baby_spncci_space.GetSubspace(baby_spncci_index_bra);
                const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
                  =baby_spncci_space.GetSubspace(baby_spncci_index_ket);
                
                int irrep_family_index_bra=baby_spncci_subspace_bra.irrep_family_index();
                int irrep_family_index_ket=baby_spncci_subspace_ket.irrep_family_index();
                
                const u3::U3& sigmap=baby_spncci_subspace_bra.sigma();
                const u3::U3& sigma=baby_spncci_subspace_ket.sigma();

                // (dimp,dim)->size of subsector
                int dimp=baby_spncci_subspace_bra.size();
                int dim=baby_spncci_subspace_ket.size();

                //Keys for looking up subsector in unit tensor cache
                std::pair<int,int> lgi_pair(irrep_family_index_bra,irrep_family_index_ket);
                std::pair<int,int> NnpNn(int(omegap.N()-sigmap.N()),int(omega.N()-sigma.N()));              
                spncci::UnitTensorU3Sector unit_sector(omegap,omega,tensor_u3st,rho0);
                
                // Get cache containing unit tensor sector 
                if(not unit_tensor_sector_cache.count(lgi_pair))
                  continue;
                if(not unit_tensor_sector_cache[lgi_pair].count(NnpNn)) 
                  continue;
                spncci::UnitTensorSectorsCache& cache=unit_tensor_sector_cache[lgi_pair][NnpNn];
                
                // // Diagonistic: check if expected unit tensor sectors are found
                // if(not cache.count(unit_sector))
                //   for(auto t=cache.begin(); t!=cache.end(); ++t)
                //     std::cout<<t->first.Str()<<std::endl;
  
                #pragma omp critical
                {
                  if(cache.count(unit_sector))
                    if(cache[unit_sector].cols()!=0)
                    {
                      matrix_vector[s].block(indexp,index,dimp,dim)+=interaction_rme*cache[unit_sector];
                      std::cout<<"unit sector "<<unit_sector.Str()<<std::endl;
                      std::cout<<cache[unit_sector]<<std::endl;
                    }
                }
              }
        }
      }
  }// end function

}  // namespace
