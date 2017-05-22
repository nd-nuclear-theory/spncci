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
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceU3S& target_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const std::vector<spncci::SectorLabelsU3S>& target_sectors_u3s,
      basis::OperatorBlocks<double>& target_blocks_u3s
    )
  {
    // Initial target sectors to zero matrices
    // target_blocks_u3s.resize(target_sectors_u3s.size());
    // for(int s=0; s<target_sectors_u3s.size(); ++s)
    //   {
    //     const spncci::SectorLabelsU3S& sector=target_sectors_u3s[s];
    //     const spncci::SubspaceU3S& ket_subspace=target_space.GetSubspace(sector.ket_index());
    //     const spncci::SubspaceU3S& bra_subspace=target_space.GetSubspace(sector.bra_index());
    //     int sector_dim_bra=bra_subspace.sector_dim();
    //     int sector_dim_ket=ket_subspace.sector_dim();
    //     target_blocks_u3s.emplace_back(Eigen::MatrixXd::Zero(sector_dim_bra,sector_dim_ket));
    //   }

    // iterate over interaction get unit tensor family index,kappa0,L0
    // iterate over U3S sectors to get target sectors
    // get corresponding unit tensor sectors 
    // Contract interaction with unit tensor blocks and accumulate in U3S blocks.   
    // std::cout<<"iterating over interation"<<std::endl;
    for(auto it=relative_observable.begin(); it!=relative_observable.end(); ++it)
      {
        int unit_tensor_subspace_index,kappa0,L0;
        std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;
        // std::cout<<"hi "<<unit_tensor_subspace_index<<"  "<<kappa0<<"  "<<L0<<std::endl;

        const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
          =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        const std::vector<double>& relative_rmes=it->second;

        int etap,eta;
        u3::SU3 x0; 
        HalfInt S0;
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
        // only the N0<=0 unit tensor hypersectors are compute.  
        // The N0>0 hypersectors are obtained via conjugation.
        bool conjugate_hypersector=((etap-eta)>0)?true:false;
        int unit_tensor_subspace_index_conj=-1;
        
        if(conjugate_hypersector)
        {
          u3shell::UnitTensorSubspaceLabels unit_tensor_labels_conj(u3::Conjugate(x0),S0,eta,etap);
          unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels_conj);
        }

        // OpenMP parallelize for look over target sectors to avoid race conditions 
        // std::cout<<"for each target sector "<<std::endl;
        for(int target_sector_index=0; target_sector_index<target_sectors_u3s.size(); ++target_sector_index)
          {
            // std::cout<<"target sector "<<target_sector_index<<"  of  "<<target_sectors_u3s.size()<<std::endl;
            const auto& target_sector=target_sectors_u3s[target_sector_index];
            bool allowed=(target_sector.x0()==x0);
            allowed&=(target_sector.S0()==S0);
            allowed&=((etap-eta)==target_sector.N0());
            allowed&=(target_sector.kappa0()==kappa0);
            allowed&=(target_sector.L0()==L0);
            if(not allowed)
              continue;

            const spncci::SubspaceU3S& ket_subspace=target_space.GetSubspace(target_sector.ket_index());
            const spncci::SubspaceU3S& bra_subspace=target_space.GetSubspace(target_sector.bra_index());
            int rho0=target_sector.rho0();
            
            // the states in U3S subspaces are baby spncci subspaces
            // std::cout<<"for each baby spncci bra and ket "<<std::endl;
            for(int bra_state_index=0; bra_state_index<bra_subspace.size(); ++bra_state_index)
              for(int ket_state_index=0; ket_state_index<ket_subspace.size(); ++ket_state_index)
                {
                  // index in courser grain u3s block
                  int block_index_u3s_bra=bra_subspace.sector_index(bra_state_index);
                  int block_index_u3s_ket=ket_subspace.sector_index(ket_state_index);

                  // extracting baby spncci information
                  int baby_spncci_index_bra,baby_spncci_index_ket;
                  std::tie(baby_spncci_index_bra)=bra_subspace.GetStateLabels(bra_state_index);
                  std::tie(baby_spncci_index_ket)=ket_subspace.GetStateLabels(ket_state_index);
                
                  const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
                          =baby_spncci_space.GetSubspace(baby_spncci_index_bra);
                  const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
                          =baby_spncci_space.GetSubspace(baby_spncci_index_ket);

                  // (dimp,dim)->size of baby spncci hyperblock
                  int dimp=baby_spncci_subspace_bra.size();
                  int dim=baby_spncci_subspace_ket.size();

                  
                  int baby_spncci_hypersector_index
                        =conjugate_hypersector?
                        baby_spncci_hypersectors.LookUpHypersectorIndex(
                            baby_spncci_index_ket,baby_spncci_index_bra,
                            unit_tensor_subspace_index_conj,rho0
                          ):                        
                        baby_spncci_hypersectors.LookUpHypersectorIndex(
                            baby_spncci_index_bra,baby_spncci_index_ket, 
                            unit_tensor_subspace_index,rho0
                          );
                  // if(baby_spncci_hypersector_index!=-1)
                  //   std::cout<<"hyper sector index "<<baby_spncci_hypersector_index<<std::endl;
                  // assert(baby_spncci_hypersector_index!=-1);
                  if(baby_spncci_hypersector_index==-1)
                    continue;

                  const basis::OperatorBlocks<double>& unit_tensor_blocks
                      =unit_tensor_hyperblocks[baby_spncci_hypersector_index];

                  // std::cout<<"for each unit tensor "<<std::endl;
                  for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
                  {
                    if(conjugate_hypersector)
                      {
                        int T0, S,T,Sp,Tp;
                        std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(unit_tensor_index);
                        std::tuple<int,int,int,int,int> conjugate_state(T0,S,T,Sp,Tp);
                        int unit_tensor_index_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj).LookUpStateIndex(conjugate_state);
                        
                        u3::U3 omegap(baby_spncci_subspace_bra.omega());
                        u3::U3 omega(baby_spncci_subspace_ket.omega());

                        int conjugation_grade=ParitySign(
                          u3::ConjugationGrade(omega)
                          +baby_spncci_subspace_ket.S()
                          +u3::ConjugationGrade(omegap)
                          +baby_spncci_subspace_bra.S()
                          );
                        
                        double conjugation_factor
                                =sqrt(
                                  1.*u3::dim(u3::SU3(etap,0))*u3::dim(omega)
                                  *am::dim(baby_spncci_subspace_ket.S())
                                  *am::dim(Sp)*am::dim(Tp)
                                  /u3::dim(u3::SU3(eta,0))/u3::dim(omegap)
                                  /am::dim(baby_spncci_subspace_bra.S())
                                  /am::dim(S)/am::dim(T)
                                  );

                        target_blocks_u3s[target_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim)
                          +=conjugation_grade*conjugation_factor
                            *relative_rmes[unit_tensor_index]
                            *unit_tensor_blocks[unit_tensor_index_conj].transpose();

                      }

                    else
                      target_blocks_u3s[target_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim)
                        +=relative_rmes[unit_tensor_index]*unit_tensor_blocks[unit_tensor_index];
                  }
                  // std::cout<<"finished sector"<<std::endl;
                }
              // std::cout<<"finshed baby spncci "<<std::endl;
          }
        // std::cout<<"finished target sectors "<<std::endl;
      }
      //   // Iterate over U3 sectors to get target sectors
      //   #pragma omp parallel for schedule(runtime)
      //   for(int s=0; s<sector_labels_vector.size(); ++s)
      //     {
      //       const spncci::SectorLabelsU3S& sector=sector_labels_vector[s];
          
      //       // Checking if sector is target sector
      //       bool allowed=sector.operator_labels()==u3shell::OperatorLabelsU3S(tensor_u3st.operator_labels());
      //       allowed&=sector.kappa0()==kappa0;
      //       allowed&=(sector.L0()==L0);
      //       if(not allowed)
      //         continue;
          
      //       //get subspace labels
      //       const spncci::SubspaceU3S& ket_subspace=target_space.GetSubspace(sector.ket_index());
      //       const spncci::SubspaceU3S& bra_subspace=target_space.GetSubspace(sector.bra_index());
      //       const u3::U3& omegap=bra_subspace.GetSubspaceLabels().U3();
      //       const u3::U3& omega=ket_subspace.GetSubspaceLabels().U3();
      //       int rho0=sector.rho0();

      //       // Iterating through U3 subspace "states" which correspond to baby spncci subspaces.
      //       for(int i=0; i<bra_subspace.size(); ++i)
      //         for(int j=0; j<ket_subspace.size(); ++j)
      //           {
      //             // (indexp,index)-> position of upper left corner of subsector in sector
      //             int indexp=bra_subspace.sector_index(i);
      //             int index=ket_subspace.sector_index(j);

      //             // extracting baby spncci information
      //             int baby_spncci_index_bra,baby_spncci_index_ket;
      //             std::tie(baby_spncci_index_bra)=bra_subspace.GetStateLabels(i);
      //             std::tie(baby_spncci_index_ket)=ket_subspace.GetStateLabels(j);
                
      //             const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
      //               =baby_spncci_space.GetSubspace(baby_spncci_index_bra);
      //             const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
      //               =baby_spncci_space.GetSubspace(baby_spncci_index_ket);
                
      //             int irrep_family_index_bra=baby_spncci_subspace_bra.irrep_family_index();
      //             int irrep_family_index_ket=baby_spncci_subspace_ket.irrep_family_index();
                
      //             const u3::U3& sigmap=baby_spncci_subspace_bra.sigma();
      //             const u3::U3& sigma=baby_spncci_subspace_ket.sigma();

      //             // (dimp,dim)->size of subsector
      //             int dimp=baby_spncci_subspace_bra.size();
      //             int dim=baby_spncci_subspace_ket.size();

      //             //Keys for looking up subsector in unit tensor cache
      //             std::pair<int,int> lgi_pair(irrep_family_index_bra,irrep_family_index_ket);
      //             std::pair<int,int> NnpNn(int(omegap.N()-sigmap.N()),int(omega.N()-sigma.N()));              
      //             spncci::UnitTensorU3Sector unit_sector(omegap,omega,tensor_u3st,rho0);
                
      //             // Get cache containing unit tensor sector 
      //             if(not unit_tensor_sector_cache.count(lgi_pair))
      //               continue;
      //             if(not unit_tensor_sector_cache.at(lgi_pair).count(NnpNn)) 
      //               continue;
      //             const spncci::UnitTensorSectorsCache& cache=unit_tensor_sector_cache.at(lgi_pair).at(NnpNn);
                
      //             // // Diagonistic: check if expected unit tensor sectors are found
      //             // if(not cache.count(unit_sector))
      //             //   for(auto t=cache.begin(); t!=cache.end(); ++t)
      //             //     std::cout<<t->first.Str()<<std::endl;
  
      //             // #pragma omp critical
      //             {
      //               if(cache.count(unit_sector))
      //                 if(cache.at(unit_sector).cols()!=0)
      //                   {
      //                     matrix_vector[s].block(indexp,index,dimp,dim)+=interaction_rme*cache.at(unit_sector);
      //                     // std::cout<<"unit sector "<<unit_sector.Str()<<std::endl;
      //                     // std::cout<<cache[unit_sector]<<std::endl;
      //                   }
      //             }
      //           }
      //     }
      // }
    std::cout<<"exiting contraction "<<std::endl;
  }// end function

}  // namespace
