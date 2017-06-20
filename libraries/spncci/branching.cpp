/****************************************************************
  branching.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching.h"

#include <set>
#include <sstream>
#include <iostream>

#include "cppformat/format.h"
#include "mcutils/parsing.h"
#include "spncci/unit_tensor.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to U3S level
  ////////////////////////////////////////////////////////////////  

  SubspaceSpU3S::SubspaceSpU3S(const u3::U3S& omegaS,const BabySpNCCISpace& baby_spncci_space)
  {

    // save labels
    labels_ = omegaS;

    // scan BabySpNCCISpace for states to accumulate
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace = baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // short circuit if subspace not relevant
        if (!(omegaS==baby_spncci_subspace.omegaS()))
          continue;

        // push state
        PushStateLabels(baby_spncci_subspace.sigmaSPN(),baby_spncci_subspace.size());

        // record auxiliary state information
        state_gamma_max_.push_back(baby_spncci_subspace.gamma_max());
        state_baby_spncci_subspace_index_.push_back(baby_spncci_subspace_index);
        state_irrep_family_index_.push_back(baby_spncci_subspace.irrep_family_index());
      }

  }

  std::string SubspaceSpU3S::LabelStr() const
  {
    return omegaS().Str();
  }

  std::string SubspaceSpU3S::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpU3S state(*this,state_index);

        os << fmt::format(
            "  index {} omegaS {} sigmaSPN {} degeneracy {} offset {}",
            state_index,
            state.omegaS().Str(),state.sigmaSPN().Str(),
            state.degeneracy(),state.offset()
          ) << std::endl;
      }

    return os.str();
  }

  SpaceSpU3S::SpaceSpU3S(const BabySpNCCISpace& baby_spncci_space)
  {
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // create new subspace -- only if not already constructed for this (omega,S)
        u3::U3S omegaS = u3::U3S(baby_spncci_subspace.omega(),baby_spncci_subspace.S());
        if(ContainsSubspace(omegaS))
          continue;
        PushSubspace(SubspaceSpU3S(omegaS,baby_spncci_space));
      }
  }

  std::string SpaceSpU3S::DebugStr(bool show_subspaces) const
  {
    std::ostringstream os;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // set up alias
        const SubspaceSpU3S& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "index {} omegaS {} size {} full_dimension {}",
            subspace_index,
            subspace.omegaS().Str(),subspace.size(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }


  std::string SectorLabelsSpU3S::Str() const
  {
    return fmt::format("( {} {}  {}{} {} : {} {}  {}", 
                       bra_index(),ket_index(), N0(), x0().Str(),S0(),kappa0(),L0(),rho0());
  }

  void GetSectorsSpU3S(
      const spncci::SpaceSpU3S& space, 
      const std::vector<u3shell::IndexedOperatorLabelsU3S>& relative_tensor_labels,
      std::vector<spncci::SectorLabelsSpU3S>& sector_vector
    )
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
            u3::U3S omegapSp(space.GetSubspace(i).labels());
            u3::U3S omegaS(space.GetSubspace(j).labels());
            int rho0_max=u3::OuterMultiplicity(omegaS.SU3(), op_labels.x0(),omegapSp.SU3());
            // Check if allowed U(1) coupling
            if(omegaS.U3().N()+op_labels.N0()!=omegapSp.U3().N())
              continue;
            // Check if allowed SU(2) coupling 
            if(not am::AllowedTriangle(omegaS.S(), op_labels.S0(), omegapSp.S()))
              continue;
            for(int rho0=1; rho0<=rho0_max; ++rho0)
              {
                spncci::SectorLabelsSpU3S sector(i,j,op_labels, kappa0,L0,rho0);
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






  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to LS level
  ////////////////////////////////////////////////////////////////  

  SubspaceSpLS::SubspaceSpLS(const spncci::LSLabels& ls_labels, const SpaceSpU3S& spu3s_space)
  {

    // save labels
    labels_ = ls_labels;

    // scan SpU3S space for states to accumulate
    for(int spu3s_subspace_index=0; spu3s_subspace_index<spu3s_space.size(); ++spu3s_subspace_index)
      {

        // set up alias
        const SubspaceSpU3S& spu3s_subspace = spu3s_space.GetSubspace(spu3s_subspace_index);

        // determine branching multiplicity to the specified L
        int kappa_max=u3::BranchingMultiplicitySO3(spu3s_subspace.x(),L());

        // short circuit if U3S subspace not relevant to current LS subspace
        if ((kappa_max==0)||(S()!=spu3s_subspace.S()))
          continue;

        // accumulate LS-branched states
        //
        // To obtain the mandated state ordering, we must pass through
        // the U3S subspace repeatedly, with kappa as the outermost
        // index.
        for (int kappa=1; kappa<=kappa_max; ++kappa)
          for (int spu3s_state_index=0; spu3s_state_index<spu3s_subspace.size(); ++spu3s_state_index)
            {
              // push state
              StateSpU3S spu3s_state(spu3s_subspace,spu3s_state_index);
              StateLabelsSpLS spls_state_labels = StateLabelsSpLS(
                  spu3s_state.omega(),
                  kappa,
                  spu3s_state.sigmaSPN()
                );
              PushStateLabels(spls_state_labels,spu3s_state.degeneracy());

              // record auxiliary state information
              state_gamma_max_.push_back(spu3s_state.gamma_max());
              state_baby_spncci_subspace_index_.push_back(spu3s_state.baby_spncci_subspace_index());
              state_spu3s_subspace_index_.push_back(spu3s_subspace_index);
            }

      }
  }

  std::string SubspaceSpLS::LabelStr() const
  {
    return fmt::format("({},{})",L(),S().Str());
  }

  std::string SubspaceSpLS::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpLS state(*this,state_index);

        os << fmt::format(
            "  index {} omegaS {} kappa {} sigmaSPN {} degeneracy {} offset {}",
            state_index,
            state.omegaS().Str(),state.kappa(),state.sigmaSPN().Str(),
            state.degeneracy(),state.offset()
          ) << std::endl;
      }

    return os.str();
  }


  SpaceSpLS::SpaceSpLS(const SpaceSpU3S& spu3s_space)
  {

    // collect (L,S) branchings
    std::set<LSLabels> ls_labels_set;
    for(int spu3s_subspace_index=0; spu3s_subspace_index<spu3s_space.size(); ++spu3s_subspace_index)
      {

        // set up alias
        const SubspaceSpU3S& spu3s_subspace = spu3s_space.GetSubspace(spu3s_subspace_index);

        // find branching L values
        u3::SU3 x = spu3s_subspace.omega().SU3();
        HalfInt S = spu3s_subspace.S();
        MultiplicityTagged<int>::vector branching = u3::BranchingSO3(x);

        // accumulate (L,S) pairs from U3S subspace
        for (const MultiplicityTagged<int>& l_kappa_max : branching)
          {
            int L = l_kappa_max.irrep;
            ls_labels_set.insert(LSLabels(L,S));
          }
      }
    
    // create subspaces
    for (const LSLabels& ls_labels : ls_labels_set)
      PushSubspace(SubspaceSpLS(ls_labels,spu3s_space));
  }

  SpaceSpLS::SpaceSpLS(const SpaceSpU3S& spu3s_space, HalfInt J)
  {

    // collect (L,S) branchings
    std::set<LSLabels> ls_labels_set;
    for(int spu3s_subspace_index=0; spu3s_subspace_index<spu3s_space.size(); ++spu3s_subspace_index)
      {

        // set up alias
        const SubspaceSpU3S& spu3s_subspace = spu3s_space.GetSubspace(spu3s_subspace_index);

        // find branching L values
        u3::SU3 x = spu3s_subspace.omega().SU3();
        HalfInt S = spu3s_subspace.S();
        HalfInt::pair l_range = am::ProductAngularMomentumRange(J,S);
        MultiplicityTagged<int>::vector branching = u3::BranchingSO3Constrained(x,l_range);

        // accumulate (L,S) pairs from U3S subspace
        for (const MultiplicityTagged<int>& l_kappa_max : branching)
          {
            int L = l_kappa_max.irrep;
            ls_labels_set.insert(LSLabels(L,S));
          }
      }
    
    // create subspaces
    for (const LSLabels& ls_labels : ls_labels_set)
      PushSubspace(SubspaceSpLS(ls_labels,spu3s_space));
  }

  std::string SpaceSpLS::DebugStr(bool show_subspaces) const
  {
    std::ostringstream os;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // set up alias
        const SubspaceSpLS& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "index {} (L,S) ({},{}) size {} full_dimension {}",
            subspace_index,
            subspace.L(),
            subspace.S().Str(),subspace.size(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }


  // void 
  // ContractAndRegroupU3S(
  //     const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  //     const spncci::BabySpNCCISpace& baby_spncci_space,
  //     const spncci::SpaceSpU3S& target_space,
  //     const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
  //     const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
  //     const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
  //     const std::vector<spncci::SectorLabelsSpU3S>& target_sectors_u3s,
  //     basis::OperatorBlocks<double>& target_blocks_u3s
  //   )
  // {
  //   for(int target_sector_index=0; target_sector_index<target_sectors_u3s.size(); ++target_sector_index)
  //     {
  //       int target_bra_index,target_ket_index,target_kappa0,target_L0,target_rho0;
  //       u3shell::OperatorLabelsU3S target_tensor_labels;
  //       std::tie(target_bra_index,target_ket_index,target_tensor_labels,target_kappa0,target_L0,target_rho0)
  //         =target_sectors_u3s[target_sector_index].Key();
        
  //       // Get corresponding SpU3S subspaces 
  //       const spncci::SubspaceSpU3S& bra_target_subspace=target_space.GetSubspace(target_bra_index);
  //       const spncci::SubspaceSpU3S& ket_target_subspace=target_space.GetSubspace(target_ket_index);

  //       u3::U3 omegap=bra_target_subspace.omega();
  //       u3::U3 omega=ket_target_subspace.omega();

  //       for(auto it=relative_observable.begin(); it!=relative_observable.end(); ++it)
  //         {
  //           int unit_tensor_subspace_index,kappa0,L0;
  //           std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;

  //           if((target_L0!=L0)||(target_kappa0!=kappa0))
  //             continue;

  //           // Get rmes 
  //           const std::vector<double>& relative_rmes=it->second;
            
  //           // Get Unit Tensor subspace 
  //           const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
  //             =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

  //           u3::SU3 x0;
  //           int etap,eta;
  //           HalfInt S0;
  //           std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();

  //           // Check if operator tensor labels match target 
  //           if( (not (target_tensor_labels.x0()==x0))
  //             || (target_tensor_labels.S0()!=S0) 
  //             || (target_tensor_labels.N0()!=(etap-eta))
  //             )
  //             continue;

  //           for(int bra_state_index=0; bra_state_index<bra_target_subspace.size(); ++bra_state_index)
  //             for(int ket_state_index=0; ket_state_index<ket_target_subspace.size(); ++ket_state_index)
  //               {
  //                 // // index in courser grain u3s block
  //                 // int block_index_u3s_bra=bra_target_subspace.sector_index(bra_state_index);
  //                 // int block_index_u3s_ket=ket_target_subspace.sector_index(ket_state_index);

  //                 // extracting baby spncci information
  //                 int baby_spncci_index_bra,baby_spncci_index_ket;
  //                 auto& baby_spncci_bra=bra_target_subspace.GetStateLabels(bra_state_index);
  //                 auto& baby_spncci_ket=ket_target_subspace.GetStateLabels(ket_state_index);
                
  //                 // const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
  //                 //         =baby_spncci_space.GetSubspace(baby_spncci_index_bra);
  //                 // const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
  //                 //         =baby_spncci_space.GetSubspace(baby_spncci_index_ket);

  //                 // (dimp,dim)->size of baby spncci hyperblock
  //                 int dimp=baby_spncci_bra.degeneracy();
  //                 int dim=baby_spncci_ket.degeneracy();

  //                 // offset of baby spncci brat and ket in courser grain u3s block
  //                 int block_index_u3s_bra=baby_spncci_bra.offset();
  //                 int block_index_u3s_ket=baby_spncci_ket.offset();

  //                 // Recurrence computes Nnp<=Nn sectors for lgi_bra<=lgi_ket and Nnp>Nn sectors for lgi_bra>lgi_ket.
  //                 // The remaining unit tensor blocks are obtained by conjugation, i.e., those satisfying
  //                 //    If Nnp-Nn>0
  //                 //    If Nnp-Nn=0, lgi_bra>lgi_ket
  //                 // 
  //                 int Nnp=baby_spncci_bra.Nn();
  //                 int Nn=baby_spncci_ket.Nn();
  //                 int unit_tensor_subspace_index_conj=-1;
                  
  //                 // Check we need hypersector or conjugate hypersector
  //                 bool conjugate_hypersector;
                  
  //                 if( ((Nnp-Nn)==0)
  //                     &&(baby_spncci_bra.irrep_family_index()>baby_spncci_ket.irrep_family_index())
  //                   )
  //                   conjugate_hypersector=true;

  //                 if(conjugate_hypersector)
  //                   {
  //                     u3shell::UnitTensorSubspaceLabels unit_tensor_labels_conj(u3::Conjugate(x0),S0,eta,etap);
  //                     unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels_conj);
  //                   }

  //                 int baby_spncci_hypersector_index
  //                       =conjugate_hypersector?
  //                       baby_spncci_hypersectors.LookUpHypersectorIndex(
  //                           baby_spncci_index_ket,baby_spncci_index_bra,
  //                           unit_tensor_subspace_index_conj,target_rho0
  //                         ):                        
  //                       baby_spncci_hypersectors.LookUpHypersectorIndex(
  //                           baby_spncci_index_bra,baby_spncci_index_ket, 
  //                           unit_tensor_subspace_index,target_rho0
  //                         );

  //                 // If baby_spncci_hypersector doesn't exist, continue
  //                 if(baby_spncci_hypersector_index==-1)
  //                   continue;

  //                 // Get unit tensor hyperblocks 
  //                 const basis::OperatorBlocks<double>& unit_tensor_blocks
  //                     =unit_tensor_hyperblocks[baby_spncci_hypersector_index];

  //                 // Loop through unit tensors in unit tensor subspace, and accumulate product with relative rmes 
  //                 // in target sector.  
  //                 for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
  //                 {
  //                   if(conjugate_hypersector)
  //                     {
  //                       int T0, S,T,Sp,Tp;
  //                       std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(unit_tensor_index);
  //                       std::tuple<int,int,int,int,int> conjugate_state(T0,S,T,Sp,Tp);
  //                       int unit_tensor_index_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj).LookUpStateIndex(conjugate_state);
                        
  //                       u3::U3 omegap(baby_spncci_subspace_bra.omega());
  //                       u3::U3 omega(baby_spncci_subspace_ket.omega());

  //                       int conjugation_grade=ParitySign(
  //                         u3::ConjugationGrade(omega)
  //                         +baby_spncci_subspace_ket.S()
  //                         +u3::ConjugationGrade(omegap)
  //                         +baby_spncci_subspace_bra.S()
  //                         );
                        
  //                       double conjugation_factor
  //                               =sqrt(
  //                                 1.*u3::dim(u3::SU3(etap,0))*u3::dim(omega)
  //                                 *am::dim(baby_spncci_ket.S())
  //                                 *am::dim(Sp)*am::dim(Tp)
  //                                 /u3::dim(u3::SU3(eta,0))/u3::dim(omegap)
  //                                 /am::dim(baby_spncci_bra.S())
  //                                 /am::dim(S)/am::dim(T)
  //                                 );

  //                       target_blocks_u3s[target_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim)
  //                         +=conjugation_grade*conjugation_factor
  //                           *relative_rmes[unit_tensor_index]
  //                           *unit_tensor_blocks[unit_tensor_index_conj].transpose();

  //                     }

  //                   else
  //                     target_blocks_u3s[target_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim)
  //                       +=relative_rmes[unit_tensor_index]*unit_tensor_blocks[unit_tensor_index];
  //                 }
  //               }
  //         }
        
  //     }
  // }



}  // namespace
