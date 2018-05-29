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
#include "mcutils/eigen.h"
#include "spncci/unit_tensor.h"

namespace spncci
{

  SubspaceU3S::SubspaceU3S(const u3::U3S& omegaS,const BabySpNCCISpace& baby_spncci_space)
  {

    // save labels
    labels_ = omegaS;

    // scan BabySpNCCISpace for states to accumulate
    int substate_offset = 0;  // accumulated offset
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace = baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // short circuit if subspace not relevant
        if (!(omegaS==baby_spncci_subspace.omegaS()))
          continue;

        // push state
        PushStateLabels(StateLabelsType(baby_spncci_subspace_index));

        // record state multiplicity indexing information
        state_substate_offset_.push_back(substate_offset);
        int state_dimension = baby_spncci_subspace.size();
        state_dimension_.push_back(state_dimension);

        // store state symplectic irrep information
        state_sigmaSPN_.push_back(baby_spncci_subspace.omegaSPN());
        state_gamma_max_.push_back(baby_spncci_subspace.gamma_max());

        // accumulate offset for next state
        substate_offset += state_dimension;
        
      }

    // store final full dimension
    full_dimension_ = substate_offset;

  }

  std::string SubspaceU3S::Str() const
  {
    return fmt::format("{}",labels_.Str());
  }

  SpaceU3S::SpaceU3S(const BabySpNCCISpace& baby_spncci_space)
  {
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // create new subspace -- only if not already constructed for this (omega,S)
        u3::U3S omegaS(baby_spncci_subspace.omega(),baby_spncci_subspace.S());
        if(ContainsSubspace(omegaS))
          continue;
        SubspaceU3S subspace(omegaS,baby_spncci_space);
        PushSubspace(subspace);
      }
  }


  ObservableHypersectorsU3S::ObservableHypersectorsU3S(
    const spncci::SpaceU3S& space_u3s,
    const u3shell::ObservableSpaceU3S& observable_space
  )
  {
    int hypersector_index=0;
    for (int bra_subspace_index=0; bra_subspace_index<space_u3s.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space_u3s.size(); ++ket_subspace_index)
        {
          // retrieve subspaces
          const SubspaceU3S& bra_subspace = space_u3s.GetSubspace(bra_subspace_index);
          const SubspaceU3S& ket_subspace = space_u3s.GetSubspace(ket_subspace_index);

          // For each observable subspace, check if its an allowed observable subspace determined
          // by SU(2) and U(3) constraints.  If allowed, push multiplicity tagged hypersectors
          for(int observable_subspace_index=0; observable_subspace_index<observable_space.size(); ++observable_subspace_index)
            {
              bool allowed_subspace = true;
              const u3shell::ObservableSubspaceU3S& 
                observable_subspace=observable_space.GetSubspace(observable_subspace_index);

              // U(1)
              allowed_subspace
                &=(ket_subspace.omega().N()+observable_subspace.N0()-bra_subspace.omega().N() == 0);
              
              // spin
              //
              // Note: Basic two-body constaints can be placed on Sp
              // and Sn triangularity based on two-body nature of
              // observable, so (delta Sp)<=2 and (delta Sn)<=2.  However, in
              // general, the observable does not have sharp Sp0 or Sn0.
              allowed_subspace &= am::AllowedTriangle(ket_subspace.S(),observable_subspace.S0(),bra_subspace.S());
              // allowed_subspace &= abs(int(ket_subspace.Sp()-bra_subspace.Sp()))<=2;
              // allowed_subspace &= abs(int(ket_subspace.Sn()-bra_subspace.Sn()))<=2;
              if (!allowed_subspace)
                continue;

              // find SU(3) multiplicity and check SU(3) selection
              int multiplicity = u3::OuterMultiplicity(
                  ket_subspace.omega().SU3(),
                  observable_subspace.x0(),
                  bra_subspace.omega().SU3()
                );

              // push sectors (tagged by multiplicity)
              for (int multiplicity_index = 1; multiplicity_index <= multiplicity; ++multiplicity_index)
                {
                  PushHypersector(
                    HypersectorType(
                      bra_subspace_index,ket_subspace_index,observable_subspace_index,
                      bra_subspace, ket_subspace,observable_subspace,
                      multiplicity_index
                      )
                    );
                }
            }
        }
  }


 void InitializeU3SObservableBlocks(
      const spncci::SpaceU3S& space_u3s,
      int num_observables, int num_hw_values,
      const spncci::ObservableHypersectorsU3S& hypersectors,
      spncci::OperatorBlocks& blocks
    )
  {
    blocks.resize(hypersectors.size());
    for(int hypersector_index=0; hypersector_index<hypersectors.size(); ++hypersector_index)
      {
        const auto& hypersector=hypersectors.GetHypersector(hypersector_index);
        int rows=space_u3s.GetSubspace(hypersector.bra_subspace_index()).full_dimension();
        int cols=space_u3s.GetSubspace(hypersector.ket_subspace_index()).full_dimension();
        blocks[hypersector_index]=spncci::OperatorBlock::Zero(rows,cols);
      }
  }




 void InitializeU3SObservableBlocks(
      const spncci::SpaceU3S& space_u3s,
      int num_observables, int num_hw_values,
      const std::vector<spncci::ObservableHypersectorsU3S>& observable_hypersectors_by_observable,
      std::vector<std::vector<spncci::OperatorBlocks>>& observables_blocks_array
    )
  {
    // vector of blocks for u3 sectors for each hbar omega,for each observable
    observables_blocks_array.resize(num_hw_values);

    // For each hbar omega, zero initialize block for each observable
    // based on basis::SetOperatorToZero in operator.h
    // int total_entries = 0;
    for(int hw_index=0; hw_index<num_hw_values; ++hw_index)
      {
        std::vector<spncci::OperatorBlocks>& observables_blocks=observables_blocks_array[hw_index];
        observables_blocks.resize(num_observables);

        for(int observable_index=0; observable_index<num_observables; ++observable_index)
          {
            spncci::OperatorBlocks& blocks=observables_blocks[observable_index];
            
            const spncci::ObservableHypersectorsU3S& hypersectors
                    =observable_hypersectors_by_observable[observable_index];
            
            blocks.resize(hypersectors.size());
            // std::cout<<"num blocks "<<blocks.size()<<std::endl;
            for(int hypersector_index=0; hypersector_index<hypersectors.size(); ++hypersector_index)
              {
                const auto& hypersector=hypersectors.GetHypersector(hypersector_index);
                int rows=space_u3s.GetSubspace(hypersector.bra_subspace_index()).full_dimension();
                int cols=space_u3s.GetSubspace(hypersector.ket_subspace_index()).full_dimension();
                blocks[hypersector_index]=spncci::OperatorBlock::Zero(rows,cols);
              }
          }
      } 
  }



void 
  ContractBabySpNCCIU3S(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      basis::OperatorHyperblocks<double>& observable_hyperblocks
    )
  {
    //DEPRECATED
    // iterate over interaction get unit tensor family index,kappa0,L0
    // iterate over U3S sectors to get observable sectors
    // get corresponding unit tensor sectors 
    // Contract interaction with unit tensor blocks and accumulate in U3S blocks.   
    // std::cout<<"iterating over interation"<<std::endl;

    for(auto it=relative_observable.begin(); it!=relative_observable.end(); ++it)
      {
        const std::vector<double>& relative_rmes=it->second;

        int unit_tensor_subspace_index,kappa0,L0;
        std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;
        // std::cout<<"hi "<<unit_tensor_subspace_index<<"  "<<kappa0<<"  "<<L0<<std::endl;

        const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
          =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        int etap,eta;
        u3::SU3 x0; 
        HalfInt S0;
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();

        // Look up observable index
        int observable_subspace_index=observable_space.LookUpSubspaceIndex(u3shell::ObservableSubspaceLabels(etap-eta,x0,S0,kappa0,L0));
        // std::cout<<"observable_subspace_index "<<observable_subspace_index<<std::endl;

        // std::cout<<"for each target sector "<<std::endl;
        for(int target_sector_index=0; target_sector_index<observable_hypersectors.size(); ++target_sector_index)
          {
            
            const auto& target_hypersector=observable_hypersectors.GetHypersector(target_sector_index);
            // std::cout<<"target sector "<<target_sector_index<<" "<<target_hypersector.operator_subspace_index()
            //     <<"  "<<observable_subspace_index<<"  "<<(target_hypersector.operator_subspace_index()!=observable_subspace_index)<<std::endl;
            
            if(target_hypersector.operator_subspace_index()!=observable_subspace_index)
              continue;

            // std::cout<<"Get target hypersector block"<<std::endl;
            basis::OperatorBlock<double>& target_block=observable_hyperblocks[target_sector_index][0];

            int baby_spncci_index_ket=target_hypersector.ket_subspace_index();
            int baby_spncci_index_bra=target_hypersector.bra_subspace_index();

             int rho0=target_hypersector.multiplicity_index();
            
            const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
                    =baby_spncci_space.GetSubspace(baby_spncci_index_bra);

            const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
                    =baby_spncci_space.GetSubspace(baby_spncci_index_ket);

             // Recurrence computes Nnp>=Nn sectors for lgi_bra=>lgi_ket.
            // The remaining unit tensor blocks are obtained by conjugation, i.e., those satisfying
            //    Nnp-Nn<0                  
            int Nnp=baby_spncci_subspace_bra.Nn();
            int Nn=baby_spncci_subspace_ket.Nn();
            int unit_tensor_subspace_index_conj=-1;
            bool conjugate_hypersector=(Nnp<Nn);

            // If we need to get conjugate hyper sector, then need conjugate labels 
            if(conjugate_hypersector)
              {
                u3shell::UnitTensorSubspaceLabels unit_tensor_labels_conj(u3::Conjugate(x0),S0,eta,etap);
                unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels_conj);
              }

            // Get baby_spncci_hypersector_index for hypersector or conjugate hypersector
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

            // std::cout<<"hypersector index "<<baby_spncci_hypersector_index<<std::endl;
            // If hypersector not found, continue
            if(baby_spncci_hypersector_index==-1)
              continue;

            const basis::OperatorBlocks<double>& unit_tensor_blocks
                =unit_tensor_hyperblocks[baby_spncci_hypersector_index];

            // std::cout<<"for each unit tensor "<<std::endl;
            for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
              {
                if(conjugate_hypersector)
                  {
                    // get unit tensor index for conjugate unit tensor
                    int T0, S,T,Sp,Tp;
                    std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(unit_tensor_index);
                    std::tuple<int,int,int,int,int> conjugate_state(T0,S,T,Sp,Tp);
                    int unit_tensor_index_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj).LookUpStateIndex(conjugate_state);
                    
                    u3::U3 omegap(baby_spncci_subspace_bra.omega());
                    u3::U3 omega(baby_spncci_subspace_ket.omega());

                    // Get conjugation grade
                    int conjugation_grade=ParitySign(
                      u3::ConjugationGrade(omega)
                      -baby_spncci_subspace_ket.S()
                      +u3::ConjugationGrade(omegap)
                      +baby_spncci_subspace_bra.S()
                      );
                    
                    // and conjugation factor
                    double conjugation_factor
                            =sqrt(
                              1.*u3::dim(u3::SU3(etap,0))*u3::dim(omega)
                              *am::dim(baby_spncci_subspace_ket.S())
                              *am::dim(Sp)*am::dim(Tp)
                              /u3::dim(u3::SU3(eta,0))/u3::dim(omegap)
                              /am::dim(baby_spncci_subspace_bra.S())
                              /am::dim(S)/am::dim(T)
                              );

                    target_block
                      +=conjugation_grade*conjugation_factor
                        *relative_rmes[unit_tensor_index]
                        *unit_tensor_blocks[unit_tensor_index_conj].transpose();
                    // std::cout<<target_block<<std::endl<<relative_rmes[unit_tensor_index]<<std::endl<<unit_tensor_blocks[unit_tensor_index_conj]<<std::endl;
                  }

                else
                {
                  target_block
                    +=relative_rmes[unit_tensor_index]*unit_tensor_blocks[unit_tensor_index];
                }
              
              }// unit tensor index
          } //target_sector_index
        // std::cout<<"finished target sectors "<<std::endl;
      }//relative observable
  }

  void 
  RegroupU3S(
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::SpaceU3S& u3s_space,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      const basis::OperatorHyperblocks<double>& observable_hyperblocks,
      const spncci::ObservableHypersectorsU3S& observable_sectors,
      spncci::OperatorBlocks& observable_blocks
    )

    {
      // std::cout<<"entering regroup"<<std::endl;
      // #pragma omp parallel for schedule(dynamic)
      for(int observable_hypersector_index=0; observable_hypersector_index<observable_hypersectors.size(); ++observable_hypersector_index)
        {
          // std::cout<<observable_hyperblocks[observable_hypersector_index][0]<<std::endl;
          // std::cout<<"Look up hypersector and extract bra, ket and operator subspace"<<std::endl;
          const auto& baby_spncci_hypersector=observable_hypersectors.GetHypersector(observable_hypersector_index);

          int baby_spncci_index_bra, baby_spncci_index_ket, operator_subspace_index, rho0;
          std::tie(baby_spncci_index_bra,baby_spncci_index_ket,operator_subspace_index,rho0)=baby_spncci_hypersector.Key();

          const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
          const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);
          
          // std::cout<<"Get labels to look up corresponding U3S subspaces"<<std::endl;
          const u3::U3& omegap=baby_spncci_subspace_bra.omega();
          HalfInt Sp=baby_spncci_subspace_bra.S();
          const u3::U3& omega=baby_spncci_subspace_ket.omega();
          HalfInt S=baby_spncci_subspace_ket.S();

          u3::U3S omegaS_bra(omegap,Sp);
          u3::U3S omegaS_ket(omega,S);

          // std::cout<<"Get U3S subspace and state indices"<<std::endl;
          int bra_subspace_index=u3s_space.LookUpSubspaceIndex(omegaS_bra);
          int ket_subspace_index=u3s_space.LookUpSubspaceIndex(omegaS_ket);

          const spncci::SubspaceU3S& bra_subspace=u3s_space.GetSubspace(bra_subspace_index);
          const spncci::SubspaceU3S& ket_subspace=u3s_space.GetSubspace(ket_subspace_index);

          int bra_state_index=bra_subspace.LookUpStateIndex(std::tuple<int>(baby_spncci_index_bra));
          int ket_state_index=ket_subspace.LookUpStateIndex(std::tuple<int>(baby_spncci_index_ket));

          int observable_sector_index
            =observable_sectors.LookUpHypersectorIndex(bra_subspace_index,ket_subspace_index,operator_subspace_index,rho0);

          // std::cout<<"Get index in courser grain u3s block"<<std::endl;
          int block_index_u3s_bra=bra_subspace.sector_index(bra_state_index);
          int block_index_u3s_ket=ket_subspace.sector_index(ket_state_index);

          // std::cout<<"Get (dimp,dim)->size of baby spncci hyperblock"<<std::endl;
          int dimp=baby_spncci_subspace_bra.size();
          int dim=baby_spncci_subspace_ket.size();

          // std::cout<<"Regroup "<<observable_sector_index<<"  "<<observable_hypersector_index<<"     "
          // <<block_index_u3s_bra<<"  "<<block_index_u3s_ket<<"  "<<dimp<<"  "<<dim<<std::endl;
          auto matrix_block=observable_blocks[observable_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim);

          // std::cout<<observable_blocks[observable_sector_index]<<std::endl<<matrix_block<<std::endl<<std::endl;
          observable_blocks[observable_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim)
            =observable_hyperblocks[observable_hypersector_index][0];
          // std::cout<<"hi"<<std::endl;
        }//observable hypersector index
    }


  void 
  ContractAndRegroupU3S(
      const std::pair<int,int>& lgi_pair,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::SpaceU3S& u3s_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableHypersectorsU3S& observable_sectors,
      spncci::OperatorBlocks& observable_blocks
    )
  {
    // DEPRECATED

    // generate hypersectors ...
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
    spncci::ObservableBabySpNCCIHypersectors observable_hypersectors(
        baby_spncci_space,observable_space,
        irrep_family_index_bra,irrep_family_index_ket
      );
    
    basis::OperatorHyperblocks<double> observable_hyperblocks;
    basis::SetHyperoperatorToZero(observable_hypersectors,observable_hyperblocks);
    
    // Contract over baby spnci observable sectors 
    spncci::ContractBabySpNCCIU3S(
      unit_tensor_space,baby_spncci_space,observable_space,
      relative_observable,baby_spncci_hypersectors,
      unit_tensor_hyperblocks,observable_hypersectors,
      observable_hyperblocks
    );

    // Regroup into u3s sectors
    spncci::RegroupU3S(
        baby_spncci_space,observable_space,u3s_space,
        observable_hypersectors,observable_hyperblocks,
        observable_sectors,observable_blocks
      );
  }// end function



void 
  ContractBabySpNCCIU3S2(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      basis::OperatorHyperblocks<double>& observable_hyperblocks
    )
  {
    // iterate over interaction get unit tensor family index,kappa0,L0
    // iterate over baby spncci unit tensors
    // Contract into corresponding observable hyperblock
    // Contract conjugated unit tensor into corresponding hyperblock
    for(auto it=relative_observable.begin(); it!=relative_observable.end(); ++it)
      {
        const std::vector<double>& relative_rmes=it->second;

        int unit_tensor_subspace_index,kappa0,L0;
        std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;
        // std::cout<<"hi "<<unit_tensor_subspace_index<<"  "<<kappa0<<"  "<<L0<<std::endl;

        const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        int etap,eta;
        u3::SU3 x0; 
        HalfInt S0;
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();

        // Look up observable index
        int observable_subspace_index=observable_space.LookUpSubspaceIndex(u3shell::ObservableSubspaceLabels(etap-eta,x0,S0,kappa0,L0));
        // std::cout<<"observable_subspace_index "<<observable_subspace_index<<std::endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // std::cout<<"Look up conjugate unit tensor subspace index"<<std::endl;
        // std::cout<<u3::Conjugate(x0).Str()<<"  "<<S0<<"  "<<eta<<"  "<<etap<<std::endl;
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
        int unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);

        // std::cout<<"Get conjugate unit tensor subspace"<<std::endl;
        const auto& unit_tensor_subspace_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj);

        // std::cout<<"Get conjugate relative observables"<<std::endl;
        std::tuple<int,int,int> rme_labels_conj(unit_tensor_subspace_index_conj,kappa0,L0);
        const std::vector<double>& relative_rmes_conj=relative_observable.at(rme_labels_conj);

        // std::cout<<"Look up conjugate observable subspace index"<<std::endl;
        int observable_subspace_index_conj
          =observable_space.LookUpSubspaceIndex(u3shell::ObservableSubspaceLabels(eta-etap,u3::Conjugate(x0),S0,kappa0,L0));

        // std::cout<<"Iterate over baby spncci hypersectors."<<std::endl; 
        for(int baby_spncci_hypersector_index=0; baby_spncci_hypersector_index<baby_spncci_hypersectors.size(); ++baby_spncci_hypersector_index)
          {
            const auto& baby_spncci_hypersector=baby_spncci_hypersectors.GetHypersector(baby_spncci_hypersector_index);

            // std::cout<<"Check if hypersector contains desired unit tensor subspace"<<std::endl;
            if(baby_spncci_hypersector.operator_subspace_index()!=unit_tensor_subspace_index)
              continue;

            // std::cout<<"Get baby_spncci subspace indices"<<std::endl;
            int baby_spncci_index_bra=baby_spncci_hypersector.bra_subspace_index();
            int baby_spncci_index_ket=baby_spncci_hypersector.ket_subspace_index();
            int rho0=baby_spncci_hypersector.multiplicity_index();

            // std::cout<<"Look up baby spncci subspaces"<<std::endl;
            const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
            const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);

            // std::cout<<"Check if hypersector belons to diagonal sector"<<std::endl;
            // Note: recurrence computes Nnp>=Nn sectors for all (irrep1,irrep2).  Nn>Nnp sectors obtained from (irrep2,irrep1)
            // rmes by conjugation 
            int Nnp=baby_spncci_subspace_bra.Nn();
            int Nn=baby_spncci_subspace_ket.Nn();
            bool diagonal_hypersector=(Nnp==Nn);

            // Look up target observable hypersector and conjugate hypersector
            int observable_hypersector_index
                =observable_hypersectors.LookUpHypersectorIndex(baby_spncci_index_bra,baby_spncci_index_ket,observable_subspace_index,rho0);

            // std::cout<<"Get observable hyperblock and unit tensor hyperblock"<<std::endl;
            const basis::OperatorBlocks<double>& unit_tensor_blocks=unit_tensor_hyperblocks[baby_spncci_hypersector_index];
            basis::OperatorBlock<double>& observable_block=observable_hyperblocks[observable_hypersector_index][0];

            // std::cout<<"Loop over unit tensors and contract"<<std::endl;
            for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
                observable_block+=relative_rmes[unit_tensor_index]*unit_tensor_blocks[unit_tensor_index];                


            // If baby spncci hypersector is not diagonal, then we need to contract over conjugate unit tensor hypersectors to get Nn>Nnp sectors 
            if(not diagonal_hypersector)
              {
                int observable_hypersector_index_conj
                    =observable_hypersectors.LookUpHypersectorIndex(baby_spncci_index_ket,baby_spncci_index_bra,observable_subspace_index_conj,rho0);

                // std::cout<<"get conjugate block "<<observable_hypersector_index_conj<<std::endl;
                basis::OperatorBlock<double>& observable_block_conj
                    =observable_hyperblocks[observable_hypersector_index_conj][0];

                // std::cout<<"conjugation factor base"<<std::endl;
                // 
                // get bra and ket subspace labels
                const u3::U3& omegap=baby_spncci_subspace_bra.omega();
                const HalfInt& Sp=baby_spncci_subspace_bra.S();
                const u3::U3& omega=baby_spncci_subspace_ket.omega();
                const HalfInt& S=baby_spncci_subspace_ket.S();

                double conjugation_factor_base
                    =ParitySign(u3::ConjugationGrade(omegap)+Sp-u3::ConjugationGrade(omega)-S)
                      *sqrt(
                          1.*u3::dim(omegap)*am::dim(Sp)*u3::dim(u3::SU3(eta,0))
                          /u3::dim(omega)/am::dim(S)/u3::dim(u3::SU3(etap,0))
                        );
                // std::cout<<"iterate over untit tensors"<<std::endl;
                for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
                  {
                    // std::cout<<"Get state labels"<<std::endl;
                    int T0, S,T,Sp,Tp;
                    std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(unit_tensor_index);


                    // std::cout<<"Get conjugate state index"<<std::endl;
                    std::tuple<int,int,int,int,int> conjugate_state(T0,S,T,Sp,Tp);
                    int unit_tensor_index_conj=unit_tensor_subspace_conj.LookUpStateIndex(conjugate_state);
                    
                    // calculate conjugation factor
                    double conjugation_factor
                        =sqrt(am::dim(S)*am::dim(T)/am::dim(Sp)/am::dim(Tp))*conjugation_factor_base;

                    // std::cout<<unit_tensor_index_conj<<"  "<<unit_tensor_index_conj<<"  "<<unit_tensor_blocks.size()<<"  "<<relative_rmes_conj.size()<<std::endl;
                    observable_block_conj
                        +=conjugation_factor
                          *relative_rmes_conj[unit_tensor_index_conj]
                          *unit_tensor_blocks[unit_tensor_index_conj].transpose();
                  }

              }

          }
      }
  }


}  // namespace
