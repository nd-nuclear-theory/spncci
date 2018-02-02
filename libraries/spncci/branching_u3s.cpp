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

  std::string SectorLabelsU3S::Str() const
  {
    return fmt::format("( {} {}  {}{} {} : {} {}  {}", 
                       bra_index(),ket_index(), N0(), x0().Str(),S0(),kappa0(),L0(),rho0());
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
                spncci::SectorLabelsU3S sector(i,j,op_labels, kappa0,L0,rho0);
                // should not alread be a sector in map
                // assert(not u3s_sectors.count(sector));
                // if(not u3s_sectors.count(sector))
                //   {
                //     u3s_sectors[sector]=u3s_sector_vector_index;
                //     ++u3s_sector_vector_index;
                //   }
                // std::cout<<omegaS.Str()<<std::endl;
                
                sector_vector.push_back(sector);
              }
          }
    // sector_vector.resize(u3s_sectors.size());
    // for(auto it=u3s_sectors.begin(); it!=u3s_sectors.end(); ++it)
    //   sector_vector[it->second]=it->first;
  }

void InitializeU3SSectors(
      const spncci::SpaceU3S& space_u3s,
      int num_observables, 
      const std::vector<std::vector<u3shell::IndexedOperatorLabelsU3S>>& observable_symmetries_u3s,
      std::vector<std::vector<spncci::SectorLabelsU3S>>& observables_sectors_u3s
    )
  {
    // vector of sectors for each observable
    observables_sectors_u3s.resize(num_observables);

    // for each observable, enumerate sectors 
    for(int observable_index=0; observable_index<num_observables; ++observable_index) 
      {
        std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];
        spncci::GetSectorsU3S(space_u3s,observable_symmetries_u3s[observable_index],sectors_u3s);
      }
  }

  void InitializeU3SBlocks(
      const spncci::SpaceU3S& space_u3s,
      int num_observables, 
      const std::vector<double>& hw_values,
      const std::vector<std::vector<spncci::SectorLabelsU3S>>& observables_sectors_u3s,
      std::vector<std::vector<spncci::OperatorBlocks>>& observables_blocks_u3s
    )
  {
    // vector of blocks for u3 sectors for each hbar omega,for each observable
    observables_blocks_u3s.resize(hw_values.size());

    // For each hbar omega, zero initialize block for each observable
    // based on basis::SetOperatorToZero in operator.h
    // int total_entries = 0;
    for(int hw_index=0; hw_index<hw_values.size(); ++hw_index)
      {
        std::vector<spncci::OperatorBlocks>& observables_blocks=observables_blocks_u3s[hw_index];
        observables_blocks.resize(num_observables);

        for(int observable_index=0; observable_index<num_observables; ++observable_index)
          {
            // Note: this would all be a one-liner call to basis::SetOperatorToZero,
            // if sectors_u3s were a proper Sectors object
            spncci::OperatorBlocks& blocks=observables_blocks[observable_index];
            const std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];
            blocks.resize(sectors_u3s.size());
            for(int sector_index=0; sector_index<sectors_u3s.size(); ++sector_index)
              {
                int rows=space_u3s.GetSubspace(sectors_u3s[sector_index].bra_index()).full_dimension();
                int cols=space_u3s.GetSubspace(sectors_u3s[sector_index].ket_index()).full_dimension();
                blocks[sector_index]=spncci::OperatorBlock::Zero(rows,cols);
              }

            // // tally allocated matrix elements
            // int entries = basis::AllocatedEntries(blocks);
            // total_entries += entries;
            // std::cout
            //   << fmt::format(
            //       "  hw_index {:2d} observable_index {:2d} sectors {:4d} entries {:d} = {:e}",
            //       hw_index,observable_index,sectors_u3s.size(),entries,double(entries)
            //     )
            //   << std::endl;
          }

      } 
  }

  void PrintU3SSector(
      const std::vector<double>& hw_values,
      const std::vector<std::vector<spncci::SectorLabelsU3S>>& observables_sectors_u3s,
      std::vector<std::vector<spncci::OperatorBlocks>>& observables_blocks_u3s, //can't be constant because of chop function
      const spncci::SpaceU3S& space_u3s,
      int num_observables
    )
  // Prints out U3SSectors and blocks 
  {
    for(int observable_index=0; observable_index<num_observables; ++observable_index)
      for(int hw_index=0; hw_index<hw_values.size(); ++hw_index)
        {
          std::cout<<"observable "<<observable_index<<" hw "<<hw_values[hw_index]<<std::endl;
          const std::vector<spncci::SectorLabelsU3S>& sectors_u3s=observables_sectors_u3s[observable_index];
          spncci::OperatorBlocks& blocks_u3s=observables_blocks_u3s[hw_index][observable_index];
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
  ContractBabySpNCCIU3S(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      basis::OperatorHyperblocks<double>& observable_hyperblocks
    )
  {
    // TODO: add option to only keep lower triangle sectors if hermitian operator
    // iterate over interaction get unit tensor family index,kappa0,L0
    // iterate over target hypersectors to identify corresponding target
    // identify baby_spncci_hypersector
    // Contract interaction with unit tensor blocks and accumulate in
    //  observable_hyperblocks

    // zero initialize observable hyperblocks 
    basis::SetHyperoperatorToZero(observable_hypersectors,observable_hyperblocks);

    for(auto it=relative_observable.begin(); it!=relative_observable.end(); ++it)
      {
        // Relative rmes by unit tensor state index 
        const std::vector<double>& relative_rmes=it->second;
        
        int unit_tensor_subspace_index,kappa0,L0;
        std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;

        // Look up unit tensor subspace to get N0, x0 and S0 
        const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
            =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        int N0=unit_tensor_subspace.N0();
        const u3::SU3& x0=unit_tensor_subspace.x0();
        HalfInt S0=unit_tensor_subspace.S0();
        int etap=unit_tensor_subspace.etap();
        int eta=unit_tensor_subspace.eta();

        // Look up conjugate unit tensor subspace index
        u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
        int unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);

        //Get conjugate unit tensor subspace
        const auto& unit_tensor_subspace_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj);

        // Get conjugate relative observables 
        std::tuple<int,int,int> rme_labels_conj(unit_tensor_subspace_index_conj,kappa0,L0);
        const std::vector<double>& relative_rmes_conj=relative_observable.at(rme_labels_conj);

        //Look up observable subspace index
        int ob_subspace_index
            =observable_space.LookUpSubspaceIndex(u3shell::ObservableSubspaceLabels(N0,x0,S0,kappa0,L0));
        
        // Look up conjugate observable subspace index
        int ob_subspace_index_conj
          =observable_space.LookUpSubspaceIndex(u3shell::ObservableSubspaceLabels(-N0,u3::Conjugate(x0),S0,kappa0,L0));
        
        // Loop over unit tensor baby_spncci hypersectors 
        for(int baby_spncci_hypersector_index=0; baby_spncci_hypersector_index<baby_spncci_hypersectors.size(); ++baby_spncci_hypersector_index)
          {
            const auto& baby_spncci_hypersector=baby_spncci_hypersectors.GetHypersector(baby_spncci_hypersector_index);
            
            // Check if hypersector contains desired unit tensor subspace
            // Note: only Nnp>=Nn hypersectors exist, so will need to also contract over conjugate
            if(baby_spncci_hypersector.operator_subspace_index()!=unit_tensor_subspace_index)
              continue;

            // get bra and ket indices
            int bra_index=baby_spncci_hypersector.bra_subspace_index();
            int ket_index=baby_spncci_hypersector.ket_subspace_index();
            int rho0=baby_spncci_hypersector.multiplicity_index();

            // Look up target observable hypersector and conjugate hypersector
            int observable_hypersector_index
                =observable_hypersectors.LookUpHypersectorIndex(bra_index,ket_index,ob_subspace_index,rho0);

            int observable_hypersector_index_conj
                =observable_hypersectors.LookUpHypersectorIndex(ket_index,bra_index,ob_subspace_index_conj,rho0);

            // Get observable hyperblock and conjugate block
            basis::OperatorBlock<double>& observable_block=observable_hyperblocks[observable_hypersector_index][1];
            basis::OperatorBlock<double>& observable_block_conj
                =observable_hyperblocks[observable_hypersector_index_conj][1];

            // Get unit tensor hyperblocks
            const basis::OperatorBlocks<double>& unit_tensor_blocks=unit_tensor_hyperblocks[baby_spncci_hypersector_index];

            // conjugation factor base
            // 
            // get bra and ket subspace labels
            const auto& bra_subspace=baby_spncci_space.GetSubspace(bra_index);
            const auto& ket_subspace=baby_spncci_space.GetSubspace(ket_index);
            const u3::U3& omegap=bra_subspace.omega();
            const HalfInt& Sp=bra_subspace.S();
            const u3::U3& omega=ket_subspace.omega();
            const HalfInt& S=ket_subspace.S();

            double conjugation_factor_base
                =ParitySign(u3::ConjugationGrade(omegap)+Sp-u3::ConjugationGrade(omega)-S)
                  *sqrt(
                      1.*u3::dim(omegap)*am::dim(Sp)*u3::dim(u3::SU3(eta,0))
                      /u3::dim(omega)/am::dim(S)/u3::dim(u3::SU3(etap,0))
                    );

            // Loop over unit tensors and contract
            for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
              {
                observable_block+=relative_rmes[unit_tensor_index]*unit_tensor_blocks[unit_tensor_index];

                // Get state labels 
                int T0, S,T,Sp,Tp;
                std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(unit_tensor_index);

                // Get conjugate state index
                std::tuple<int,int,int,int,int> conjugate_state(T0,S,T,Sp,Tp);
                int unit_tensor_index_conj=unit_tensor_subspace_conj.LookUpStateIndex(conjugate_state);
                
                // calculate conjugation factor
                double conjugation_factor
                    =sqrt(am::dim(S)*am::dim(T)/am::dim(Sp)/am::dim(Tp))*conjugation_factor_base;

                observable_block_conj
                    +=conjugation_factor
                      *relative_rmes_conj[unit_tensor_index_conj]
                      *unit_tensor_blocks[unit_tensor_index_conj].transpose();
              }

          }
      }          
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
      spncci::OperatorBlocks& target_blocks_u3s
    )
  {
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


        // std::cout<<fmt::format("unit tensor subspace {} {} {} {} ",x0.Str(),S0,etap,eta)<<std::endl;
        // for(auto rme :relative_rmes)
        //   std::cout<<"  "<<rme<<std::endl;

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

                  // Recurrence computes Nnp>=Nn sectors for lgi_bra=>lgi_ket.
                  // The remaining unit tensor blocks are obtained by conjugation, i.e., those satisfying
                  //    Nnp-Nn<0
                  // or
                  //    Nnp-Nn=0, lgi_bra<lgi_ket [NOT at present]
                  // 
                  int Nnp=baby_spncci_subspace_bra.Nn();
                  int Nn=baby_spncci_subspace_ket.Nn();
                  int unit_tensor_subspace_index_conj=-1;
                  
                  // Check first condition 
                  bool conjugate_hypersector=(Nnp<Nn);

                  // Check second condition
                  // if(
                  //     ((Nnp-Nn)==0)
                  //     &&
                  //     (baby_spncci_subspace_bra.irrep_family_index()<baby_spncci_subspace_ket.irrep_family_index())
                  //   )
                  //   conjugate_hypersector=true;

                  // If we need to get conjugate hyper sector, then need conjugate labels 
                  if(conjugate_hypersector)
                    {
                      u3shell::UnitTensorSubspaceLabels unit_tensor_labels_conj(u3::Conjugate(x0),S0,eta,etap);
                      unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels_conj);
                    }

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
  }// end function




}  // namespace
