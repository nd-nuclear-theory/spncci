/****************************************************************
  branching_u3lsj.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching_u3lsj.h"

#include <fstream>
#include <iostream>

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "cppformat/format.h"
#include "mcutils/parsing.h"


namespace spncci
{

  SubspaceLS::SubspaceLS(const spncci::LSLabels& ls_labels, const SpaceU3S& u3s_space)
  {

    // std::cout << fmt::format("[Building LS subpsace: (L,S) = ({},{})]",L,S) << std::endl;

    // save labels
    labels_ = ls_labels;

    // scan SpaceU3S for states to accumulate
    int substate_offset = 0;  // accumulated offset
    for(int u3s_subspace_index=0; u3s_subspace_index<u3s_space.size(); ++u3s_subspace_index)
      {

        // set up alias
        const SubspaceU3S& u3s_subspace = u3s_space.GetSubspace(u3s_subspace_index);

        // determine branching multiplicity to the specified L
        int kappa_max=u3::BranchingMultiplicitySO3(u3s_subspace.x(),L());

        // short circuit if U3S subspace not relevant to current LS subspace
        if ((kappa_max==0)||(S()!=u3s_subspace.S()))
          continue;

        // push state
        PushStateLabels(StateLabelsType(u3s_subspace_index));

        // record state multiplicity indexing information
        state_substate_offset_.push_back(substate_offset);
        int state_dimension = u3s_subspace.full_dimension();
        state_dimension_.push_back(state_dimension);

        // store state U3 irrep information
        state_omega_.push_back(u3s_subspace.omega());

        // store state symplectic irrep information -- NOT WELL DEFINED
        // state_sigmaSPN_.push_back(baby_spncci_subspace.omegaSPN());
        // state_gamma_max_.push_back(baby_spncci_subspace.gamma_max());

        // accumulate offset for next state
        substate_offset += kappa_max*state_dimension;
        
      }

    // store final full dimension
    full_dimension_ = substate_offset;

    // std::cout << fmt::format("Subspace: size {}, full_dimension {}",size(),full_dimension()) << std::endl;

  }

  std::string SubspaceLS::Str() const
  {
    return fmt::format("[{} {}]",L(),S());
  }


  // (mac): clean implementation, with sorted LS spaces
  SpaceLS::SpaceLS(const SpaceU3S& u3s_space, HalfInt J)
  {

    // collect (L,S) branchings
    std::set<LSLabels> ls_labels_set;
    for(int u3s_subspace_index=0; u3s_subspace_index<u3s_space.size(); ++u3s_subspace_index)
      {

        // set up alias
        const SubspaceU3S& u3s_subspace = u3s_space.GetSubspace(u3s_subspace_index);

        // find branching L values yielding given J
        u3::SU3 x = u3s_subspace.omega().SU3();
        HalfInt S = u3s_subspace.S();
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
      {
        PushSubspace(SubspaceLS(ls_labels,u3s_space));
      }
  }



  std::string SectorLabelsLS::Str() const
  {
    return fmt::format("( {} {}  {} {}", 
      bra_index(),ket_index(), L0(),S0());
  }

  void GenerateOperatorLabelsLS(
    const HalfInt& J0, 
    std::vector<OperatorLabelsLS>& tensor_labels
    )
  {
    int L0_min=std::max(int(J0-2),0);
    for(int L0=L0_min; L0<=int(J0+2); L0++)
      for(int S0=0; S0<=2; S0++)
        {
        if(am::AllowedTriangle(L0,S0,J0))
          tensor_labels.emplace_back(L0,S0);    
        }  
  }

  void GetSectorsLS(
    const spncci::SpaceLS& space_bra, 
    const spncci::SpaceLS& space_ket,
    const std::vector<OperatorLabelsLS>& tensor_labels,
    std::vector<spncci::SectorLabelsLS>& sector_labels
    )
  {
    for(auto tensor_label: tensor_labels)
      for(int j=0; j<space_ket.size(); ++j)
        for(int i=0; i<space_bra.size(); ++i)
          {
            int L,Lp,L0;
            HalfInt S, Sp, S0;
            std::tie(Lp,Sp)=space_bra.GetSubspace(i).labels();
            std::tie(L,S)=space_ket.GetSubspace(j).labels();
            std::tie(L0,S0)=tensor_label;
            if(am::AllowedTriangle(L,Lp,L0) && am::AllowedTriangle(S,Sp,S0))
              sector_labels.emplace_back(i,j,tensor_label);
          }
  }



  void 
  ContractAndRegroupLSJ(
        const HalfInt& Jp,const HalfInt& J0, const HalfInt& J,
        u3::WCoefCache& w_cache,
        const spncci::SpaceU3S& u3s_space,
        const u3shell::ObservableSpaceU3S& observable_space,
        const spncci::ObservableHypersectorsU3S& source_hypersectors,
        const spncci::OperatorBlocks& source_hyperblocks,
        const spncci::SpaceLS& target_space_bra,
        const spncci::SpaceLS& target_space_ket,
        const std::vector<spncci::SectorLabelsLS>& target_sector_labels,
        spncci::OperatorBlocks& target_blocks
    )
  {
    //NEW VERSION
    // For a given Jp,J0,J sector
    // spncci::SpaceLS ls_space(u3s_space, J);
    //for each target, get sources and multiply by appropriate coefficient.
    // std::cout<<"target sector "<<Jp<<"  "<<J0<<"  "<<J<<"  ............................................................."<<std::endl;
    target_blocks.resize(target_sector_labels.size());
    for(int t=0; t<target_sector_labels.size(); ++t)
      {
        // std::cout<<"t "<<t<<std::endl;

        const spncci::SectorLabelsLS& sector_labels=target_sector_labels[t];
        const spncci::SubspaceLS& 
          ket_subspace=target_space_ket.GetSubspace(sector_labels.ket_index());
        const spncci::SubspaceLS& 
          bra_subspace=target_space_bra.GetSubspace(sector_labels.bra_index());
        
        int target_dim_bra=bra_subspace.full_dimension();
        int target_dim_ket=ket_subspace.full_dimension();
        // Zero initialize
        spncci::OperatorBlock& target_block=target_blocks[t];
        target_block=spncci::OperatorBlock::Zero(target_dim_bra,target_dim_ket);
       


        // Extract target labels 
        int L0(sector_labels.L0());
        HalfInt S0(sector_labels.S0());
        
        int L, Lp;
        HalfInt S,Sp;

        std::tie(Lp,Sp)=bra_subspace.labels();
        std::tie(L,S)=ket_subspace.labels();

        spncci::MatrixFloatType Jcoef=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp);
        // std::cout<<fmt::format("{} {} {}  {} {} {}  {} {} {}    {}",L,S,J,L0,S0,J0,Lp,Sp,Jp,Jcoef)<<std::endl;
        
        // std::cout<<"starting loop over sources "<<std::endl;
        for(int s=0; s<source_hypersectors.size(); ++s)
          {
            const auto& source_hypersector=source_hypersectors.GetHypersector(s);

            // const spncci::SectorLabelsU3S& source_labels=source_sector_labels[s];
            const spncci::OperatorBlock& source_block=source_hyperblocks.at(s);

            int source_index_ket=source_hypersector.ket_subspace_index();
            int source_index_bra=source_hypersector.bra_subspace_index();
            int source_index_operator=source_hypersector.operator_subspace_index();
            

            // Check if sector is source for target sector
            if(not bra_subspace.ContainsState(std::tuple<int>(source_index_bra)))
              continue;

            if(not ket_subspace.ContainsState(std::tuple<int>(source_index_ket)))
              continue;

            const auto& observable_subspace=observable_space.GetSubspace(source_index_operator);

            if(L0!=observable_subspace.L0())
              continue;
            
            if(S0!=observable_subspace.S0())
              continue;

            // (indexp,index)->position of upper left corner of subsector
            int indexp=bra_subspace.sector_index(bra_subspace.LookUpStateIndex(std::tuple<int>(source_index_bra)));
            int index=ket_subspace.sector_index(ket_subspace.LookUpStateIndex(std::tuple<int>(source_index_ket)));

            //Extract source hypersector labels 
            const u3::SU3& x0=observable_subspace.x0();
            const HalfInt& S0=observable_subspace.S0();
            int kappa0=observable_subspace.kappa0();
            int rho0=source_hypersector.multiplicity_index();

            const spncci::SubspaceU3S& 
              u3s_subspace_bra=u3s_space.GetSubspace(source_index_bra);
            const spncci::SubspaceU3S& 
              u3s_subspace_ket=u3s_space.GetSubspace(source_index_ket);

            //source sector dimensions 
            int source_dimp=u3s_subspace_bra.full_dimension();
            int source_dim=u3s_subspace_ket.full_dimension();

            // Extract source state labels 
            const u3::U3S& omegaSp=u3s_subspace_bra.labels();
            const u3::U3S& omegaS=u3s_subspace_ket.labels();
            assert(omegaSp.S()==Sp);
            assert(omegaS.S()==S);
            u3::SU3 xp(omegaSp.U3().SU3());
            u3::SU3 x(omegaS.U3().SU3());

            // Get branching multiplicity
            int kappa_max_p=u3::BranchingMultiplicitySO3(xp,Lp);
            int kappa_max=u3::BranchingMultiplicitySO3(x,L);
            
            // std::cout<<"starting loop over kappap"<<std::endl;
            for(int kappa_p=1; kappa_p<=kappa_max_p; ++kappa_p)
              for(int kappa=1; kappa<=kappa_max; ++kappa)
                {
                  // Generate coefficient for each kappa and kappa_p and accumulate in
                  // target sector. Source sector dimensions are source_dimp x source_dim
                  // starting position given by :
                  //    ((kappa_p-1)*source_dimp+indexp, (kappa-1)*source_dim+index) 
                  spncci::MatrixFloatType Wcoef=u3::WCached(w_cache,x,kappa,L,x0,kappa0,L0,xp,kappa_p,Lp,rho0);
                  // std::cout<<x.Str()<<"  "<<kappa<<"  "<<L<<"  "<<x0.Str()<<"  "<<kappa0
                  //           <<"  "<<L0<<"  "<<xp.Str()<<"  "<<kappa_p<<"  "<<Lp<<"  "<<rho0<<std::endl;
                  int start_indexp=(kappa_p-1)*source_dimp+indexp;
                  int start_index=(kappa-1)*source_dim+index;

                  target_block.block(start_indexp,start_index,source_dimp,source_dim)+=Jcoef*Wcoef*source_block;
                  // std::cout<<target_block<<std::endl<<std::endl;
                }
          }
      }
  }


  void
  ConstructOperatorMatrix(
    const spncci::SpaceLS& bra_source_space,
    const spncci::SpaceLS& ket_source_space,
    const std::vector<spncci::SectorLabelsLS>& source_sector_labels,
    const spncci::OperatorBlocks& source_blocks,
    spncci::OperatorBlock& operator_matrix
    )
  {
    // Get size of matrix
    // Generate look up table for sector indices
    int bra_index=0;
    std::map<int,int> bra_matrix_index_lookup;
    for(int s=0; s<bra_source_space.size(); ++s)
      {
        const spncci::SubspaceLS& subspace=bra_source_space.GetSubspace(s);
        int full_dimension=subspace.full_dimension();
        bra_matrix_index_lookup[s]=bra_index;
        bra_index+=full_dimension;
      }
    
    int bra_matrix_dim=bra_index;

    int ket_index=0;
    std::map<int,int> ket_matrix_index_lookup;
    for(int s=0; s<ket_source_space.size(); ++s)
      {
        const spncci::SubspaceLS& subspace=ket_source_space.GetSubspace(s);
        int full_dimension=subspace.full_dimension();
        ket_matrix_index_lookup[s]=ket_index;
        ket_index+=full_dimension;
      }

    int ket_matrix_dim=ket_index;
    
    operator_matrix=spncci::OperatorBlock::Zero(bra_matrix_dim,ket_matrix_dim);
    
    // For each sector in source sectors, get bra and ket indices, look-up matrix index
    // and accumulate in full matrix 
    for(int s=0; s<source_blocks.size(); ++s)
      {
        const spncci::SectorLabelsLS& sector_labels=source_sector_labels[s];
        int subspace_index_bra=sector_labels.bra_index();
        int subspace_index_ket=sector_labels.ket_index();

        int matrix_index_bra=bra_matrix_index_lookup[subspace_index_bra];
        int matrix_index_ket=ket_matrix_index_lookup[subspace_index_ket];

        const spncci::OperatorBlock& block=source_blocks[s];
        operator_matrix.block(matrix_index_bra,matrix_index_ket,block.rows(),block.cols())
          +=block;
      }
  }

  void ConstructBranchedBlock(
      u3::WCoefCache& w_cache,
      const spncci::SpaceU3S& space_u3s,
      const u3shell::ObservableSpaceU3S& observable_space,
      const spncci::ObservableHypersectorsU3S& observable_hypersectors,
      const spncci::OperatorBlocks& observable_blocks,
      std::map<HalfInt,spncci::SpaceLS>& spaces_lsj,
      int J0,
      const typename spncci::SectorsSpJ::SectorType& sector_spj,
      spncci::OperatorBlock& block_spj
    )
  {
    // read off J values
    const HalfInt bra_J = sector_spj.bra_subspace().J();
    const HalfInt ket_J = sector_spj.ket_subspace().J();

    // determine allowed LS sectors
    const spncci::SpaceLS& bra_space_lsj = spaces_lsj[bra_J];
    const spncci::SpaceLS& ket_space_lsj = spaces_lsj[ket_J];

    // determine set of (L0,S0) labels for this observable (triangular with J0)
    std::vector<spncci::OperatorLabelsLS> operator_labels_ls;
    spncci::GenerateOperatorLabelsLS(J0,operator_labels_ls);

    // ...
    std::vector<spncci::SectorLabelsLS> sectors_lsj;
    spncci::GetSectorsLS(bra_space_lsj,ket_space_lsj,operator_labels_ls,sectors_lsj);

    // branch LS sectors to LSJ
    spncci::OperatorBlocks matrices_lsj;  
    // std::cout<<"contract and regroup"<<std::endl;
    ContractAndRegroupLSJ(
      bra_J,J0,ket_J,w_cache,space_u3s,
      observable_space,observable_hypersectors,observable_blocks,
      bra_space_lsj,ket_space_lsj,sectors_lsj,matrices_lsj
    );

    // std::cout<<"get operator matrix "<<std::endl;
    // collect LSJ sectors into J matrix
    spncci::ConstructOperatorMatrix(
        bra_space_lsj,ket_space_lsj,sectors_lsj,matrices_lsj,
        block_spj
      );

  }

}  // namespace
