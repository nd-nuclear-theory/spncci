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

    // // old aem code:
    //
    // int sector_index=0;
    // int state_index=0;
    // // iterate over U(3)xSU(2) irreps
    // for(int u3s_subspace_index=0; u3s_subspace_index<u3s_space.size(); ++u3s_subspace_index)
    //   {
    //     SubspaceU3S u3s_subspace=u3s_space.GetSubspace(u3s_subspace_index);
    //     u3::U3 omega(u3s_subspace.omega());
    //     int kappa_max=u3::BranchingMultiplicitySO3(omega.SU3(),L);
    //     // if space contains S and omega can branch to L
    //     if(kappa_max>0 && u3s_subspace.S()==S)
    //       {
    //         //Construct subspace
    //         int dim=u3s_subspace.full_dimension();
    //         PushStateLabels(StateLabelsType(u3s_subspace_index));
    //         sector_index_lookup_[state_index]=sector_index;
    //         // increment index 
    //         ++state_index;
    //         sector_index+=kappa_max*dim;
    //       }
    //   }

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

#if 0
  // (mac): aem implementation, with ad hoc sorting order and risk of empty LS spaces
  SpaceLS::SpaceLS(const SpaceU3S& u3s_space, HalfInt J)
  {
 
    // std::cout << fmt::format("[Building LS space: J={}]",J.Str()) << std::endl;
 
    // iterate over U(3)xSU(2) subspaces
    for(int u3s_subspace_index=0; u3s_subspace_index<u3s_space.size(); ++u3s_subspace_index)
      // for(auto u3s_subspace : u3s_space)
      {
        const SubspaceU3S& u3s_subspace=u3s_space.GetSubspace(u3s_subspace_index);
        HalfInt S = u3s_subspace.S();
        // iterate through omega space
        u3::U3 omega = u3s_subspace.omega();
 
        // std::cout << fmt::format("U3S subspace {}",u3s_subspace.Str()) << std::endl;
 
        // interate over possible L values
        //
        // CAUTION (mac): I believe we risk creating empty LS spaces.
        // In the following code, the (L,S) pairs used in creating LS
        // subspaces are determined purely by triangularity JxS->L
        // without regard to whether or not this L exists in the
        // branching of any U3 subspace obtained for that S.
      
        for(int L=int(abs(S-J)); L<=int(S+J); ++L)
          {
            //std::cout << fmt::format("Trying (L,S) = ({},{})",L,S) << std::endl;
            if(ContainsSubspace(LSLabels(L,S)))
              continue;
            SubspaceLS ls_subspace(LSLabels(L,S),u3s_space);
            PushSubspace(ls_subspace);
          }
      }
    // // get dimension and starting index of last subspace
    // const SubspaceType& subspace=GetSubspace(subspaces_.size()-1);
    // HalfInt S=subspace.S();
    // u3::U3 omega=std::get<0>(subspace.GetStateLabels(subspace.size()-1));
    // int index=std::get<2>(subspace.GetStateLabels(subspace.size()-1));
    // u3::U3S omegaS(omega,S);
    // int dim=u3s_space.LookUpSubspace(omegaS).full_dimension();
    // dimension_=dim+index;
  }
#endif

#if 1  
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
#endif


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
        const std::vector<spncci::SectorLabelsU3S>& source_sector_labels,
        const basis::MatrixVector& source_sectors,
        const spncci::SpaceLS& target_space_bra,
        const spncci::SpaceLS& target_space_ket,
        const std::vector<spncci::SectorLabelsLS>& target_sector_labels,
        basis::MatrixVector& target_sectors
    )
  {
    // For a given Jp,J0,J sector
    // spncci::SpaceLS ls_space(u3s_space, J);
    //for each target, get sources and multiply by appropriate coefficient.
    target_sectors.resize(target_sector_labels.size());
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
        Eigen::MatrixXd& target_sector=target_sectors[t];
        target_sector=Eigen::MatrixXd::Zero(target_dim_bra,target_dim_ket);
       
        // Extract target labels 
        int L0(sector_labels.L0());
        HalfInt S0(sector_labels.S0());
        
        int L, Lp;
        HalfInt S,Sp;

        std::tie(Lp,Sp)=bra_subspace.labels();
        std::tie(L,S)=ket_subspace.labels();

        double Jcoef=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp);
        // std::cout<<fmt::format("{} {} {}  {} {} {}  {} {} {}    {}",L,S,J,L0,S0,J0,Lp,Sp,Jp,Jcoef)<<std::endl;
        
        // std::cout<<"starting loop over sources "<<std::endl;
        for(int s=0; s<source_sector_labels.size(); ++s)
          {
            const spncci::SectorLabelsU3S& source_labels=source_sector_labels[s];
            const Eigen::MatrixXd& source_sector=source_sectors.at(s);

            if(L0!=source_labels.L0())
              continue;
            
            int source_index_ket=source_labels.ket_index();
            int source_index_bra=source_labels.bra_index();
            
            // Check if sector is source for target sector
            if(not bra_subspace.ContainsState(std::tuple<int>(source_index_bra)))
              continue;
            if(not ket_subspace.ContainsState(std::tuple<int>(source_index_ket)))
              continue;

            // (indexp,index)->position of upper left corner of subsector
            int indexp=bra_subspace.sector_index(bra_subspace.LookUpStateIndex(std::tuple<int>(source_index_bra)));
            int index=ket_subspace.sector_index(ket_subspace.LookUpStateIndex(std::tuple<int>(source_index_ket)));

            //Extract source operator labels 
            const u3::SU3& x0=source_labels.x0();
            const HalfInt& S0=source_labels.S0();
            int kappa0=source_labels.kappa0();
            int rho0=source_labels.rho0();

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
                  double Wcoef=u3::WCached(w_cache,x,kappa,L,x0,kappa0,L0,xp,kappa_p,Lp,rho0);
                  // std::cout<<x.Str()<<"  "<<kappa<<"  "<<L<<"  "<<x0.Str()<<"  "<<kappa0
                            // <<"  "<<L0<<"  "<<xp.Str()<<"  "<<kappa_p<<"  "<<Lp<<"  "<<rho0<<std::endl;
                  int start_indexp=(kappa_p-1)*source_dimp+indexp;
                  int start_index=(kappa-1)*source_dim+index;
                  // std::cout<<"branching"<<std::endl
                  //   <<"W "<<Wcoef<<"  Jcoef "<<Jcoef<<std::endl
                  //   <<source_sector<<std::endl<<std::endl;
                  // std::cout<<" target "<<start_indexp<<"  "<< start_index<<"  "<< source_dimp<<"  "<< source_dim<<std::endl
                  // <<target_sector.rows()<<"  "<<target_sector.cols()<<std::endl<<std::endl;
                  target_sector.block(start_indexp,start_index,source_dimp,source_dim)+=Jcoef*Wcoef*source_sector;
                }
          }
      }
  }

  void
  ConstructOperatorMatrix(
    const spncci::SpaceLS& bra_source_space,
    const spncci::SpaceLS& ket_source_space,
    std::vector<spncci::SectorLabelsLS>& source_sector_labels,
    basis::MatrixVector& source_sectors,
    Eigen::MatrixXd& operator_matrix
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
    
    operator_matrix=Eigen::MatrixXd::Zero(bra_matrix_dim,ket_matrix_dim);
    
    // For each sector in source sectors, get bra and ket indices, look-up matrix index
    // and accumulate in full matrix 
    for(int s=0; s<source_sectors.size(); ++s)
      {
        const spncci::SectorLabelsLS& sector_labels=source_sector_labels[s];
        int subspace_index_bra=sector_labels.bra_index();
        int subspace_index_ket=sector_labels.ket_index();

        int matrix_index_bra=bra_matrix_index_lookup[subspace_index_bra];
        int matrix_index_ket=ket_matrix_index_lookup[subspace_index_ket];

        Eigen::MatrixXd& sector=source_sectors[s];
        operator_matrix.block(matrix_index_bra,matrix_index_ket,sector.rows(),sector.cols())
          +=sector;
      }
  }


}  // namespace
