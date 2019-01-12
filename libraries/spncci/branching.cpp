/****************************************************************
  branching.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching.h"

#include <set>
#include <sstream>
#include <iostream>

#include "am/halfint.h"
#include "am/halfint_fmt.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"
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
        state_irrep_family_index_.push_back(baby_spncci_subspace.irrep_family_index());
        state_baby_spncci_subspace_index_.push_back(baby_spncci_subspace_index);
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
            "  index {} omegaS {} sigmaSPN {} multiplicity {} offset {}",
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

#if 0
  // FOR FUTURE USE: constructor with sorted (omega,S) labels

  SpaceSpU3S::SpaceSpU3S(const BabySpNCCISpace& baby_spncci_space)
  {

    // collect (omega,S) labels
    std::set<u3::U3S> omegaS_set;
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // accumulate (omega,S) label
        u3::U3S omegaS = u3::U3S(baby_spncci_subspace.omega(),baby_spncci_subspace.S());
        omegaS_set.insert(omegaS);
      }

    // create subspaces
    for (const u3::U3S& omegaS : omegaS_set)
        PushSubspace(SubspaceSpU3S(omegaS,baby_spncci_space));
  }
#endif

  std::string SpaceSpU3S::DebugStr(bool show_subspaces) const
  {
    std::ostringstream os;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // set up alias
        const SubspaceType& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "subspace_index {} labels {} size {} full_dimension {}",
            subspace_index,subspace.LabelStr(),subspace.size(),subspace.full_dimension()
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
              state_irrep_family_index_.push_back(spu3s_state.irrep_family_index());
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
            "  index {} omegaS {} kappa {} sigmaSPN {} multiplicity {} offset {}",
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
        const SubspaceType& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "subspace_index {} labels {} size {} full_dimension {}",
            subspace_index,subspace.LabelStr(),subspace.size(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }

  std::string SectorLabelsSpLS::Str() const
  {
    return fmt::format("( {} {}  {} {}",
      bra_index(),ket_index(), L0(),S0());
  }

  void GenerateOperatorLabelsSpLS(
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


  void GetSectorsSpLS(
    const spncci::SpaceSpLS& space_bra,
    const spncci::SpaceSpLS& space_ket,
    const std::vector<OperatorLabelsLS>& tensor_labels,
    std::vector<spncci::SectorLabelsSpLS>& sector_labels
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
  ContractAndRegroupSpU3S(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceSpU3S& target_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const std::vector<spncci::SectorLabelsSpU3S>& target_sectors_u3s,
      basis::OperatorBlocks<double>& target_blocks_u3s
    )
  {
    for(int target_sector_index=0; target_sector_index<target_sectors_u3s.size(); ++target_sector_index)
      {
        int target_bra_index,target_ket_index,target_kappa0,target_L0,target_rho0;
        u3shell::OperatorLabelsU3S target_tensor_labels;
        std::tie(target_bra_index,target_ket_index,target_tensor_labels,target_kappa0,target_L0,target_rho0)
          =target_sectors_u3s[target_sector_index].Key();

        // Get corresponding SpU3S subspaces
        const spncci::SubspaceSpU3S& bra_target_subspace=target_space.GetSubspace(target_bra_index);
        const spncci::SubspaceSpU3S& ket_target_subspace=target_space.GetSubspace(target_ket_index);

        u3::U3 omegap=bra_target_subspace.omega();
        u3::U3 omega=ket_target_subspace.omega();

        // std::cout<<"begin observables for target "<<target_sector_index<<std::endl;
        for(auto it=relative_observable.begin(); it!=relative_observable.end(); ++it)
          {
            int unit_tensor_subspace_index,kappa0,L0;
            std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;

            if((target_L0!=L0)||(target_kappa0!=kappa0))
              continue;

            // Get rmes
            const std::vector<double>& relative_rmes=it->second;

            // Get Unit Tensor subspace
            const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace
              =unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

            u3::SU3 x0;
            int etap,eta;
            HalfInt S0;
            std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();

            // Check if operator tensor labels match target
            if( (not (target_tensor_labels.x0()==x0))
              || (target_tensor_labels.S0()!=S0)
              || (target_tensor_labels.N0()!=(etap-eta))
              )
              continue;

            for(int bra_state_index=0; bra_state_index<bra_target_subspace.size(); ++bra_state_index)
              for(int ket_state_index=0; ket_state_index<ket_target_subspace.size(); ++ket_state_index)
                {
                  // if(target_sector_index==316)
                  //   std::cout<<"made it here 1"<<std::endl;
                  // extracting baby spncci state information
                  spncci::StateSpU3S baby_spncci_bra(bra_target_subspace, bra_state_index);
                  spncci::StateSpU3S baby_spncci_ket(ket_target_subspace, ket_state_index);

                  // Index of baby spncci subspace in baby_spncci_space for looking up hypersector
                  int baby_spncci_index_bra=baby_spncci_bra.baby_spncci_subspace_index();
                  int baby_spncci_index_ket=baby_spncci_ket.baby_spncci_subspace_index();

                  // Shouldn't need to get baby spncci subspace, all relevant information should be stored in spu3s states
                  // TODO: debug irrep_family_index() and Nn() for StateSpU3S and implement here

                  const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra
                          =baby_spncci_space.GetSubspace(baby_spncci_index_bra);
                  const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket
                          =baby_spncci_space.GetSubspace(baby_spncci_index_ket);


                  // (dimp,dim)->size of baby spncci hyperblock
                  int dimp=baby_spncci_bra.degeneracy();
                  int dim=baby_spncci_ket.degeneracy();

                  // offset of baby spncci brat and ket in courser grain u3s block
                  int block_index_u3s_bra=baby_spncci_bra.offset();
                  int block_index_u3s_ket=baby_spncci_ket.offset();


                  /////////////////////////////////////////////////////////////////////////////////////////////////////
                  // Checking if need conjugate hypersector
                  // Recurrence computes Nnp<=Nn sectors for lgi_bra<=lgi_ket and Nnp>Nn sectors for lgi_bra>lgi_ket.
                  // The remaining unit tensor blocks are obtained by conjugation, i.e., those satisfying
                  //    If Nnp-Nn>0
                  //    If Nnp-Nn=0, lgi_bra>lgi_ket
                  //
                  int Nnp=baby_spncci_subspace_bra.Nn();
                  int Nn=baby_spncci_subspace_ket.Nn();
                  int unit_tensor_subspace_index_conj=-1;


                  bool conjugate_hypersector=(Nnp-Nn)>0;

                  if( ((Nnp-Nn)==0)
                      &&(baby_spncci_subspace_bra.irrep_family_index()>baby_spncci_subspace_ket.irrep_family_index())
                    )
                    conjugate_hypersector=true;
                  /////////////////////////////////////////////////////////////////////////////////////////////////////

                  // If conjugate hypersector, get unit tensor subspace index of conjugate of unit tensor
                  if(conjugate_hypersector)
                    {
                      u3shell::UnitTensorSubspaceLabels unit_tensor_labels_conj(u3::Conjugate(x0),S0,eta,etap);
                      unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_labels_conj);
                    }

                  // Look up baby spncci hypersector index accomodating need to look up conjugate hypersector
                  int baby_spncci_hypersector_index
                        =conjugate_hypersector?
                        baby_spncci_hypersectors.LookUpHypersectorIndex(
                            baby_spncci_index_ket,baby_spncci_index_bra,
                            unit_tensor_subspace_index_conj,target_rho0
                          ):
                        baby_spncci_hypersectors.LookUpHypersectorIndex(
                            baby_spncci_index_bra,baby_spncci_index_ket,
                            unit_tensor_subspace_index,target_rho0
                          );

                  // If baby_spncci_hypersector doesn't exist, continue
                  if(baby_spncci_hypersector_index==-1)
                    continue;

                  // Get unit tensor hyperblocks
                  const basis::OperatorBlocks<double>& unit_tensor_blocks
                      =unit_tensor_hyperblocks[baby_spncci_hypersector_index];

                  // Loop through unit tensors in unit tensor subspace, and accumulate product with relative rmes
                  // in target sector.
                  for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
                  {
                    // If conjugate hypersector, we need to get conjugation factors and phases
                    if(conjugate_hypersector)
                      {
                        int T0, S,T,Sp,Tp;
                        std::tie(T0,Sp,Tp,S,T)=unit_tensor_subspace.GetStateLabels(unit_tensor_index);
                        std::tuple<int,int,int,int,int> conjugate_state(T0,S,T,Sp,Tp);
                        int unit_tensor_index_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj).LookUpStateIndex(conjugate_state);

                        u3::U3 omegap(baby_spncci_bra.omega());
                        u3::U3 omega(baby_spncci_ket.omega());

                        int conjugation_grade=ParitySign(
                          u3::ConjugationGrade(omega)
                          +baby_spncci_ket.S()
                          +u3::ConjugationGrade(omegap)
                          +baby_spncci_bra.S()
                          );

                        double conjugation_factor
                                =sqrt(
                                  1.*u3::dim(u3::SU3(etap,0))*u3::dim(omega)
                                  *am::dim(baby_spncci_ket.S())
                                  *am::dim(Sp)*am::dim(Tp)
                                  /u3::dim(u3::SU3(eta,0))/u3::dim(omegap)
                                  /am::dim(baby_spncci_bra.S())
                                  /am::dim(S)/am::dim(T)
                                  );

                        // Accumulate sum over conjugated hypersector and realtive unit tensor
                        target_blocks_u3s[target_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim)
                          +=conjugation_grade*conjugation_factor
                            *relative_rmes[unit_tensor_index]
                            *unit_tensor_blocks[unit_tensor_index_conj].transpose();
                      }

                    else
                      target_blocks_u3s[target_sector_index].block(block_index_u3s_bra,block_index_u3s_ket,dimp,dim)
                        +=relative_rmes[unit_tensor_index]*unit_tensor_blocks[unit_tensor_index];
                   }
                }
          }

      }
  }


  // void
  // ContractAndRegroupLSJ(
  //       const HalfInt& Jp,const HalfInt& J0, const HalfInt& J,
  //       u3::WCoefCache& w_cache,
  //       const spncci::SpaceSpU3S& u3s_space,
  //       const std::vector<spncci::SectorLabelsSpU3S>& source_sector_labels,
  //       const basis::MatrixVector& source_sectors,
  //       const spncci::SpaceSpLS& target_space_bra,
  //       const spncci::SpaceSpLS& target_space_ket,
  //       const std::vector<spncci::SectorLabelsSpLS>& target_sector_labels,
  //       basis::MatrixVector& target_blocks
  //   )
  // {
  //   // For a given Jp,J0,J sector
  //   // spncci::SpaceLS ls_space(u3s_space, J);
  //   //for each target, get sources and multiply by appropriate coefficient.
  //   target_blocks.resize(target_sector_labels.size());
  //   for(int t=0; t<target_sector_labels.size(); ++t)
  //     {
  //       // std::cout<<"t "<<t<<std::endl;

  //       const spncci::SectorLabelsSpLS& sector_labels=target_sector_labels[t];
  //       const spncci::SubspaceSpLS&
  //         ket_subspace=target_space_ket.GetSubspace(sector_labels.ket_index());
  //       const spncci::SubspaceSpLS&
  //         bra_subspace=target_space_bra.GetSubspace(sector_labels.bra_index());

  //       int target_dim_bra=bra_subspace.full_dimension();
  //       int target_dim_ket=ket_subspace.full_dimension();

  //       // Zero initialize
  //       Eigen::MatrixXd& target_sector=target_blocks[t];
  //       target_sector=Eigen::MatrixXd::Zero(target_dim_bra,target_dim_ket);

  //       // Extract target labels
  //       int L0(sector_labels.L0());
  //       HalfInt S0(sector_labels.S0());

  //       int L, Lp;
  //       HalfInt S,Sp;

  //       std::tie(Lp,Sp)=bra_subspace.labels();
  //       std::tie(L,S)=ket_subspace.labels();

  //       double Jcoef=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp);
  //       // std::cout<<fmt::format("{} {} {}  {} {} {}  {} {} {}    {}",L,S,J,L0,S0,J0,Lp,Sp,Jp,Jcoef)<<std::endl;
  //       // States are actually baby spncci subspaces
  //       int bra_offset=0; //offset within u3ssubspace
  //       for(int bra_state_index=0; bra_state_index<bra_subspace.size(); ++bra_state_index)
  //         {
  //           int ket_offset=0; //offset within u3ssubspace
  //           for(int ket_state_index=0; ket_state_index<ket_subspace.size(); ++ket_state_index)
  //             {// std::cout<<"starting loop over sources "<<std::endl;
  //               for(int s=0; s<source_sector_labels.size(); ++s)
  //                 {
  //                   const spncci::SectorLabelsSpU3S& source_labels=source_sector_labels[s];
  //                   const Eigen::MatrixXd& source_sector=source_sectors.at(s);

  //                   if(L0!=source_labels.L0())
  //                     continue;

  //                   if(S0!=source_labels.S0())
  //                     continue;

  //                   spncci::StateSpLS bra_state(bra_subspace,bra_state_index);
  //                   spncci::StateSpLS ket_state(ket_subspace,ket_state_index);

  //                   if(bra_state.spu3s_subspace_index()!=source_labels.bra_index())
  //                     continue;

  //                   if(ket_state.spu3s_subspace_index()!=source_labels.ket_index())
  //                     continue;

  //                   int spu3s_subspace_index_bra=bra_state.spu3s_subspace_index();
  //                   int spu3s_subspace_index_ket=ket_state.spu3s_subspace_index();


  //                   int source_index_ket=source_labels.ket_index();
  //                   int source_index_bra=source_labels.bra_index();

  //                   // (indexp,index)->position of upper left corner of subsector (full ls subspace)
  //                   int indexp=bra_state.offset();
  //                   int index=ket_state.offset();

  //                   // size of state degeneracy (upsilon and gamma)
  //                   int dimp=bra_state.degeneracy();
  //                   int dim=ket_state.degeneracy();

  //                   //Extract source operator labels
  //                   const u3::SU3& x0=source_labels.x0();
  //                   const HalfInt& S0=source_labels.S0();
  //                   int kappa0=source_labels.kappa0();
  //                   int rho0=source_labels.rho0();

  //                   const spncci::SubspaceSpU3S&
  //                     u3s_subspace_bra=u3s_space.GetSubspace(spu3s_subspace_index_bra);
  //                   const spncci::SubspaceSpU3S&
  //                     u3s_subspace_ket=u3s_space.GetSubspace(spu3s_subspace_index_ket);

  //                   // Look up index of corresponding baby spncci state in u3s_subspace so we can look of offset in u3s_unit tensor block
  //                   int baby_spnci_index_u3s_subspace_bra=u3s_subspace_bra.LookUpStateIndex(bra_state.sigmaSPN());
  //                   int baby_spnci_index_u3s_subspace_ket=u3s_subspace_ket.LookUpStateIndex(ket_state.sigmaSPN());

  //                   // Offset in u3s_block


  //                   //source sector dimensions
  //                   int source_dimp=u3s_subspace_bra.full_dimension();
  //                   int source_dim=u3s_subspace_ket.full_dimension();

  //                   // Extract source state labels
  //                   const u3::U3S& omegaSp=u3s_subspace_bra.labels();
  //                   const u3::U3S& omegaS=u3s_subspace_ket.labels();
  //                   assert(omegaSp.S()==Sp);
  //                   assert(omegaS.S()==S);
  //                   u3::SU3 xp(omegaSp.U3().SU3());
  //                   u3::SU3 x(omegaS.U3().SU3());

  //                   // Get branching multiplicities
  //                   int kappa_p=bra_state.kappa();
  //                   int kappa=ket_state.kappa();

  //                   int kappa_max_p=u3::BranchingMultiplicitySO3(xp,Lp);
  //                   int kappa_max=u3::BranchingMultiplicitySO3(x,L);

  //                   // Generate coefficient for each kappa and kappa_p and accumulate in
  //                   // target sector. Source sector dimensions are source_dimp x source_dim
  //                   // starting position given by :
  //                   //    ((kappa_p-1)*source_dimp+indexp, (kappa-1)*source_dim+index)
  //                   double Wcoef=u3::WCached(w_cache,x,kappa,L,x0,kappa0,L0,xp,kappa_p,Lp,rho0);
  //                   // std::cout<<x.Str()<<"  "<<kappa<<"  "<<L<<"  "<<x0.Str()<<"  "<<kappa0
  //                             // <<"  "<<L0<<"  "<<xp.Str()<<"  "<<kappa_p<<"  "<<Lp<<"  "<<rho0<<std::endl;


  //                   // int start_indexp=(kappa_p-1)*source_dimp+indexp;
  //                   // int start_index=(kappa-1)*source_dim+index;

  //                   int start_indexp=bra_state.offset();
  //                   int start_index=ket_state.offset();

  //                   // std::cout<<"branching"<<std::endl
  //                   //   <<"W "<<Wcoef<<"  Jcoef "<<Jcoef<<std::endl
  //                   //   <<source_sector<<std::endl<<std::endl;
  //                   std::cout<<" target "<<t<<"  "<<start_indexp<<"  "<< start_index<<"  "<< source_dimp<<"  "<< source_dim<<std::endl
  //                   <<target_sector.rows()<<"  "<<target_sector.cols()<<"  "<<source_sector.rows()<<"  "<<source_sector.cols()<<std::endl<<std::endl;
  //                   target_sector.block(start_indexp,start_index,source_dimp,source_dim)+=Jcoef*Wcoef*source_sector;
  //                 }
  //             }
  //         }
  //     }
  // }

  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to J level
  ////////////////////////////////////////////////////////////////

  SubspaceSpJ::SubspaceSpJ(HalfInt J, const SpaceSpLS& spls_space)
  {

    // save labels
    labels_ = J;

    // scan SpLS space for states to accumulate
    //
    // If LS subspace is triangular with J, then all states in that LS
    // subspace are simply "copied" as states in this J subspace.
    for(int spls_subspace_index=0; spls_subspace_index<spls_space.size(); ++spls_subspace_index)
      {

        // set up alias
        const SubspaceSpLS& spls_subspace = spls_space.GetSubspace(spls_subspace_index);

        // short circuit if LS subspace not relevant to current J subspace
        if (!am::AllowedTriangle(spls_subspace.L(),spls_subspace.S(),J))
          continue;

        // accumulate J-branched states
        for (int spls_state_index=0; spls_state_index<spls_subspace.size(); ++spls_state_index)
          {
            // push state
            StateSpLS spls_state(spls_subspace,spls_state_index);
            StateLabelsSpJ spj_state_labels = StateLabelsSpJ(
                spls_state.LS(),
                spls_state.omega(),
                spls_state.kappa(),
                spls_state.sigmaSPN()
                );
            PushStateLabels(spj_state_labels,spls_state.degeneracy());

            // record auxiliary state information
            state_gamma_max_.push_back(spls_state.gamma_max());
            state_irrep_family_index_.push_back(spls_state.irrep_family_index());
            state_baby_spncci_subspace_index_.push_back(spls_state.baby_spncci_subspace_index());
            state_spu3s_subspace_index_.push_back(spls_state.spu3s_subspace_index());
            state_spls_subspace_index_.push_back(spls_subspace_index);
          }

      }
  }

  std::string SubspaceSpJ::LabelStr() const
  {
    return J().Str();
  }

  std::string SubspaceSpJ::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpJ state(*this,state_index);

        os << fmt::format(
            "  index {} labels {} degeneracy {} offset {}",
            state_index,state.LabelStr(),state.degeneracy(),state.offset()
          ) << std::endl;
      }

    return os.str();
  }

  std::string StateSpJ::LabelStr() const
  {
    return fmt::format("[{} {} {} {} {} {}]",L(),S().Str(),omega().Str(),kappa(),sigmaSPN().Str(),J().Str());
  }

  SpaceSpJ::SpaceSpJ(const std::vector<HalfInt>& J_values, const SpaceSpLS& spls_space)
  {
    for (HalfInt J : J_values)
      PushSubspace(SubspaceSpJ(J,spls_space));
  }

  std::string SpaceSpJ::DebugStr(bool show_subspaces) const
  {
    std::ostringstream os;

    for (int subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        // set up alias
        const SubspaceType& subspace = GetSubspace(subspace_index);

        os << fmt::format(
            "subspace_index {} labels {} size {} full_dimension {}",
            subspace_index,subspace.LabelStr(),subspace.size(),subspace.full_dimension()
          ) << std::endl;
        if (show_subspaces)
          os << subspace.DebugStr();

      }

    return os.str();
  }

  SectorsSpJ::SectorsSpJ(
        const SpaceSpJ& space,
        HalfInt J0,
        basis::SectorDirection sector_direction
    )
  {
    for (int bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
      for (int ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
	{
          // enforce canonical ordering
          if (
              (sector_direction == basis::SectorDirection::kCanonical)
              && !(bra_subspace_index<=ket_subspace_index)
            )
            continue;

          // retrieve subspaces
          const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
          const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);

          // verify angular momentum selection rule
          bool allowed = true;
          allowed &= am::AllowedTriangle(ket_subspace.J(),J0,bra_subspace.J());
          // to add parity selection: allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);

          // push sector
	  if (allowed)
            PushSector(SectorType(bra_subspace_index,ket_subspace_index,bra_subspace,ket_subspace));
        }
  }

  // void
  // ConstructOperatorMatrix(
  //   const spncci::SpaceSpLS& bra_source_space,
  //   const spncci::SpaceSpLS& ket_source_space,
  //   std::vector<spncci::SectorLabelsSpLS>& source_sector_labels,
  //   basis::MatrixVector& source_sectors,
  //   Eigen::MatrixXd& operator_matrix
  //   )
  // {
  //   // Get size of matrix
  //   // Generate look up table for sector indices
  //   int bra_index=0;
  //   std::map<int,int> bra_matrix_index_lookup;
  //   for(int s=0; s<bra_source_space.size(); ++s)
  //     {
  //       const spncci::SubspaceSpLS& subspace=bra_source_space.GetSubspace(s);
  //       int full_dimension=subspace.full_dimension();
  //       bra_matrix_index_lookup[s]=bra_index;
  //       bra_index+=full_dimension;
  //     }

  //   int bra_matrix_dim=bra_index;

  //   int ket_index=0;
  //   std::map<int,int> ket_matrix_index_lookup;
  //   for(int s=0; s<ket_source_space.size(); ++s)
  //     {
  //       const spncci::SubspaceSpLS& subspace=ket_source_space.GetSubspace(s);
  //       int full_dimension=subspace.full_dimension();
  //       ket_matrix_index_lookup[s]=ket_index;
  //       ket_index+=full_dimension;
  //     }

  //   int ket_matrix_dim=ket_index;

  //   operator_matrix=Eigen::MatrixXd::Zero(bra_matrix_dim,ket_matrix_dim);

  //   // For each sector in source sectors, get bra and ket indices, look-up matrix index
  //   // and accumulate in full matrix
  //   for(int s=0; s<source_sectors.size(); ++s)
  //     {
  //       const spncci::SectorLabelsSpLS& sector_labels=source_sector_labels[s];
  //       int subspace_index_bra=sector_labels.bra_index();
  //       int subspace_index_ket=sector_labels.ket_index();

  //       int matrix_index_bra=bra_matrix_index_lookup[subspace_index_bra];
  //       int matrix_index_ket=ket_matrix_index_lookup[subspace_index_ket];

  //       Eigen::MatrixXd& sector=source_sectors[s];
  //       operator_matrix.block(matrix_index_bra,matrix_index_ket,sector.rows(),sector.cols())
  //         +=sector;
  //     }
  // }


}  // namespace
