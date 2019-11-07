/****************************************************************
  branching.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/branching2.h"

#include <set>
#include <sstream>
#include <fstream>
#include <iostream>


#include "am/wigner_gsl.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "spncci/io_control.h"
// #include "spncci/unit_tensor.h"




namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // SpNCCI basis branched to U3S level
  ////////////////////////////////////////////////////////////////

  SubspaceSpBasis::SubspaceSpBasis(
    const HalfInt& J,
    const u3shell::U3SPN& sigmaSPN,
    int irrep_family_index,
    const BabySpNCCISpace& baby_spncci_space
  )
  {

    // save labels
    // labels_ = sigmaSPN;
    labels_ = irrep_family_index;
    sigmaSPN_=sigmaSPN;
    dimension_=0;
    // scan BabySpNCCISpace for states to accumulate
    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace = baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        if(baby_spncci_subspace.irrep_family_index()!=irrep_family_index)
          continue;

        // push state
        const u3::U3& omega=baby_spncci_subspace.omega();
        const HalfInt& S=baby_spncci_subspace.S();
        int gamma_max=baby_spncci_subspace.gamma_max();
        int upsilon_max=baby_spncci_subspace.upsilon_max();

        //Should be the same for all states with same irrep_family_index
        gamma_max_=gamma_max;

        MultiplicityTagged<int>::vector so3_irreps=u3::BranchingSO3(omega.SU3());
        for(auto& tagged_L : so3_irreps)
          {
            int L=tagged_L.irrep;
            int kappa_max=tagged_L.tag;
            if(am::AllowedTriangle(L,S,J))
              {
                int degeneracy=gamma_max*upsilon_max*kappa_max;
                omegaLLabels state_labels(omega,L);
                PushStateLabels(state_labels,degeneracy);

                // record auxiliary state information
                state_kappa_max_.push_back(kappa_max);
                state_upsilon_max_.push_back(upsilon_max);
                state_baby_spncci_subspace_index_.push_back(baby_spncci_subspace_index);

                // increment dimension of subspace 
                dimension_+=degeneracy;
              }
          }
      }
  }

  std::string SubspaceSpBasis::LabelStr() const
  {
    return sigmaSPN().Str();
  }

  std::string SubspaceSpBasis::DebugStr() const
  {
    std::ostringstream os;

    for (int state_index=0; state_index<size(); ++state_index)
      {
        const StateSpBasis state(*this,state_index);

        os << fmt::format(
            "  index {} sigmaSPN {} omega,L {}, {} multiplicity {} offset {}",
            state_index,
            state.sigmaSPN().Str(),state.omega().Str(),state.L(),
            state.degeneracy(),state.offset()
          ) << std::endl;
      }

    return os.str();
  }

  SpaceSpBasis::SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J)
  {
    J_=J;

    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {

        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);

        // create new subspace -- only if not already constructed for this (omega,S)
        const u3shell::U3SPN& sigmaSPN = baby_spncci_subspace.sigmaSPN();
        int irrep_family_index = baby_spncci_subspace.irrep_family_index();

        if(ContainsSubspace(irrep_family_index))
          continue;

        PushSubspace(SubspaceSpBasis(J,sigmaSPN,irrep_family_index,baby_spncci_space));
      }
  }


  SpaceSpBasis::SpaceSpBasis(const BabySpNCCISpace& baby_spncci_space, const HalfInt& J, std::set<int>irrep_family_subset)
  {
    J_=J;

    for(int baby_spncci_subspace_index=0; baby_spncci_subspace_index<baby_spncci_space.size(); ++baby_spncci_subspace_index)
      {
  
        // set up alias
        const BabySpNCCISubspace& baby_spncci_subspace=baby_spncci_space.GetSubspace(baby_spncci_subspace_index);
  
        // create new subspace -- only if not already constructed for this (omega,S)
        const u3shell::U3SPN& sigmaSPN = baby_spncci_subspace.sigmaSPN();
        int irrep_family_index = baby_spncci_subspace.irrep_family_index();
        
        if(not irrep_family_subset.count(irrep_family_index))
          continue;

        if(ContainsSubspace(irrep_family_index))
          continue;

        PushSubspace(SubspaceSpBasis(J,sigmaSPN,irrep_family_index,baby_spncci_space));
      }
  }






  std::string SpaceSpBasis::DebugStr(bool show_subspaces) const
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


// // TODO: FIX. DONT REALLY NEED SECTORS
//   SectorsSpBasis::SectorsSpBasis(
//         const SpaceSpBasis& space_bra,
//         const SpaceSpBasis& space_ket,
//         HalfInt J0,
//         basis::SectorDirection sector_direction
//     )
//     : BaseSectors(space)
//   {
//     HalfInt Jp=spbasis_bra.J();
//     HalfInt J=spbasis_ket.J();
//     if(not am::AllowedTriangle(J,J0,Jp))
//       return;

//     for (int bra_subspace_index=0; bra_subspace_index<space_bra.size(); ++bra_subspace_index)
//       for (int ket_subspace_index=0; ket_subspace_index<space_ket.size(); ++ket_subspace_index)
//         {
//           // enforce canonical ordering
//           if (
//               (sector_direction == basis::SectorDirection::kCanonical)
//               && !(bra_subspace_index<=ket_subspace_index)
//             )
//             continue;

//           // retrieve subspaces
//           const SubspaceType& bra_subspace = space_bra.GetSubspace(bra_subspace_index);
//           const SubspaceType& ket_subspace = space_ket.GetSubspace(ket_subspace_index);

//           // push sector
//           PushSector(bra_subspace_index,ket_subspace_index);
//         }
//   }








  void GetSpBasisOffsets(
    const spncci::SpaceSpBasis& spbasis,
    std::vector<std::vector<int>>& offsets
    )
  {
    offsets.resize(spbasis.size());
    int offset=0;
    for(int subspace_index=0; subspace_index<spbasis.size(); ++subspace_index)
      {
        const spncci::SubspaceSpBasis& subspace=spbasis.GetSubspace(subspace_index);
        offsets[subspace_index].resize(subspace.size());
        const std::vector<int>& multiplicities=subspace.state_multiplicities();
        for(int state_index=0; state_index<subspace.size(); ++state_index)
          {
            offsets[subspace_index][state_index]=offset;
            int multiplicity=multiplicities[state_index];
            offset+=multiplicity;
          }
      }
    assert(spbasis.FullDimension()==offset);
    // std::cout<<"getting offsets "<<std::endl;
    // std::cout<<"full dimension "<<spbasis.FullDimension()<<std::endl;
    // for(int subspace_index=0; subspace_index<spbasis.size(); ++subspace_index)
    //   {
    //     const spncci::SubspaceSpBasis& subspace=spbasis.GetSubspace(subspace_index);
    //     for(int state_index=0; state_index<subspace.size(); ++state_index)
    //       {
    //         std::cout<<subspace_index<<"  "<<state_index<<"  "<<offsets[subspace_index][state_index]<<std::endl;
    //       }
    //   }

  }


  void ConstructOperatorMatrix(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    const HalfInt& J0,
    // u3::WCoefCache& w_cache,
    const spncci::SpaceSpBasis& spbasis_bra, //For a given J
    const spncci::SpaceSpBasis& spbasis_ket, //For a given J
    std::vector<int>& nums_lgi_pairs,
    int observable_index, int hw_index,
    spncci::OperatorBlock& operator_matrix
  )
  {
    // Get dimension of basis
    int basis_size_bra=spbasis_bra.FullDimension();
    int basis_size_ket=spbasis_ket.FullDimension();
    HalfInt Jp=spbasis_bra.J();
    HalfInt J=spbasis_ket.J();

    std::vector<std::vector<int>> offsets_bra;
    std::vector<std::vector<int>> offsets_ket;
    spncci::GetSpBasisOffsets(spbasis_bra,offsets_bra);
    spncci::GetSpBasisOffsets(spbasis_ket,offsets_ket);


    //Set up full matrix
    operator_matrix=spncci::OperatorBlock::Zero(basis_size_bra,basis_size_ket);

    int num_files=nums_lgi_pairs.size();
    for(int thread_num=0; thread_num<num_files; ++thread_num)
      {
        // std::cout<<"thread_num"<<std::endl;
        int num_lgi_pairs=nums_lgi_pairs[thread_num];

        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Open files containing hyperblocks and hypersectors interating over thread number
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        std::ifstream hyperblocks_stream;
        // #pragma omp single
        {
          std::string filename=fmt::format(
              "hyperblocks/observable_hyperblocks_{:02d}_{:02d}_{:02d}.rmes",observable_index,hw_index,thread_num
            );

          std::ios_base::openmode mode_argument = std::ios_base::in | std::ios_base::binary;


          hyperblocks_stream.open(filename,mode_argument);

          // Check if file found
          if(not bool(hyperblocks_stream))
            {
              std::cout<<filename+" not found."<<std::endl;
              assert(hyperblocks_stream);
            }
        }

        std::ifstream hypersectors_stream;
        // #pragma omp single
        {
          std::string filename=fmt::format(
              "hyperblocks/observable_hypersectors_{:02d}_{:02d}.rmes",observable_index,thread_num
            );

          std::ios_base::openmode mode_argument = std::ios_base::in | std::ios_base::binary;

          hypersectors_stream.open(filename,mode_argument);

          // Check if file found
          if(not bool(hypersectors_stream))
            {
              std::cout<<filename+" not found."<<std::endl;
              assert(hypersectors_stream);
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Iterate over lgi pairs stored in given file
        #pragma omp parallel for schedule(dynamic) shared(hyperblocks_stream)
        for(int i=0; i<num_lgi_pairs; ++i)
          {
            // std::cout<<"lgi pair "<<i<<" of "<<num_lgi_pairs<<std::endl;
            //Private Caches
            u3::WCoefCache w_cache;

            spncci::LGIPair lgi_pair;
            std::vector<spncci::ObservableHypersectorLabels> list_baby_spncci_hypersectors;
            // spncci::LGIPair lgi_pair_test;
            basis::OperatorHyperblocks<double> baby_spncci_observable_hyperblocks;

            //Each thread, one at a time, reads in observable hypersectors
            // and hyperblocks for a given lgi pair.
            #pragma omp critical (read_observabl_hyperblocks)
              {
                int num_hypersectors;
                spncci::ReadObservableHypersectors(
                  hypersectors_stream,lgi_pair,
                  list_baby_spncci_hypersectors,num_hypersectors
                  );

                // std::cout<<"reading observable hyperblocks"<<std::endl;
                spncci::ReadObservableHyperblocks(
                  hyperblocks_stream,lgi_pair,
                  baby_spncci_observable_hyperblocks
                );

              }

            int irrep_family_index_bra, irrep_family_index_ket;
            std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
            //////////////////////////////////////////////////////////////////////////////////////////
            //Branching and insert into matrix
            //////////////////////////////////////////////////////////////////////////////////////////
            // std::cout<<"branching to J"<<std::endl;
            for(int observable_hypersector_index=0; observable_hypersector_index<list_baby_spncci_hypersectors.size(); ++observable_hypersector_index)
              {
                int baby_spncci_index_bra, baby_spncci_index_ket, operator_subspace_index, rho0;
                std::tie(baby_spncci_index_bra,baby_spncci_index_ket,operator_subspace_index,rho0)
                          =list_baby_spncci_hypersectors[observable_hypersector_index];

                spncci::OperatorBlock& block=baby_spncci_observable_hyperblocks[observable_hypersector_index][0];
                int block_rows=block.rows();
                int block_cols=block.cols();

                // std::cout<<"Basis information"<<std::endl;
                const spncci::BabySpNCCISubspace& baby_spncci_subspace_bra=baby_spncci_space.GetSubspace(baby_spncci_index_bra);
                const spncci::BabySpNCCISubspace& baby_spncci_subspace_ket=baby_spncci_space.GetSubspace(baby_spncci_index_ket);

                // std::cout<<fmt::format("baby spncci {} {}  {} {}",
                //   irrep_family_index_bra,irrep_family_index_ket,
                //   baby_spncci_subspace_bra.irrep_family_index(),
                //   baby_spncci_subspace_ket.irrep_family_index()
                //   )<<std::endl;

                // Because conjugate included, bra might equal ket and vice versa
                // assert(irrep_family_index_bra==baby_spncci_subspace_bra.irrep_family_index());
                // assert(irrep_family_index_ket==baby_spncci_subspace_ket.irrep_family_index());

                // std::cout<<"Get subspace labels for branching"<<std::endl;
                const u3::U3& omegap=baby_spncci_subspace_bra.omega();
                HalfInt Sp=baby_spncci_subspace_bra.S();
                MultiplicityTagged<int>::vector Lp_list=u3::BranchingSO3(omegap.SU3());

                const u3::U3& omega=baby_spncci_subspace_ket.omega();
                HalfInt S=baby_spncci_subspace_ket.S();
                MultiplicityTagged<int>::vector L_list=u3::BranchingSO3(omega.SU3());

                int spbasis_index_bra=spbasis_bra.LookUpSubspaceIndex(baby_spncci_subspace_bra.irrep_family_index());
                int spbasis_index_ket=spbasis_ket.LookUpSubspaceIndex(baby_spncci_subspace_ket.irrep_family_index());

                const spncci::SubspaceSpBasis& spbasis_subspace_bra=spbasis_bra.GetSubspace(spbasis_index_bra);
                const spncci::SubspaceSpBasis& spbasis_subspace_ket=spbasis_ket.GetSubspace(spbasis_index_ket);

                // const std::vector<int>& state_kappa_max_bra=spbasis_subspace_bra.state_kappa_max();
                // const std::vector<int>& state_kappa_max_ket=spbasis_subspace_ket.state_kappa_max();

                // const std::vector<int>& state_upsilon_max_bra=spbasis_subspace_bra.state_upsilon_max();
                // const std::vector<int>& state_upsilon_max_ket=spbasis_subspace_ket.state_upsilon_max();

                // const std::vector<int>& state_multiplicities_bra=spbasis_subspace_bra.state_multiplicities();
                // const std::vector<int>& state_multiplicities_ket=spbasis_subspace_ket.state_multiplicities();

                // const std::vector<int>& state_offsets_bra=spbasis_subspace_bra.state_offsets();
                // const std::vector<int>& state_offsets_ket=spbasis_subspace_ket.state_offsets();

                // std::cout<<"subspace index "<<spbasis_index_bra<<std::endl;
                // for(int i=0; i<spbasis_subspace_bra.size(); ++i)
                //   std::cout<<i<<"  "<<state_offsets_bra[i]<<std::endl;

                //Operator information
                const u3shell::ObservableSubspaceU3S& observable_subspace=observable_space.GetSubspace(operator_subspace_index);

                int N0,kappa0,L0;
                u3::SU3 x0;
                HalfInt S0;
                std::tie(N0,x0,S0,kappa0,L0)=observable_subspace.Key();
                // std::cout<<"iterate over possible Lp values "<<std::endl;
                for(auto Lp_tagged : Lp_list)
                  {
                    int Lp=Lp_tagged.irrep;

                    int kappap_max=Lp_tagged.tag;
                    // std::cout<<"kappap_max "<<kappap_max<<std::endl;
                    spncci::omegaLLabels state_labels(omegap,Lp);
                    // std::cout<<"looking up state "<<std::endl;
                    // std::cout<<"omegap Lp "<<omegap.Str()<<"  "<<Lp<<std::endl;
                    // std::cout<<spbasis_subspace_bra.DebugStr()<<std::endl;
                    // if(not spbasis_subspace_bra.ContainsState(state_labels))
                    //   continue;



                    // std::cout<<"contained state "<<std::endl;
                    int spbasis_state_index_bra=spbasis_subspace_bra.LookUpStateIndex(state_labels);
                    // std::cout<<"index "<<spbasis_state_index_bra<<std::endl;

                    if(spbasis_state_index_bra==-1)
                      continue;

                    // int upsilon_max=state_upsilon_max_bra[spbasis_state_index_bra];
                    // int multiplicity_bra=spbasis_state_bra.degeneracy()/kappap_max;
                    // std::cout<<"getting offset "<<std::endl;

                    // int index_bra=state_offsets_bra[spbasis_state_index_bra];
                    // std::cout<<"iterate over possible L values "<<std::endl;
                    for(auto L_tagged : L_list)
                      {
                        int L=L_tagged.irrep;

                        if(not am::AllowedTriangle(L,L0,Lp))
                          continue;

                        int kappa_max=L_tagged.tag;
                        int spbasis_state_index_ket=spbasis_subspace_ket.LookUpStateIndex(spncci::omegaLLabels(omega,L));
                        if(spbasis_state_index_ket==-1)
                          continue;

                        // int upsilon_max=state_upsilon_max_ket[spbasis_state_index_ket];

                        double LSJcoef=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp);

                        // std::cout<<"iterate over possible kappap values "<<std::endl;
                        int index_bra=offsets_bra[spbasis_index_bra][spbasis_state_index_bra];
                        // std::cout<<"initial bra "<<spbasis_index_bra<<"  "<<spbasis_state_index_bra
                        //   <<"  "<<index_bra<<std::endl;
                        for(int kappap=1; kappap<=kappap_max; ++kappap)
                        {

                          int index_ket=offsets_ket[spbasis_index_ket][spbasis_state_index_ket];
                          // std::cout<<"initial ket "<<spbasis_index_ket<<"  "
                          //   <<spbasis_state_index_ket<<"  "<<index_ket<<std::endl;
                          // std::cout<<"iterate over possible kappa values "<<kappa_max<<std::endl;
                          for(int kappa=1; kappa<=kappa_max; ++kappa)
                            {
                              spncci::MatrixFloatType Wcoef
                                =u3::WCached(w_cache,omega.SU3(),kappa,L,x0,kappa0,L0,omegap.SU3(),kappap,Lp,rho0);

                              // std::cout<<Wcoef<<"  "<<LSJcoef<<"  "<<std::endl;
                              // std::cout<<block<<std::endl;
                              // std::cout<<"kappas "<<kappap_max<<"  "<<kappa_max<<std::endl;
                              // std::cout<<fmt::format("{} {} {} {}",index_bra,index_ket,block_rows,block_cols)<<std::endl;
                              operator_matrix.block(index_bra,index_ket,block_rows,block_cols)
                                +=Wcoef*LSJcoef*block;

                              //increment starting index for ket
                              index_ket+=block_cols;


                            }
                          //increment starting index for bra
                          index_bra+=block_rows;
                        }

                      }
                  }


            }

            //////////////////////////////////////////////////////////////////////////////////////////

          }// End LGI pair loop
        hyperblocks_stream.close();
        hypersectors_stream.close();
    } //end thread_num

    //Iterate through baby spncci observable blocks
    //Get list of L's to branch to
    // Branch block and insert into correct position in operator matrix
  }






}  // namespace
