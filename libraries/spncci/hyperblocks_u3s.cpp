/****************************************************************
  branching_u3s.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/hyperblocks_u3s.h"

#include <fstream>
#include <iostream>
#include "am/wigner_gsl.h" 
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "mcutils/eigen.h"
// #include "spncci/branching2.h"
#include "spncci/io_control.h"
#include "spncci/unit_tensor.h"

namespace spncci
{
void 
  ContractSymmetricOpBabySpNCCI(
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const u3shell::ObservableSpaceU3S& observable_space,
      const u3shell::RelativeRMEsU3SSubspaces& relative_observable,
      const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
      const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks,
      const spncci::ObservableBabySpNCCIHypersectors& observable_hypersectors,
      basis::OperatorHyperblocks<double>& observable_hyperblocks,
      bool conjugate_hyperblocks
    )
  {
    
    // iterate over baby spncci unit tensors
    // Contract into corresponding observable hyperblock if conjugate_hyperblocks is false
    // Contract conjugated unit tensor into corresponding hyperblock if conjugate_hyperblocks is true

    // Iterate over relative observable rmes to get unit_tensor_subspace_index, kappa0, and L0
    // and use them to look up unit tensors labels.
    // Use labels to look up corresponding index in observable space
    for(auto it=relative_observable.begin(); it!=relative_observable.end(); ++it)
      {
        const std::vector<double>& relative_rmes=it->second;

        // Get unit tensor labels 
        int unit_tensor_subspace_index,kappa0,L0;
        std::tie(unit_tensor_subspace_index,kappa0,L0)=it->first;
        const u3shell::RelativeUnitTensorSubspaceU3S& unit_tensor_subspace=unit_tensor_space.GetSubspace(unit_tensor_subspace_index);

        // Look up observable index
        int etap,eta;
        u3::SU3 x0; 
        HalfInt S0;
        std::tie(x0,S0,etap,eta)=unit_tensor_subspace.labels();
        
        int observable_subspace_index,observable_subspace_index_conj;

        // If not conjugate_hyperblocks, then look up observable index, else look up conjugate index
        if(not conjugate_hyperblocks)
        {
          observable_subspace_index=observable_space.LookUpSubspaceIndex(u3shell::ObservableSubspaceLabels(etap-eta,x0,S0,kappa0,L0));
        
          // std::cout<<"Iterate over baby spncci hypersectors."<<std::endl; 
          for(int baby_spncci_hypersector_index=0; baby_spncci_hypersector_index<baby_spncci_hypersectors.size(); ++baby_spncci_hypersector_index)
            {
              // std::cout<<"get hypersector "<<baby_spncci_hypersector_index<<std::endl;
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
               int observable_hypersector_index
                  =observable_hypersectors.LookUpHypersectorIndex(baby_spncci_index_bra,baby_spncci_index_ket,observable_subspace_index,rho0);
              
              // std::cout<<"Get observable hyperblock and unit tensor hyperblock "<<observable_hypersector_index<<std::endl;
              const basis::OperatorBlocks<double>& unit_tensor_blocks=unit_tensor_hyperblocks[baby_spncci_hypersector_index];
              basis::OperatorBlock<double>& observable_block=observable_hyperblocks[observable_hypersector_index][0];

              // std::cout<<"Loop over unit tensors and contract"<<std::endl;
              for(int unit_tensor_index=0; unit_tensor_index<unit_tensor_subspace.size(); ++unit_tensor_index)
                  observable_block+=relative_rmes[unit_tensor_index]*unit_tensor_blocks[unit_tensor_index];                
            }
        }
            // int observable_hypersector_index, observable_hypersector_index_conj;

        else
          {
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // std::cout<<"Look up conjugate unit tensor subspace index"<<std::endl;
            u3shell::UnitTensorSubspaceLabels unit_tensor_subspace_labels_conj(u3::Conjugate(x0),S0,eta,etap);
            int unit_tensor_subspace_index_conj=unit_tensor_space.LookUpSubspaceIndex(unit_tensor_subspace_labels_conj);

            // std::cout<<"Get conjugate unit tensor subspace"<<std::endl;
            const auto& unit_tensor_subspace_conj=unit_tensor_space.GetSubspace(unit_tensor_subspace_index_conj);

            // std::cout<<"Get conjugate relative observables"<<std::endl;
            std::tuple<int,int,int> rme_labels_conj(unit_tensor_subspace_index_conj,kappa0,L0);
            const std::vector<double>& relative_rmes_conj=relative_observable.at(rme_labels_conj);

            // std::cout<<"Look up conjugate observable subspace index"<<std::endl;
            observable_subspace_index_conj
              =observable_space.LookUpSubspaceIndex(u3shell::ObservableSubspaceLabels(eta-etap,u3::Conjugate(x0),S0,kappa0,L0));
            // std::cout<<"index is "<<observable_subspace_index_conj<<std::endl;
           // std::cout<<"Iterate over baby spncci hypersectors. "<<baby_spncci_hypersectors.size()<<std::endl; 
            for(int baby_spncci_hypersector_index=0; baby_spncci_hypersector_index<baby_spncci_hypersectors.size(); ++baby_spncci_hypersector_index)
              {
                // std::cout<<"get hypersector "<<baby_spncci_hypersector_index<<std::endl;
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
                
                // If Nnp==Nn then sector already included from conjugate hypersectors

                int Nnp=baby_spncci_subspace_bra.Nn();
                int Nn=baby_spncci_subspace_ket.Nn();
                // std::cout<<"Nnp and Nn "<<Nnp<<"  "<<Nn<<std::endl;

                if(Nnp==Nn)
                  continue;

                const basis::OperatorBlocks<double>& unit_tensor_blocks=unit_tensor_hyperblocks[baby_spncci_hypersector_index];
                
                // std::cout<<baby_spncci_index_ket<<"  "<<baby_spncci_index_bra<<"  "<<observable_subspace_index_conj<<"  "<<rho0<<std::endl;
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

void ContractBabySpNCCISymmetricHypersectors(
  const spncci::LGIPair& lgi_pair,
  int num_observables, int num_hw_values,
  const spncci::BabySpNCCISpace& baby_spncci_space,
  const std::vector<u3shell::ObservableSpaceU3S>& observable_spaces,
  const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors1,
  const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors2,
  const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks1,
  const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks2,
  const std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>>& observables_relative_rmes,
  spncci::ObservableHypersectorsTable& observable_hypersectors_table,
  spncci::ObservableHyperblocksTable& observable_hyperblocks_table
  )
  {
    int irrep_family_index_bra,irrep_family_index_ket;
    std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
    // should always be false. may remove
    // bool conjugate_hyperblocks=irrep_family_index_bra<irrep_family_index_ket;
    bool conjugate_hyperblocks=false;
    // assert(not conjugate_hyperblocks);

    for(int observable_index=0; observable_index<num_observables; ++observable_index)
      for(int hw_index=0; hw_index<num_hw_values; ++hw_index)
        {
          
            const u3shell::RelativeRMEsU3SSubspaces& relative_observable
                =observables_relative_rmes[hw_index][observable_index];

            // by observable, by lgi pair
            // std::cout<<"get baby_spncci observable hypersectors"<<std::endl;
            spncci::ObservableBabySpNCCIHypersectors& baby_spncci_observable_hypersectors
                =observable_hypersectors_table[observable_index];

            // std::cout<<"populate hypersectors"<<std::endl;
            bool restrict_upper_triangle=true;
            baby_spncci_observable_hypersectors
              =spncci::ObservableBabySpNCCIHypersectors(
                  baby_spncci_space,observable_spaces[observable_index],
                  irrep_family_index_bra,irrep_family_index_ket,
                  restrict_upper_triangle
                );

            // Get baby spncci observable hyperblocks
            basis::OperatorHyperblocks<double>& baby_spncci_observable_hyperblocks
              =observable_hyperblocks_table[observable_index][hw_index];

            basis::OperatorHyperblocks<double> baby_spncci_observable_hyperblocks_test;
            
            //zero initalize 
            basis::SetHyperoperatorToZero(baby_spncci_observable_hypersectors,baby_spncci_observable_hyperblocks);

            // std::cout<<"Contract over baby spnci observable sectors"<<std::endl;
            // std::cout<<irrep_family_index_bra<<"  "<<irrep_family_index_ket<<std::endl;

            spncci::ContractSymmetricOpBabySpNCCI(
                unit_tensor_space,baby_spncci_space,observable_spaces[observable_index],
                relative_observable,baby_spncci_hypersectors1,unit_tensor_hyperblocks1,
                baby_spncci_observable_hypersectors,baby_spncci_observable_hyperblocks,
                conjugate_hyperblocks
              );  

            // std::cout<<"Contract over baby spnci observable sectors conjugate"<<std::endl;
              spncci::ContractSymmetricOpBabySpNCCI(
                  unit_tensor_space,baby_spncci_space,observable_spaces[observable_index],
                  relative_observable,baby_spncci_hypersectors2,unit_tensor_hyperblocks2,
                  baby_spncci_observable_hypersectors,baby_spncci_observable_hyperblocks, 
                  not conjugate_hyperblocks
                );  
          }
  }

//////////////////////////////////////////
  void ConstructSymmetricOperatorMatrix(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    const HalfInt& J0,
    const spncci::SpaceSpBasis& spbasis_bra, //For a given J
    const spncci::SpaceSpBasis& spbasis_ket, //For a given J
    const std::vector<spncci::LGIPair>& lgi_pairs,
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
    // std::cout<<fmt::format("operator matrix size {}  {}", basis_size_bra,basis_size_ket)<<std::endl;
    operator_matrix=spncci::OperatorBlock::Zero(basis_size_bra,basis_size_ket);

    // int num_files=nums_lgi_pairs.size();
    //Private Caches
    u3::WCoefCache w_cache;
    //TODO break up parallel region and loop to bring w_cache inside parallel region
    // #pragma omp parallel for schedule(dynamic) private(w_cache)
    for(const spncci::LGIPair& lgi_pair : lgi_pairs)
      { 
        int irrep_family_index_bra, irrep_family_index_ket;
        std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
        // std::cout<<"new pair "<<irrep_family_index_bra<<"  "<<irrep_family_index_ket<<std::endl;
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Identify target tile in full matrix 
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        const int spbasis_index_bra=spbasis_bra.LookUpSubspaceIndex(irrep_family_index_bra);
        const int spbasis_index_ket=spbasis_ket.LookUpSubspaceIndex(irrep_family_index_ket);
        
        const spncci::SubspaceSpBasis& spbasis_subspace_bra=spbasis_bra.GetSubspace(spbasis_index_bra);
        const spncci::SubspaceSpBasis& spbasis_subspace_ket=spbasis_ket.GetSubspace(spbasis_index_ket);
        
        const int tile_dimension_bra=spbasis_subspace_bra.dimension();
        const int tile_dimension_ket=spbasis_subspace_ket.dimension();
        
        // std::cout<<"Tile defined as chunk of matrix between states belonging to irrep pair"<<std::endl;
        spncci::OperatorBlock tile=spncci::OperatorBlock::Zero(tile_dimension_bra,tile_dimension_ket);
        const int start_index_bra=tile_dimension_bra>0?offsets_bra[spbasis_index_bra][0]:0;
        const int start_index_ket=tile_dimension_ket>0?offsets_ket[spbasis_index_ket][0]:0;
        

        // std::cout<<"Conjugate block (only lower triangle operator rmes stored in file)"<<std::endl;
        const int spbasis_index_bra_conj=spbasis_bra.LookUpSubspaceIndex(irrep_family_index_ket);
        const int spbasis_index_ket_conj=spbasis_ket.LookUpSubspaceIndex(irrep_family_index_bra);
        
        const spncci::SubspaceSpBasis& spbasis_subspace_bra_conj=spbasis_bra.GetSubspace(spbasis_index_bra_conj);
        const spncci::SubspaceSpBasis& spbasis_subspace_ket_conj=spbasis_ket.GetSubspace(spbasis_index_ket_conj);
        
        const int tile_dimension_bra_conj=spbasis_subspace_bra_conj.dimension();
        const int tile_dimension_ket_conj=spbasis_subspace_ket_conj.dimension();
  
        const int start_index_bra_conj=tile_dimension_bra_conj>0?offsets_bra[spbasis_index_bra_conj][0]:0;
        const int start_index_ket_conj=tile_dimension_ket_conj>0?offsets_ket[spbasis_index_ket_conj][0]:0;


        spncci::OperatorBlock tile_conj=spncci::OperatorBlock::Zero(tile_dimension_bra_conj,tile_dimension_ket_conj);
        spncci::OperatorBlock temp_tile;
  

        // If Jp not equal to J, then accumulate in transposed sector and tranpose to conjugate tile after accumulation finished
        if(Jp!=J)
        {
          temp_tile=spncci::OperatorBlock::Zero(tile_dimension_ket_conj,tile_dimension_bra_conj);
        }

        // std::cout<<fmt::format("({}, {}): ({},{})  ({}, {})  {}x{}",
        //   irrep_family_index_bra,irrep_family_index_ket,spbasis_index_bra,spbasis_index_ket,
        //   start_index_bra,start_index_ket, tile_dimension_bra,tile_dimension_ket
        //   )<<std::endl;

        /////////////////////////////////////////////////////////////////////////////////////////////////
        //// Open files containing hyperblocks and hypersectors interating over thread number
        /////////////////////////////////////////////////////////////////////////////////////////////////
        std::vector<spncci::ObservableHypersectorLabels> list_baby_spncci_hypersectors;
        int num_hypersectors; //Read in from hyperblocks 
        {
          std::ifstream hypersectors_stream;
          std::string filename
            =fmt::format("hyperblocks/observable_hypersectors_{:02d}_{:04d}_{:04d}.rmes",
                observable_index,irrep_family_index_bra,irrep_family_index_ket
              );          
          std::ios_base::openmode mode_argument = std::ios_base::in | std::ios_base::binary;
          hypersectors_stream.open(filename,mode_argument);
          
          // Check if file found
          if(not bool(hypersectors_stream))
            {
              std::cout<<filename+" not found."<<std::endl;
              assert(hypersectors_stream);
            }
          spncci::LGIPair lgi_pair_in;
          spncci::ReadObservableHypersectors(hypersectors_stream,lgi_pair_in,list_baby_spncci_hypersectors,num_hypersectors);
          assert(lgi_pair_in == lgi_pair);
          hypersectors_stream.close();
        }
        
        basis::OperatorHyperblocks<double> baby_spncci_observable_hyperblocks;
        {
          std::ifstream hyperblocks_stream;
          int irrep_family_index_bra_in,irrep_family_index_ket_in;
          std::string filename
            =fmt::format("hyperblocks/observable_hyperblocks_{:02d}_{:02d}_{:04d}_{:04d}.rmes",
              observable_index,hw_index,irrep_family_index_bra,irrep_family_index_ket
            );
              
          std::ios_base::openmode mode_argument = std::ios_base::in | std::ios_base::binary;
          
           
          hyperblocks_stream.open(filename,mode_argument);
          
          // Check if file found
          if(not bool(hyperblocks_stream))
            {
              std::cout<<filename+" not found."<<std::endl;
              assert(hyperblocks_stream);
            }
          
          // std::cout<<"reading observable hyperblocks"<<std::endl;
          spncci::ReadObservableHyperblocks(hyperblocks_stream,lgi_pair,baby_spncci_observable_hyperblocks);
          hyperblocks_stream.close();
        }

        //////////////////////////////////////////////////////////////////////////////////////////
        //Branch and insert into matrix 
        //////////////////////////////////////////////////////////////////////////////////////////
        // std::cout<<"branching to J"<<std::endl;
        for(int observable_hypersector_index=0; 
            observable_hypersector_index<list_baby_spncci_hypersectors.size();
            ++observable_hypersector_index
          )
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
            // std::cout<<fmt::format("irrep family indices  ({}, {})   ({}, {})",irrep_family_index_bra,irrep_family_index_ket,
            //   baby_spncci_subspace_bra.irrep_family_index(),baby_spncci_subspace_ket.irrep_family_index())
            //   <<std::endl;

            // assert(irrep_family_index_bra==baby_spncci_subspace_bra.irrep_family_index());
            // assert(irrep_family_index_ket==baby_spncci_subspace_ket.irrep_family_index());

            // std::cout<<"Get subspace labels for branching"<<std::endl;
            const u3::U3& omegap=baby_spncci_subspace_bra.omega();
            HalfInt Sp=baby_spncci_subspace_bra.S();
            MultiplicityTagged<int>::vector Lp_list=u3::BranchingSO3(omegap.SU3());
            
            const u3::U3& omega=baby_spncci_subspace_ket.omega();
            HalfInt S=baby_spncci_subspace_ket.S();
            MultiplicityTagged<int>::vector L_list=u3::BranchingSO3(omega.SU3());

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
                spncci::omegaLLabels state_labels(omegap,Lp);

                //target state lookup
                int spbasis_state_index_bra=spbasis_subspace_bra.LookUpStateIndex(state_labels);
                int spbasis_state_index_bra_conj=spbasis_subspace_bra_conj.LookUpStateIndex(state_labels);

                if(spbasis_state_index_bra==-1)
                  continue;
                
                // std::cout<<"iterate over possible L values "<<std::endl;
                for(auto L_tagged : L_list)
                  {
                    int L=L_tagged.irrep;

                    if(not am::AllowedTriangle(L,L0,Lp))
                      continue;

                    int kappa_max=L_tagged.tag;

                    //target state lookup
                    int spbasis_state_index_ket=spbasis_subspace_ket.LookUpStateIndex(spncci::omegaLLabels(omega,L));  
                    int spbasis_state_index_ket_conj=spbasis_subspace_ket_conj.LookUpStateIndex(spncci::omegaLLabels(omega,L));  

                    if(spbasis_state_index_ket==-1)
                      continue;

                    double LSJcoef=am::Unitary9J(L,S,J,L0,S0,J0,Lp,Sp,Jp);

                    // offsets_bra returns global indexing.  For accumulation, indexing is relative to
                    int index_bra=offsets_bra[spbasis_index_bra][spbasis_state_index_bra]-start_index_bra;
                    for(int kappap=1; kappap<=kappap_max; ++kappap)
                      {
                        
                        int index_ket=offsets_ket[spbasis_index_ket][spbasis_state_index_ket]-start_index_ket;
                        for(int kappa=1; kappa<=kappa_max; ++kappa)
                          {
                            // std::cout<<"accumulate "<<std::endl;
                            spncci::MatrixFloatType Wcoef
                              =u3::WCached(w_cache,omega.SU3(),kappa,L,x0,kappa0,L0,omegap.SU3(),kappap,Lp,rho0);
   
                             // operator_matrix.block(index_bra,index_ket,block_rows,block_cols)
                              // std::cout<<fmt::format("{}  {}  {}x{}   {}x{}",
                              //   index_bra,index_ket,block_rows,block_cols,tile.rows(),tile.cols()
                              //   )<<std::endl;
                              tile.block(index_bra,index_ket,block_rows,block_cols)+=Wcoef*LSJcoef*block;

                            // std::cout<<"successfully added block"<<std::endl;
                            //increment starting index for ket
                            index_ket+=block_cols;

                          }

                        //increment starting index for bra
                        index_bra+=block_rows;
                      }
                  }
              }//end Lp
          }// end hypersectors 

        // If diagonal sector, only upper triangle was stored in hypersector files
        // need to populate lower triangle of tile
        // std::cout<<tile<<std::endl<<std::endl;
        if(irrep_family_index_bra==irrep_family_index_ket)
        {
          for (int j=0; j<tile.cols(); ++j)
            for(int i=0; i<j; ++i)
              {
                tile(i,j)=tile(j,i);
                // std::cout<<tile(j,i)<<"  "<<tile(i,j)<<std::endl;
              }
        }

        // std::cout<<"Put tile in matrix (or destribute if MPI)"<<std::endl;
        // std::cout<<tile<<std::endl<<std::endl;
        // std::cout<<fmt::format("{} {}   {} {}",start_index_bra,start_index_ket,tile_dimension_bra,tile_dimension_ket)<<std::endl;
        // std::cout<<operator_matrix.rows()<<"  "<<operator_matrix.cols()<<std::endl;
        operator_matrix.block(start_index_bra,start_index_ket,tile_dimension_bra,tile_dimension_ket)
          =tile;
        // std::cout<<"finished plopping"<<std::endl;
        if(Jp==J && irrep_family_index_bra!=irrep_family_index_ket)
        {
          operator_matrix.block(start_index_ket,start_index_bra,tile_dimension_ket,tile_dimension_bra)
            =tile.transpose();
        }

      } //lgi_pair
    // std::cout<<"matrix "<<std::endl<<operator_matrix<<std::endl;
    // for (int i=0; i<basis_size_bra; ++i)
    //   for(int j=0; j<i; ++j)
    //     {
    //       operator_matrix(j,i)=operator_matrix(i,j);
    //     }
  }


}  // namespace
