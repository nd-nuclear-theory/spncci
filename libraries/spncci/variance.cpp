/****************************************************************
  variance.cpp

  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/
#include "spncci/variance.h"

#include <fstream>

#include "fmt/format.h"
#include "spncci/hyperblocks_u3s.h"
#include "spncci/eigenproblem.h"
namespace spncci
{
//////////////////////////////////////////

  void GetVarianceBlock(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    const HalfInt& J0,
    const spncci::SpaceSpBasis& spbasis_bra, //For a given J
    const spncci::SpaceSpBasis& spbasis_ket, //For a given J
    const std::vector<spncci::LGIPair>& lgi_pairs, //Defines tiles to get computed 
    int observable_index, int hw_index,
    spncci::OperatorBlock& operator_matrix
  )
  // For give lgi pairs, compute operator tiles and accumulate in operator matrix
  {
    // Get dimension of basis 
    int basis_size_bra=spbasis_bra.FullDimension();
    int basis_size_ket=spbasis_ket.FullDimension();
    HalfInt Jp=spbasis_bra.J();
    HalfInt J=spbasis_ket.J();

    // Get offsets
    std::vector<std::vector<int>> offsets_bra;
    std::vector<std::vector<int>> offsets_ket;
    spncci::GetSpBasisOffsets(spbasis_bra,offsets_bra);
    spncci::GetSpBasisOffsets(spbasis_ket,offsets_ket);

    //Set up full matrix 
    operator_matrix=spncci::OperatorBlock::Zero(basis_size_bra,basis_size_ket);

    //Private Caches
    u3::WCoefCache w_cache;
    
    // For each pair of lgi, compute operator tile
    #pragma omp parallel for schedule(dynamic) private(w_cache)
    for(int i=0; i<lgi_pairs.size(); ++i)
      { 
        const spncci::LGIPair& lgi_pair=lgi_pairs[i];
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        // Generate irrep pair tile
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        int irrep_family_index_bra, irrep_family_index_ket;
        std::tie(irrep_family_index_bra,irrep_family_index_ket)=lgi_pair;
        const int spbasis_index_bra=spbasis_bra.LookUpSubspaceIndex(irrep_family_index_bra);
        const int spbasis_index_ket=spbasis_ket.LookUpSubspaceIndex(irrep_family_index_ket);
    
        const spncci::SubspaceSpBasis& spbasis_subspace_bra=spbasis_bra.GetSubspace(spbasis_index_bra);
        const spncci::SubspaceSpBasis& spbasis_subspace_ket=spbasis_ket.GetSubspace(spbasis_index_ket);
    
        const std::vector<int>& offsets_bra_subspace=offsets_bra[spbasis_index_bra];
        const std::vector<int>& offsets_ket_subspace=offsets_ket[spbasis_index_ket];
    
        const int tile_dimension_bra=spbasis_subspace_bra.dimension();
        const int tile_dimension_ket=spbasis_subspace_ket.dimension();
        
        // std::cout<<"tile dimensions "<<tile_dimension_bra<<"  "<<tile_dimension_ket<<std::endl;
        if(tile_dimension_ket==0 || tile_dimension_bra==0)
          continue;

        const int start_index_bra=offsets_bra_subspace[0];
        const int start_index_ket=offsets_ket_subspace[0];

        // If bra<ket, then need to first compute adjoint and then get tile by 
        // transposing when inserting into matrix
        bool get_adjoint=irrep_family_index_bra<irrep_family_index_ket;
        if(get_adjoint)
          {
            spncci::LGIPair lgi_pair_adjoint(irrep_family_index_ket,irrep_family_index_bra);
            std::vector<spncci::ObservableHypersectorLabels> list_baby_spncci_hypersectors;
            basis::OperatorHyperblocks<double> baby_spncci_observable_hyperblocks;

            spncci::GetBabySpNCCIHyperBlocks(
              observable_index,hw_index,lgi_pair_adjoint,
              list_baby_spncci_hypersectors,
              baby_spncci_observable_hyperblocks
              );

            spncci::OperatorBlock adjoint_tile;
            spncci::GetOperatorTile(
              baby_spncci_space,observable_space,
              spbasis_subspace_ket,spbasis_subspace_bra,
              offsets_ket_subspace,offsets_bra_subspace,
              J0,J,Jp,hw_index,observable_index,
              lgi_pair_adjoint,w_cache,
              list_baby_spncci_hypersectors,
              baby_spncci_observable_hyperblocks,
              adjoint_tile
            );
            
            operator_matrix.block(start_index_bra,start_index_ket,tile_dimension_bra,tile_dimension_ket)
              =adjoint_tile.transpose();
          }
        // otherwise just compute tile and insert into full matrix 
        else
          {
            std::vector<spncci::ObservableHypersectorLabels> list_baby_spncci_hypersectors;
            basis::OperatorHyperblocks<double> baby_spncci_observable_hyperblocks;
            spncci::GetBabySpNCCIHyperBlocks(
              observable_index,hw_index,lgi_pair,
              list_baby_spncci_hypersectors,
              baby_spncci_observable_hyperblocks
              );

            spncci::OperatorBlock tile;
            spncci::GetOperatorTile(
              baby_spncci_space,observable_space,spbasis_subspace_bra,spbasis_subspace_ket,
              offsets_bra_subspace,offsets_ket_subspace,J0,Jp,J,hw_index,observable_index,
              lgi_pair,w_cache,list_baby_spncci_hypersectors,baby_spncci_observable_hyperblocks,
              tile
            );
          
             operator_matrix.block(start_index_bra,start_index_ket,tile_dimension_bra,tile_dimension_ket)=tile;
          }
      } //lgi_pair
  }

void CalculateVariance(
  const spncci::OperatorBlock& eigenvectors,
  const spncci::OperatorBlock& operator_matrix,
  spncci::OperatorBlock& variance_block
  )
  // Given the operator matrix with rows corresponding to H space and columns corresponding to V spaces
  // and eigenvectors in H space, compute variance given by <psi|H12*H21|psi>
  {
    variance_block=eigenvectors.transpose()*operator_matrix*operator_matrix.transpose()*eigenvectors;
  }


void GetEigensystemH(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    int hw_index, const spncci::RunParameters& run_parameters,
    const std::set<int>& irrep_families_H,
    const std::vector<spncci::SpaceSpBasis>& spbasis_H_byJ,
    std::vector<spncci::Vector>& eigenvalues,  // eigenvalues by J subspace
    std::vector<spncci::Matrix>& eigenvectors  // eigenvectors by J subspace
  )
  // Calculate and diagonalize Hamiltonian matrix in H space
  {
    // std::cout<<"create list of lgi pairs for constructing the Hamiltonian matrix in H space"<<std::endl;
    std::vector<spncci::LGIPair>lgi_pairs_H;
    for(int irrep_family_index_bra : irrep_families_H)
      for(int irrep_family_index_ket : irrep_families_H)
        {
          if (irrep_family_index_bra>=irrep_family_index_ket)
            lgi_pairs_H.emplace_back(irrep_family_index_bra,irrep_family_index_ket);
        }
    
    // observable_index and J0 for Hamiltonian are both zero
    int observable_index_H=0;
    HalfInt J0_H=0;

    // std::cout<<"Set up eigenvector and eigenvalue containers"<<std::endl;
    const std::vector<HalfInt>& Jvalues=run_parameters.J_values;
    eigenvalues.resize(Jvalues.size());  // eigenvalues by J subspace
    eigenvectors.resize(Jvalues.size());  // eigenvectors by J subspace
    
    // std::cout<<"For each J, construct Hamiltonian and get eigensystem"<<std::endl;
    for(int j=0; j<Jvalues.size(); ++j)
      {
        const HalfInt& J=Jvalues[j]; 
        spncci::Vector& eigenvalues_J = eigenvalues[j];
        spncci::Matrix& eigenvectors_J = eigenvectors[j];   
        
        // Get truncated basis branched to J
        const spncci::SpaceSpBasis& spbasis_H=spbasis_H_byJ[j];
        // std::cout<<"dimension of H space "<<spbasis_H.FullDimension()<<" for J="<<J<<std::endl;
        
        // std::cout<<"construct Hamiltonian"<<std::endl;
        spncci::OperatorBlock hamiltonian_matrix;
        spncci::ConstructSymmetricOperatorMatrix(
            baby_spncci_space,observable_space,
            J0_H,spbasis_H,spbasis_H,lgi_pairs_H,
            observable_index_H, hw_index,
            hamiltonian_matrix
          );

        // std::cout<<"Solve Hamiltonian"<<std::endl;
        bool verbose=false;
        spncci::SolveHamiltonian(
            hamiltonian_matrix,
            run_parameters.num_eigenvalues,
            run_parameters.eigensolver_num_convergence,
            run_parameters.eigensolver_max_iterations,
            run_parameters.eigensolver_tolerance,
            eigenvalues_J,eigenvectors_J,verbose
          );
      }//end j
  }//end function


void GetVariances(
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    int observable_index, int hw_index, const HalfInt& J0,
    const std::vector<std::pair<int,int>>& sectors_J,
    const spncci::RunParameters& run_parameters,
    const std::set<int>& irrep_families_H,
    const std::set<int>& irrep_families_V,
    const std::vector<spncci::SpaceSpBasis>& spbasis_H_byJ,
    std::vector<spncci::Matrix>& eigenvectors,  // eigenvectors by J subspace
    std::vector<std::vector<double>>& variances
  )
  // Given H and V spaces and eigenvectors obtained in the H space, get variances for given
  // observable for each valueof J.  
  {
    // Resize variance container
    const std::vector<HalfInt>& Jvalues=run_parameters.J_values;
    variances.resize(Jvalues.size());

    // Create list of lgi pairs for HV block
    std::vector<spncci::LGIPair>lgi_pairs_HV;
    for(int irrep_family_index_bra : irrep_families_H)
      for(int irrep_family_index_ket : irrep_families_V)
          lgi_pairs_HV.emplace_back(irrep_family_index_bra,irrep_family_index_ket);
        
    
    // std::cout<<"Iterate over Jp,J sectors of operator"<<std::endl;
    // Compute variance for each eigenstate with the given Jp
    for(const std::pair<int,int>& J_pair : sectors_J)
      {
        int jp,j; 
        std::tie(jp,j)=J_pair;
        const HalfInt& Jp=Jvalues[jp];
        const HalfInt& J=Jvalues[j];

        // std::cout<<"get J branched basis for H and V space"<<std::endl<<Jp<<"  "<<J<<std::endl;
        const spncci::SpaceSpBasis& spbasis_H=spbasis_H_byJ[jp];
        spncci::SpaceSpBasis spbasis_V=SpaceSpBasis(baby_spncci_space, J, irrep_families_V);
        
        // std::cout<<spbasis_V.DebugStr()<<std::endl;

       //Calculate block of matrix used to calculate variance <H|Op|V>
        spncci::OperatorBlock operator_block;
        spncci::GetVarianceBlock(
          baby_spncci_space, observable_space,J0,
          spbasis_H, spbasis_V, lgi_pairs_HV, 
          observable_index, hw_index,operator_block
        );

        //Calculate the variance
        spncci::OperatorBlock variance_block;
        const spncci::OperatorBlock& eigenvectors_J = eigenvectors[jp];
        CalculateVariance(eigenvectors_J,operator_block,variance_block);

        // Store variances in vector.  Variances are diagonal elements of variance_block
        std::vector<double>& variances_J=variances[jp];
        variances_J.resize(eigenvectors_J.cols());
        for(int r=0; r<eigenvectors_J.cols(); ++r)
          variances_J[r]=variance_block(r,r);        
      }
  }

void GetVariancesForIrrepFamilies(
    const std::vector<int>& irrep_families,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const u3shell::ObservableSpaceU3S& observable_space,
    int observable_index, int hw_index,
    const spncci::RunParameters& run_parameters,
    const std::set<int>& irrep_families_H,
    std::vector<std::vector<std::vector<double>>>& variance_table,
    std::vector<int>& list_irrep_families_V
  )
  // Iterate over all irrep families not in H and compute the variance in the new space
  // formed by H+irrep_family.  Varience stored in variance table index by index of 
  // irrep_family in list_irrep_families_V then by J index
  {
    // Set up SpNCCI basis and sectors for H for each J
    const std::vector<HalfInt>& Jvalues=run_parameters.J_values;
    std::vector<std::pair<int,int>> sectors_J;
    std::vector<spncci::SpaceSpBasis> spbasis_H_byJ(Jvalues.size());
    for(int j=0; j<Jvalues.size(); ++j)
      {
        spbasis_H_byJ[j]=SpaceSpBasis(baby_spncci_space, Jvalues[j], irrep_families_H);
        sectors_J.emplace_back(j,j);
      }

    //Get eigenvalues and eigenvectors for H space
    int J0=0; //for Hamiltonian
    std::vector<spncci::Vector> eigenvalues;  // eigenvalues by J subspace
    std::vector<spncci::Matrix> eigenvectors;  // eigenvectors by J subspace
    spncci::GetEigensystemH(
      baby_spncci_space,observable_space,hw_index,run_parameters,
      irrep_families_H,spbasis_H_byJ,eigenvalues, eigenvectors 
    );

    // Set up V space for each irrep family not in H
    for(int irrep_index: irrep_families)
      {
        if(irrep_families_H.count(irrep_index))
          continue;

        list_irrep_families_V.push_back(irrep_index);                
      }
      
    // std::cout<<"Initializing variance table"<<std::endl;
    variance_table.resize(list_irrep_families_V.size()); // by irrep family, by J
    for(int i=0; i<list_irrep_families_V.size(); ++i)
      variance_table[i].resize(Jvalues.size());

    // std::cout<<"Iterate over different sets of irrep familiy indices defining difference variance subspaces"<<std::endl;
    for(int i=0; i<list_irrep_families_V.size(); ++i)
      {
        int irrep_family_index=list_irrep_families_V[i];
        std::set<int> irrep_families_V;
        irrep_families_V.insert(irrep_family_index);
        std::vector<std::vector<double>>& variances=variance_table[i];
        spncci::GetVariances(
          baby_spncci_space,observable_space,observable_index, hw_index,J0,
          sectors_J,run_parameters,irrep_families_H,irrep_families_V,
          spbasis_H_byJ,eigenvectors, variances
        );
      }
}


void SortIrrepFamiliesByVariance(
    const std::vector<std::vector<std::vector<double>>>& variances,
    const std::vector<int>& irrep_families_V,
    int J_index, int eigenvalue_index,
    std::vector<int>& irrep_families_by_variance,
    double variance_threshold
  )
  // Orders irrep families from largest to smallest variance 
  {
    // Define set which will be used to sort irrep families from largest to smallest variance
    // Set takes pair of values <variance,irrep_family_index>
    //
    // Insert variance, index pairs into set
    std::set<std::pair<double,int>,std::greater<std::pair<double,int>>> irrep_family_sorter;
    for(int i=0; i<irrep_families_V.size(); ++i)
      {
        int irrep_family_index=irrep_families_V[i];
        double variance=variances[i][J_index][eigenvalue_index]; //Only 1 J values 
        
        // If variance is above threshold, include in sorter
        if(variance>variance_threshold)
          {
            std::pair<double,int>key_value(variance,irrep_family_index);
            irrep_family_sorter.insert(key_value); 
          }
      }

    // Use sorter to define list of irrep family indices ordered according to variance
    for(auto& key_value :irrep_family_sorter)
      {
        int irrep_family_index;
        double variance;
        std::tie(variance,irrep_family_index)=key_value;
        irrep_families_by_variance.push_back(irrep_family_index);
        // std::cout<<"family: "<<irrep_family_index<<"  variance: "<<variance<<std::endl;
      }
  }



void SortIrrepFamiliesByVariance(
    const std::vector<std::vector<std::vector<double>>>& variances,
    const std::vector<int>& irrep_families_V,
    int J_index, int eigenvalue_index,
    std::vector<int>& irrep_families_by_variance,
    std::vector<double>& variances_sorted,
    double variance_threshold
  )
  // Orders irrep families from largest to smallest variance 
  {
    // Define set which will be used to sort irrep families from largest to smallest variance
    // Set takes pair of values <variance,irrep_family_index>
    //
    // Insert variance, index pairs into set
    std::set<std::pair<double,int>,std::greater<std::pair<double,int>>> irrep_family_sorter;
    for(int i=0; i<irrep_families_V.size(); ++i)
      {
        int irrep_family_index=irrep_families_V[i];
        double variance=variances[i][J_index][eigenvalue_index]; //Only 1 J values 
        
        // If variance is above threshold, include in sorter
        if(variance>variance_threshold)
          {
            std::pair<double,int>key_value(variance,irrep_family_index);
            irrep_family_sorter.insert(key_value); 
          }
      }

    // Use sorter to define list of irrep family indices ordered according to variance
    for(auto& key_value :irrep_family_sorter)
      {
        int irrep_family_index;
        double variance;
        std::tie(variance,irrep_family_index)=key_value;
        irrep_families_by_variance.push_back(irrep_family_index);
        variances_sorted.push_back(variance);
        // std::cout<<"family: "<<irrep_family_index<<"  variance: "<<variance<<std::endl;
      }
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Removed from spncci main.  TODO: sort through and decide what to keep. 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

void SortIrrepFamiliesByNex(
  const lgi::MultiplicityTaggedLGIVector& lgi_families,
  std::vector<std::vector<int>>& irrep_families_by_Nex,
  int Nmax
  )
  {
    irrep_families_by_Nex.resize(Nmax/2+1);
    for(int index=0; index<lgi_families.size(); ++index)
      {
        int Nex=lgi_families[index].irrep.Nex();
        irrep_families_by_Nex[Nex/2].push_back(index);
      }

    for(auto Nex_set: irrep_families_by_Nex)
    {
      std::cout<<"-----------------"<<std::endl;
      for(auto index : Nex_set)
      {
        std::cout<<"  index: "<<index<<std::endl;
      }
    }
  }


void SortIrrepFamiliesByNex(
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    std::vector<int>& irrep_families,
    std::vector<std::vector<int>>& irrep_families_by_Nex,
    int Nmax
  )
  {
    irrep_families_by_Nex.resize(Nmax/2+1);
    for(int index : irrep_families)
      {
        int Nex=lgi_families[index].irrep.Nex();
        irrep_families_by_Nex[Nex/2].push_back(index);
      }

    for(auto Nex_set: irrep_families_by_Nex)
    {
      std::cout<<"-----------------"<<std::endl;
      for(auto index : Nex_set)
      {
        std::cout<<"  index: "<<index<<std::endl;
      }
    }
  }


void DefineVarianceTruncatedSpace(
    const std::vector<std::vector<std::vector<double>>>& variances,
    const std::vector<std::set<int>>& list_irrep_families_V,
    int eigenvalue_index,
    std::vector<std::vector<int>>& irrep_families_by_variance
)
//DEPRECATED
  {
    irrep_families_by_variance.resize(8); 
    for(int i=0; i<list_irrep_families_V.size(); ++i)
      {
        int irrep_family_index=*(list_irrep_families_V[i].begin()); //Only one irrep_family_index in set 
        double variance=variances[i][0][eigenvalue_index]; //Only 1 J values 
        if(variance>=100)
          irrep_families_by_variance[0].push_back(irrep_family_index);
        else if (variance>=50)
          irrep_families_by_variance[1].push_back(irrep_family_index); 
        else if (variance>=10)
          irrep_families_by_variance[2].push_back(irrep_family_index); 
        else if (variance>=5)
          irrep_families_by_variance[3].push_back(irrep_family_index); 
        else if (variance>=1)          
          irrep_families_by_variance[4].push_back(irrep_family_index);
        else if (variance>=5e-1)
          irrep_families_by_variance[5].push_back(irrep_family_index); 

        else if (variance>1e-1)
          irrep_families_by_variance[6].push_back(irrep_family_index);
        else 
          irrep_families_by_variance[7].push_back(irrep_family_index);

        std::cout<<fmt::format("irrep family variance {:2d}  {:8f}",irrep_family_index, variance)<<std::endl;
      }
    // for(auto& vector :irrep_families_by_variance)
    //   std::cout<<"num irrep families "<<vector.size()<<std::endl;

  }


void TestingVariances(
    const spncci::RunParameters& run_parameters, 
    int hw_index, int J_index, int eigenvalue_index,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const std::vector<u3shell::ObservableSpaceU3S>& observable_spaces,
    const std::vector<int>& irrep_families,
    const std::set<int>& reference_H,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    double variance_threshold
  )
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Testing variance calculation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Get variance for each irrep family with reference wavefunction defined by projection onto dominant irrep 

    // //Set up H space 
    // int dominant_irrep_family_index=3;
    // std::set<int> reference_H; 
    // reference_H.insert(dominant_irrep_family_index);
    

    // Setting up Hamiltonian eigenproblem in H space
    int observable_index=0;
    // int J_index=0;
    // int eigenvalue_index=0;

    std::vector<std::pair<int,int>> sectors_J;
    sectors_J.emplace_back(J_index,J_index);

    //Observable space for Hamiltonian
    const u3shell::ObservableSpaceU3S& observable_space=observable_spaces[observable_index];

    // Compute variance of each irrep family outside of H space 
    std::vector<std::vector<std::vector<double>>> variances;
    std::vector<int> irrep_families_V;
    spncci::GetVariancesForIrrepFamilies(
        irrep_families,baby_spncci_space,observable_space,
        observable_index,hw_index,run_parameters,
        reference_H,variances,irrep_families_V
    );

    //Returns list of irrep families ordered by variance with repsect to reference_H
    // double variance_threshold=10;
    std::vector<int> irrep_families_by_variance_initial;
    std::vector<double> variances_sorted;
    spncci::SortIrrepFamiliesByVariance(
      variances,irrep_families_V,J_index,
      eigenvalue_index,irrep_families_by_variance_initial,
      variances_sorted,variance_threshold
    );

    //TODO: Remove for now.  May want to re-include later
    //Sort by Nsex
    std::vector<std::vector<int>> irrep_families_by_Nex;
    spncci::SortIrrepFamiliesByNex(
      lgi_families,irrep_families_by_variance_initial,
      irrep_families_by_Nex,run_parameters.Nmax
    );

    std::vector<int> irrep_families_by_variance;
    for(auto list : irrep_families_by_Nex)
      for(int index : list)
        irrep_families_by_variance.push_back(index);

    // // Resorting irrep families by Nex, then by variance if reference space increased each iteration
    // std::vector<int> irrep_families_by_variance;
    // std::set<int> irrep_families_H=reference_H;
    // std::vector<int> irrep_families_by_variance_temp;
    // // for(std::vector<int> irrep_families_by_variance_temp : irrep_families_by_Nex)
    // //   { 
    //     int num_iterations=irrep_families_by_variance_initial.size();
    //     irrep_families_by_variance_temp=irrep_families_by_variance_initial;
    //     for(int i=0; i<num_iterations; ++i)
    //       {
    //         int irrep_family_index=irrep_families_by_variance_temp[0];
            
    //         // std::cout<<"Adding family"<<irrep_family_index<<std::endl;
    //         //Add to reference subspaces
    //         irrep_families_H.insert(irrep_family_index);

    //         irrep_families_by_variance.push_back(irrep_family_index);
          
    //         std::vector<std::vector<std::vector<double>>> variances_temp;
    //         std::vector<int> individual_irrep_families_V_temp;
    //         spncci::GetVariancesForIrrepFamilies(
    //           irrep_families,baby_spncci_space,observable_space,
    //           observable_index,hw_index,run_parameters,
    //           irrep_families_H,variances_temp,individual_irrep_families_V_temp
    //         );

    //         // std::cout<<"Sorting by variance "<<std::endl;
    //         irrep_families_by_variance_temp.resize(0);
            
    //         spncci::SortIrrepFamiliesByVariance(
    //           variances_temp, individual_irrep_families_V_temp,J_index,
    //           eigenvalue_index,irrep_families_by_variance_temp
    //         );
    //       }
    //   // }
  

    //   ////////////////////////////////////////////////////////////////////////////////////////////
    //   // reinitializing irrep_families_H with just the dominant irrep
    //   irrep_families_H=reference_H; 
    //   std::cout<<"Starting with irrep families:";
    //   for (auto it=irrep_families_H.begin(); it != irrep_families_H.end(); ++it) 
    //     std::cout << ' ' << *it; 
    //   std::cout<<std::endl;
    //   // std::cout<<"Dominant irrep family index "<<dominant_irrep_family_index<<std::endl;
    //   // irrep_families_H.insert(dominant_irrep_family_index);

      // For each irrep family index in with non-zero variance, add to H space one by one
      // and compute energies and variances
      // std::vector<std::pair<double,double>>variences_for_irrep_families(irrep_families_by_variance.size());
    std::vector<std::pair<double,double>>variences_for_irrep_families(irrep_families_by_variance_initial.size());
      
      //Setting irrep families H to intial reference H 
      std::set<int> irrep_families_H=reference_H;
      for(int i=0; i<irrep_families_by_variance.size(); ++i)
      // for(int i=0; i<irrep_families_by_variance_initial.size(); ++i)
        {
          // int irrep_family_index=irrep_families_by_variance_initial[i];
          int irrep_family_index=irrep_families_by_variance[i];
          irrep_families_H.insert(irrep_family_index);
          // std::cout<<"---------------------------------------------"<<std::endl;
          // std::cout<<"irrep family index "<<irrep_family_index<<std::endl;
          // set up up V space. In this case, there is only 1 V space. 

          const std::vector<HalfInt>& Jvalues=run_parameters.J_values;
          std::vector<spncci::SpaceSpBasis> spbasis_H_byJ(Jvalues.size());
          for(int j=0; j<Jvalues.size(); ++j)
            spbasis_H_byJ[j]=spncci::SpaceSpBasis(baby_spncci_space, Jvalues[j], irrep_families_H);


          //Get eigenvalues and vectors for Hamiltonian in H subspaces
          std::vector<spncci::Vector> eigenvalues;  // eigenvalues by J subspace
          std::vector<spncci::Matrix> eigenvectors;  // eigenvectors by J subspace
          spncci::GetEigensystemH(
            baby_spncci_space,observable_space,hw_index,run_parameters,
            irrep_families_H,spbasis_H_byJ,eigenvalues, eigenvectors 
          );

          std::set<int> irrep_families_V;
          for(int index=0; index<lgi_families.size(); ++index)
            {
              // int index=*(subspace.begin()); //Only one irrep_family_index in set 
              if( not irrep_families_H.count(index))
                irrep_families_V.insert(index);
            }

          if(irrep_families_V.size()==0)
            {
              std::cout<<"H space is full space "<<std::endl;
              continue;
            }
          int J0=0;
          std::vector<std::vector<double>> variances(run_parameters.J_values.size());
          spncci::GetVariances(
            baby_spncci_space,observable_space,observable_index, hw_index,J0,
            sectors_J,run_parameters,irrep_families_H,irrep_families_V,
            spbasis_H_byJ,eigenvectors, variances
          );


          double eigenvalue=eigenvalues[J_index][eigenvalue_index];
          double variance=variances[J_index][eigenvalue_index];
          variences_for_irrep_families[i]=std::pair<double,double>(eigenvalue,variance);

        }
      std::ofstream out_file;
      out_file.open("variances.dat");

      // for(int i=0; i<irrep_families_by_variance_initial.size(); ++i)
      for(int i=0; i<irrep_families_by_variance.size(); ++i)
        {
          // int irrep_family_index=irrep_families_by_variance_initial[i];
          int irrep_family_index=irrep_families_by_variance[i];
          double individual_variance=variances_sorted[i];
          double eigenvalue,variance;
          std::tie(eigenvalue,variance)=variences_for_irrep_families[i];
          //When doing Nex sorted variances, individual variance not correct
          out_file<<fmt::format("{:3d}  {:8.4f}  {:8.4f}  {:8.4f}",irrep_family_index,individual_variance,eigenvalue,variance)<<std::endl;
          std::cout<<fmt::format("{:3d}  {:8.4f}  {:8.4f}  {:8.4f}",irrep_family_index,individual_variance,eigenvalue,variance)<<std::endl;
        }
      out_file.close();
  }



}//end namespace