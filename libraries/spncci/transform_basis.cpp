/****************************************************************
  transform_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

****************************************************************/
#include "spncci/transform_basis.h"

#include <experimental/random> //For DefineIrrepFamilyRotation
#include <fstream>
#include <omp.h>  

#include "lgi/lgi_unit_tensors.h"
#include "mcutils/eigen.h"
#include "mcutils/io.h"
#include "spncci/parameters.h"
#include "spncci/spncci_common.h"


namespace spncci
{
void RegroupIntoIrrepFamilies(
    const spncci::SpaceSpJ& spj_space,
    const std::vector<spncci::Matrix>& eigenvectors,
    int num_irrep_families,
    std::vector<spncci::OperatorBlocks>& irrep_family_blocks
  )
  {
    std::cout<<"regrouping"<<std::endl;
    //By J, by irrep family index: dimensions of subspace (gamma_max, num_states)
    std::vector<std::vector<std::pair<int,int>>> irrep_family_subspaces(spj_space.size());

    std::cout<<"geting dimensions"<<std::endl;
    // iterate over j spaces and count number of states in each subspace and get gamma_max
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
      // for each J subspace
      {
        // set up aliases for current J subspace
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        //aliase for J-irrepfamily subspaces 
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
        irrep_family_subspacesJ.resize(num_irrep_families);
        for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
          {
            // retrieve basis state information
            StateSpJ spj_state(spj_subspace,spj_state_index);
            int degeneracy = spj_state.degeneracy();
            int irrep_family_index = spj_state.irrep_family_index();
            int gamma_max=spj_state.gamma_max();
            int upsilon_max=int(degeneracy/gamma_max);
            std::cout<<degeneracy<<"  "<<irrep_family_index<<"  "<<gamma_max<<"  "<<upsilon_max<<std::endl;
            auto& dimensions=irrep_family_subspacesJ[irrep_family_index];
            // assert(dimensions.first==gamma_max);
            irrep_family_subspacesJ[irrep_family_index]
              =std::pair<int,int>(gamma_max,dimensions.second+upsilon_max);
          }
      }

    std::cout<<"seting up new blocks"<<std::endl;
    //initialize irrep family subspace blocks 
    irrep_family_blocks.resize(spj_space.size());
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)    
      {
        // Aliase subspaces
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];

        //Resize vector of blocks 
        spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_subspace_index];
        blocks.resize(irrep_family_subspacesJ.size());

        for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
          {
            int gamma_max,num_states;
            std::tie(gamma_max,num_states)=irrep_family_subspacesJ[irrep_family_index];
            blocks[irrep_family_index]=spncci::OperatorBlock(num_states,gamma_max);
          }        
      }

    // Populate blocks 
    std::cout<<"populating blocks"<<std::endl;
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
      {
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        const spncci::Matrix& eigenvectors_J = eigenvectors[spj_subspace_index];
        const int num_eigenvectors = eigenvectors_J.cols();
        spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_subspace_index];
        
        // initialize offsets to zero for each irrep family
        std::vector<int> offsets(irrep_family_subspacesJ.size(),0);
        
        for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
          {
            // retrieve basis state information
            StateSpJ spj_state(spj_subspace,spj_state_index);
            int degeneracy = spj_state.degeneracy();
            int irrep_family_index = spj_state.irrep_family_index();
            int gamma_max=spj_state.gamma_max();
            int upsilon_max=degeneracy/gamma_max;
            int eigen_offset=spj_state.offset();
            int offset=offsets[irrep_family_index];
            for(int gamma=1; gamma<=gamma_max; ++gamma)
              {
                //Taking the lowest energy eigenvector for each J
                //TODO: make J, and which eigenvector of J optional
                blocks[irrep_family_index].block(0,gamma-1,upsilon_max,1)
                  =eigenvectors_J.block(eigen_offset,0,upsilon_max,1);
                
                // Increment offsets
                eigen_offset+=upsilon_max;
                offsets[irrep_family_index]+=upsilon_max;
              }
            
          }
      }
  }

void  DefineIrrepFamilyTransformation(
  const spncci::SpaceSpJ& spj_space,
  std::vector<spncci::OperatorBlocks>& irrep_family_blocks
)
{
  std::cout<<"defining rotation"<<std::endl;
  int max_iterations=5;
  int num_random_test=10;
  for( int j_index=0; j_index<spj_space.size(); ++j_index)
    {
      spncci::OperatorBlocks& blocks=irrep_family_blocks[j_index];
      std::cout<<"blocks size "<<blocks.size()<<std::endl;
      for(int irrep_family_index=0; irrep_family_index<blocks.size(); ++irrep_family_index)
        {
          spncci::OperatorBlock block=blocks[irrep_family_index].transpose();
          int gamma_max=block.rows();
          int irrep_dim=block.cols();
          
          if(gamma_max==1)
            continue;

          // Normalize irrep family block
          double norm=block.squaredNorm();
          block=block/std::sqrt(norm);
          
          // Maximum probablity of a single irrep
          int max_probability=block.rowwise().squaredNorm().maxCoeff();
          spncci::OperatorBlock transformation_matrix=Eigen::MatrixXd::Identity(gamma_max,gamma_max);
          
          std::cout<<"constructing matrix"<<std::endl;
          // Construction transformation matrix
          for(int iteration=0; iteration<=max_iterations; ++iteration)
            {
              //for each iteration, find best transformation
              spncci::OperatorBlock max_eigenvectors=Eigen::MatrixXd::Identity(gamma_max,gamma_max);
              for(int test=0; test<=num_random_test; ++test)
                {
                  //Set up text sub-block
                  spncci::OperatorBlock sub_block(gamma_max,gamma_max);
                  for(int gamma=1; gamma<=gamma_max; ++gamma)
                    {
                      int column=std::experimental::randint(0,irrep_dim-1); //I think its inclusive
                      sub_block.block(0,gamma-1,gamma_max,1)=block.block(0,column,gamma_max,1);

                    }
                  
                  // Get eigen-vectors of sub-block
                  // May need to switch to complex eigenvectors and square
                  Eigen::EigenSolver<spncci::OperatorBlock> eigensystem(sub_block);
                  spncci::OperatorBlock eigenvectors=eigensystem.pseudoEigenvectors();                  

                  spncci::OperatorBlock temp_block=eigenvectors*block;
                  // Normalize irrep family block
                  double temp_norm=temp_block.squaredNorm();
                  temp_block=temp_block/std::sqrt(temp_norm);
                  std::cout<<"norm "<<temp_norm<<std::endl;
                  // Maximum probablity of a single irrep
                  int temp_max_probability=block.rowwise().squaredNorm().maxCoeff();

                  if(temp_max_probability>max_probability)
                    {
                      max_probability=temp_max_probability;
                      max_eigenvectors=eigenvectors/temp_norm;
                    }

                }

              //Transform block
              block=max_eigenvectors*block;
              
              //accumulate transformations
              transformation_matrix=max_eigenvectors*transformation_matrix;
            }
        max_probability=block.rowwise().squaredNorm().maxCoeff();
        std::cout<<"max probability "<<max_probability<<std::endl;

        }
    }
}
}