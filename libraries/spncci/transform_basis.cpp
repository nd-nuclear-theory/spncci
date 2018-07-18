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
    int num_irrep_families,
    int num_eigenvalues,
    const std::vector<spncci::Matrix>& eigenvectors,
    std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks
  )
	// irrep_family_blocks: by J, by n, by irrep family
  {
    std::cout<<"regrouping"<<std::endl;
    //By J, by irrep family index: dimensions of subspace (gamma_max, num_states)
    std::vector<std::vector<std::pair<int,int>>> irrep_family_subspaces(spj_space.size());

    std::cout<<"geting dimensions"<<std::endl;
    // iterate over j spaces and count number of states in each subspace and get gamma_max
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
      // for each J subspace
      {
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
       	
       	//resize container
        irrep_family_subspacesJ.resize(num_irrep_families);
        
        for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
          {
            // retrieve basis state information
            StateSpJ spj_state(spj_subspace,spj_state_index);
            int degeneracy = spj_state.degeneracy();
            int irrep_family_index = spj_state.irrep_family_index();
            int gamma_max=spj_state.gamma_max();
            int upsilon_max=int(degeneracy/gamma_max);
            // std::cout<<degeneracy<<"  "<<irrep_family_index<<"  "<<gamma_max<<"  "<<upsilon_max<<std::endl;
            auto& dimensions=irrep_family_subspacesJ[irrep_family_index];
            // assert(dimensions.first==gamma_max);
            irrep_family_subspacesJ[irrep_family_index]
              =std::pair<int,int>(gamma_max,dimensions.second+upsilon_max);
          }
      }

    std::cout<<"initialize irrep family subspace blocks "<<std::endl;
    irrep_family_blocks.resize(spj_space.size());
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)    
      {
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
        
      	// resize by number of eigenvalues at each J
      	irrep_family_blocks[spj_subspace_index].resize(num_eigenvalues);
      	
      	// For each degenerate eigenstate of J
      	for(int n=0; n<num_eigenvalues; ++n)
      		{
      			//Resize block
						spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_subspace_index][n];
        		blocks.resize(irrep_family_subspacesJ.size());

        		//Zero initilize each block if gamma_max > 1
		        for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
		          {
		            int gamma_max,num_states;
		            std::tie(gamma_max,num_states)=irrep_family_subspacesJ[irrep_family_index];
		            // std::cout<<irrep_family_index<<"  "<<gamma_max<<"  "<<num_states<<std::endl;
		            if(gamma_max>1)
		  	          blocks[irrep_family_index]=Eigen::MatrixXd::Zero(num_states,gamma_max);
		          }  
		      }      
      }

    std::cout<<"populating blocks"<<std::endl;
    for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
      {
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
        const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
        const spncci::Matrix& eigenvectors_J = eigenvectors[spj_subspace_index];
        const int num_eigenvectors = eigenvectors_J.cols();

        for(int n=0; n<num_eigenvalues; ++n)
        	{
		        spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_subspace_index][n];
		        
		        // initialize offsets to zero for each irrep family
		        std::vector<int> offsets(irrep_family_subspacesJ.size(),0);
		        
		        for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
		          {
		            // retrieve basis state information
		            StateSpJ spj_state(spj_subspace,spj_state_index);
		            int gamma_max=spj_state.gamma_max();

		            if(gamma_max<2)
		            	continue;

		            int degeneracy = spj_state.degeneracy();
		            int irrep_family_index = spj_state.irrep_family_index();
		            
		            int upsilon_max=degeneracy/gamma_max;
		            int eigen_offset=spj_state.offset();
		            int offset=offsets[irrep_family_index];
		            for(int gamma=1; gamma<=gamma_max; ++gamma)
		              {
		                //Taking the lowest energy eigenvector for each J
		                //TODO: make J, and which eigenvector of J optional
		                blocks[irrep_family_index].block(offset,gamma-1,upsilon_max,1)
		                  =eigenvectors_J.block(eigen_offset,0,upsilon_max,1);
		                
		                // std::cout<<eigenvectors_J.block(eigen_offset,0,upsilon_max,1)<<std::endl<<std::endl;
		                // std::cout<<blocks[irrep_family_index]<<std::endl<<std::endl;

		                // Increment offset in eigenvector
		                eigen_offset+=upsilon_max;
		              }

		            // Increment offset in irrep family block 
		            offsets[irrep_family_index]+=upsilon_max;
		            
		           //  if(irrep_family_index==6)
		           //  {
		           //  std::cout<<blocks[irrep_family_index]<<std::endl<<std::endl;
		           //  std::cout<<"--------------------------------"<<std::endl<<std::endl;
		          	// }
		          }
		      }
      }
  }


  void WriteIrrepFamilyBlocks(  
    const spncci::SpaceSpJ& spj_space,
    int num_irrep_families,
    int num_eigenvalues,
    const std::vector<int>& lgi_full_space_index_lookup,
    const std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
    const std::string& filename
  )
  {
    std::ios_base::openmode mode_argument = std::ios_base::out | std::ios::app | std::ios_base::binary;
    std::ofstream out_file;
    out_file.open(filename,mode_argument);

    if (!out_file)
      {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
      }

    // Number of J values, num of eigenstates with given J eigenvalue 
    // and num irrep families in full space
    mcutils::WriteBinary<int>(out_file,lgi::binary_float_precision);
    mcutils::WriteBinary<int>(out_file,spj_space.size());
    mcutils::WriteBinary<int>(out_file,num_eigenvalues);
    mcutils::WriteBinary<int>(out_file,num_irrep_families);

    //for each irrep family
		for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
			{
		    // Write irrep family index in full space 
    		int full_space_irrep_family_index=lgi_full_space_index_lookup[irrep_family_index];
    		mcutils::WriteBinary<int>(out_file,full_space_irrep_family_index);
				
		    //for reach J eigenvalue
		    for( int j_index=0; j_index<spj_space.size(); ++j_index)
			    {
						// get J
						HalfInt J = spj_space.GetSubspace(j_index).J();
		      	
		      	// get number of rows and columns
						const spncci::OperatorBlock& block=irrep_family_blocks[j_index][0][irrep_family_index];
						int rows=block.rows();
						int cols=block.cols();

						// Write block information
						mcutils::WriteBinary<int>(out_file,TwiceValue(J));
		    		mcutils::WriteBinary<int>(out_file,rows);
		    		mcutils::WriteBinary<int>(out_file,cols);

		    		// Only write if gamma_max (corresponding to rows) is >1
						if(rows<2)
							continue;

		    		int size=rows*cols;

		      	// For each n value
		      	for(int n=0; n<num_eigenvalues; ++n)
		      		{
				        const spncci::OperatorBlock& block=irrep_family_blocks[j_index][n][irrep_family_index];
				        // write matrix.  Order is column major (Eigen default)
				        if(lgi::binary_float_precision==4)
				          {
				            Eigen::MatrixXf buffer_matrix=block.cast<float>();
				            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);
				            
				          }  
				          
				        else if (lgi::binary_float_precision==8)
				          {
				            Eigen::MatrixXd buffer_matrix=block;
				            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);

				          }
				      }
					}
			}
    out_file.close();    
  }



void ReadIrrepFamilyBlocks(
	std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,//by irrep family, by J, by n
    const std::string& filename
  )
{

  std::ios_base::openmode mode_argument = std::ios_base::in | std::ios::app | std::ios_base::binary;
  std::ifstream in_stream;
  in_stream.open(filename,mode_argument);

  int num_J_values,num_irrep_families,num_eigenvalues, binary_float_precision;
  mcutils::ReadBinary<int>(in_stream,binary_float_precision);  
  mcutils::ReadBinary<int>(in_stream,num_J_values);  
  mcutils::ReadBinary<int>(in_stream,num_eigenvalues);
  mcutils::ReadBinary<int>(in_stream,num_irrep_families);
	
	//for each irrep family
	for(int i=0; i<num_irrep_families; ++i)
		{
			//Read irrep_family index (corresponds to full space)
			int irrep_family_index;
  		mcutils::ReadBinary<int>(in_stream,irrep_family_index);
		
			std::vector<spncci::OperatorBlocks>& blocks=irrep_family_blocks[irrep_family_index];
		  for( int j_index=0; j_index<num_J_values; ++j_index)
		  	{
					int twice_J, rows, cols;
					mcutils::ReadBinary<int>(in_stream,twice_J);
				  mcutils::ReadBinary<int>(in_stream,rows);
					mcutils::ReadBinary<int>(in_stream,cols);

						if(rows<2)
							continue;

		      	// For each n value
		      	for(int n=0; n<num_eigenvalues; ++n)
		      		{
				        spncci::OperatorBlock& block=irrep_family_blocks[j_index][n][irrep_family_index];
				        
				        // Read matrix.  Order is column major (Eigen default)
					      if(lgi::binary_float_precision==4)
					        {
					          float buffer[rows*cols];
					          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
					          block=Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
					        }
					      else if (lgi::binary_float_precision==8)
					        {
					          double buffer[rows*cols];
					          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
					          block=Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
					        }
					    }
		  	}
    }
// assert(0);
}




void  DefineIrrepFamilyTransformation(
  const spncci::SpaceSpJ& spj_space,
  std::vector<spncci::OperatorBlocks>& irrep_family_blocks,
  std::vector<spncci::OperatorBlocks>& transformations
)
{
  std::cout<<"defining rotation"<<std::endl;
  int max_iterations=5;
  int num_random_test=100;
  transformations.resize(spj_space.size());
  for( int j_index=0; j_index<spj_space.size(); ++j_index)
    {
      spncci::OperatorBlocks& blocks=irrep_family_blocks[j_index];
      std::cout<<"blocks size "<<blocks.size()<<std::endl;
      spncci::OperatorBlocks& transformationsJ=transformations[j_index];
      transformationsJ.resize(blocks.size());

      for(int irrep_family_index=0; irrep_family_index<blocks.size(); ++irrep_family_index)
        {
          spncci::OperatorBlock block=blocks[irrep_family_index].transpose();
          int gamma_max=block.rows();
          int irrep_dim=block.cols();
          
          if(gamma_max<2)
            continue;

          // Normalize irrep family block
          double norm=block.squaredNorm();
          block=block/std::sqrt(norm);
          std::cout<<"block norm "<<block.squaredNorm()<<std::endl;
         // if(irrep_family_index!=6)
         // 	continue;
          
          // Maximum probablity of a single irrep
          double max_probability=block.rowwise().squaredNorm().maxCoeff();
          std::cout<<"irrep_family_index  "<<irrep_family_index<<"  "<<gamma_max<<"  "<<norm<<std::endl;
          std::cout<<"initial max probability "<<max_probability<<std::endl;
          
          spncci::OperatorBlock& transformation_matrix=transformationsJ[irrep_family_index];
          transformation_matrix=Eigen::MatrixXd::Identity(gamma_max,gamma_max);
          
          // std::cout<<"constructing matrix"<<std::endl;
          // Construction transformation matrix
          for(int iteration=0; iteration<=max_iterations; ++iteration)
            {
              //for each iteration, find best transformation
              spncci::OperatorBlock max_eigenvectors=Eigen::MatrixXd::Identity(gamma_max,gamma_max);
              for(int test=0; test<=num_random_test; ++test)
                {
                  //Set up text sub-block
                  spncci::OperatorBlock sub_block=Eigen::MatrixXd::Zero(gamma_max,gamma_max);
                  for(int gamma=1; gamma<=gamma_max; ++gamma)
                    {
                      int column=std::experimental::randint(0,irrep_dim-1); //I think its inclusive
                      sub_block.block(0,gamma-1,gamma_max,1)=block.block(0,column,gamma_max,1);

                    }
                  
                  // Get eigen-vectors of sub-block
                  // May need to switch to complex eigenvectors and square
                  Eigen::EigenSolver<spncci::OperatorBlock> eigensystem(sub_block);
                  spncci::OperatorBlock eigenvectors=eigensystem.pseudoEigenvectors(); 
                  // std::cout<<"eigenvectors "<<std::endl;
                  for(int i=0; i<eigenvectors.cols(); ++i)
                  	eigenvectors.col(i).normalize();
                  
                  // std::cout<<eigenvectors.colwise().squaredNorm()<<std::endl;                 
                  spncci::OperatorBlock temp_block=eigenvectors*eigenvectors.transpose()*block;
                  // Normalize irrep family block
                  double temp_norm=temp_block.squaredNorm();
                  // std::cout<<"temp_norm "<<temp_norm<<std::endl;
                  temp_block=temp_block/std::sqrt(temp_norm);
                  // std::cout<<"norm "<<temp_norm<<std::endl;
                  // Maximum probablity of a single irrep
                  double temp_max_probability=temp_block.rowwise().squaredNorm().maxCoeff();
                  // std::cout<<temp_block<<std::endl;
                  // std::cout<<"max probablity temp "<<temp_max_probability<<"  "<<max_probability<<std::endl;
                  // std::cout<<temp_block.rowwise().squaredNorm()<<std::endl;
                  // std::cout<<"-----------------------------"<<std::endl<<std::endl;
                  if(temp_max_probability>max_probability)
                    {
                       // std::cout<<"new max"<<std::endl;
                      max_probability=temp_max_probability;
                      max_eigenvectors=eigenvectors*eigenvectors.transpose()/std::sqrt(temp_norm);
                      // std::cout<<temp_block<<std::endl<<std::endl;
                      // std::cout<<max_eigenvectors<<std::endl<<std::endl;
                      // std::cout<<max_eigenvectors*block<<std::endl<<std::endl;
                    }

                }

              //Transform block
              block=max_eigenvectors*block;
              // std::cout<<"**********"<<std::endl<<block.rowwise().squaredNorm().maxCoeff()<<std::endl;
              // std::cout<<block<<std::endl;
              //accumulate transformations
              transformation_matrix=max_eigenvectors*transformation_matrix;
            }

        	max_probability=block.rowwise().squaredNorm().maxCoeff();
        	std::cout<<"final max probability "<<max_probability<<std::endl<<std::endl;
        	std::cout<<block.rowwise().squaredNorm()<<std::endl;
        	std::cout<<"---------------------------"<<std::endl<<std::endl;

        }
    }
  // assert(0);
}

  void WriteTransformationMatrices(  
  	const spncci::SpaceSpJ& spj_space,
  	int num_irrep_families,
  	const std::vector<spncci::OperatorBlocks>& transformations,
    const std::string& filename
    )
  {
    std::ios_base::openmode mode_argument = std::ios_base::out | std::ios::app | std::ios_base::binary;
    std::ofstream out_file;
    out_file.open(filename,mode_argument);

    if (!out_file)
      {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
      }

    for( int j_index=0; j_index<spj_space.size(); ++j_index)
			for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
			{
				const spncci::OperatorBlock& transformation_matrix=transformations[j_index][irrep_family_index];
				int rows=transformation_matrix.rows();
				if(rows<2)
					continue;

		    // Write irrep family index and rows (should be same as number of columns)
    		mcutils::WriteBinary<int>(out_file,irrep_family_index);
    		mcutils::WriteBinary<int>(out_file,rows);

    		int size=rows*rows;

        // write matrix.  Order is column major (Eigen default)
        if(lgi::binary_float_precision==4)
          {
            Eigen::MatrixXf buffer_matrix=transformation_matrix.cast<float>();
            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);
            
          }  
          
        else if (lgi::binary_float_precision==8)
          {
            Eigen::MatrixXd buffer_matrix=transformation_matrix;
            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);

          }
			}

    out_file.close();    
    

  }

// void ReadObservableHyperblocks(
// 	  int num_irrep_families,
//   	const std::vector<spncci::OperatorBlocks>& transformations,
//     const std::string& filename
//   )

// {
  
//   // Read lgi family 
//   int irrep_family_index_bra, irrep_family_index_ket, num_hyperblocks;
//   mcutils::ReadBinary<int>(in_stream,irrep_family_index_bra);
//   mcutils::ReadBinary<int>(in_stream,irrep_family_index_ket);
//   lgi_pair=spncci::LGIPair(irrep_family_index_bra,irrep_family_index_ket);

//   // std::cout<<irrep_family_index_bra<<"  "<<irrep_family_index_ket<<"  "<<omp_get_thread_num()<<std::endl;
//   // Read number of hyperblocks 
//   mcutils::ReadBinary<int>(in_stream,num_hyperblocks);
//   baby_spncci_obserable_hyperblocks.resize(num_hyperblocks);
//   // std::cout<<" Read in each hyperblock "<< num_hyperblocks<<std::endl;
//   for(int hypersector_index=0; hypersector_index<num_hyperblocks; ++hypersector_index)
//     {
//       baby_spncci_obserable_hyperblocks[hypersector_index].resize(1);
//       spncci::ReadBlock(in_stream, baby_spncci_obserable_hyperblocks[hypersector_index][0]);
//     }
// assert(0);
// }

}