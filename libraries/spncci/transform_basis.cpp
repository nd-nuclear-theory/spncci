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
		                //Taking the nth eigenvector for the given J value
		                blocks[irrep_family_index].block(offset,gamma-1,upsilon_max,1)
		                  =eigenvectors_J.block(eigen_offset,n,upsilon_max,1);

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
    std::ios_base::openmode mode_argument = std::ios_base::out | std::ios_base::binary;
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

						// Write block information, col and rows are transposed
						mcutils::WriteBinary<int>(out_file,TwiceValue(J));
		    		mcutils::WriteBinary<int>(out_file,cols);
		    		mcutils::WriteBinary<int>(out_file,rows);

		    		// Only write if gamma_max (corresponding to rows) is >1
						if(cols<2)
							continue;

		    		int size=rows*cols;

		      	// For each n value
		      	// Transpose matrix and write to file
		      	for(int n=0; n<num_eigenvalues; ++n)
		      		{
				        const spncci::OperatorBlock& block=irrep_family_blocks[j_index][n][irrep_family_index];
				        // write matrix.  Order is column major (Eigen default)
				        if(lgi::binary_float_precision==4)
				          {
				            Eigen::MatrixXf buffer_matrix=block.transpose().cast<float>();
				            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);
				            
				          }  
				          
				        else if (lgi::binary_float_precision==8)
				          {
				            Eigen::MatrixXd buffer_matrix=block.transpose();
				            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);

				          }
				      }
					}
			}
    out_file.close();    
  }



void ReadIrrepFamilyBlocks(
	std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,//by irrep family, by J, by n
	std::map<int,std::map<int,int>>& J_index_lookup_table,
  const std::string& filename
  )
{

  std::ios_base::openmode mode_argument = std::ios_base::in | std::ios_base::binary;
  std::ifstream in_stream;
  in_stream.open(filename,mode_argument);

  if(not bool(in_stream))
      std::cout<<filename+" not found."<<std::endl;

  int num_J_values,num_irrep_families,num_eigenvalues, binary_float_precision;
  mcutils::ReadBinary<int>(in_stream,binary_float_precision);  
  mcutils::ReadBinary<int>(in_stream,num_J_values);  
  mcutils::ReadBinary<int>(in_stream,num_eigenvalues);
  mcutils::ReadBinary<int>(in_stream,num_irrep_families);
	
	//for each irrep family
	for(int i=0; i<num_irrep_families; ++i)
		{
			// std::cout<<"Read irrep_family index (corresponds to full space)"<<std::endl;
			int irrep_family_index;
  		mcutils::ReadBinary<int>(in_stream,irrep_family_index);
		
			std::vector<spncci::OperatorBlocks>& blocks=irrep_family_blocks[irrep_family_index];
			blocks.resize(num_J_values);
		  
		  for( int j_index=0; j_index<num_J_values; ++j_index)
		  	{
					int twice_J, rows, cols;
					mcutils::ReadBinary<int>(in_stream,twice_J);
				  mcutils::ReadBinary<int>(in_stream,rows);
					mcutils::ReadBinary<int>(in_stream,cols);

					if(rows<2)
					{
						irrep_family_blocks.erase(irrep_family_index);
						continue;
					}
					
					blocks[j_index].resize(num_eigenvalues);
					J_index_lookup_table[irrep_family_index][twice_J]=j_index;

	      	for(int n=0; n<num_eigenvalues; ++n)
	      		{
			        spncci::OperatorBlock& block=blocks[j_index][n];
			        // std::cout<<"Read matrix.  Order is column major (Eigen default)"<<std::endl;
				      if(binary_float_precision==4)
				        {
				          float buffer[rows*cols];
				          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
				          block=Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
				        }
				      else if (binary_float_precision==8)
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

void  RegroupBlocks(
  const std::vector<std::pair<int,int>>& Jn_set,
  const std::vector<spncci::OperatorBlocks>& blocks,
  std::map<int,int>& J_index_table,
  spncci::OperatorBlock& irrep_family_block
)
//
//Regroup different Jn blocks for a given irrep family into a single block 
//which is gamma_max x sum(Jn_subspaces)
{
  // std::cout<<"Regrouping irrep family blocks "<<blocks.size()<<std::endl;
	int num_rows=blocks[0][0].rows();

	//Counting pass to get num columns of composite Jn block
	int num_cols=0;
	for(auto Jn : Jn_set)
		{
			int twice_J,n;
			std::tie(twice_J,n)=Jn;
			int j_index=J_index_table[twice_J];
			num_cols+=blocks[j_index][n].cols();
		}

	// std::cout<<"Create superblock containing all Jn blocks"<<std::endl;
	irrep_family_block=Eigen::MatrixXd::Zero(num_rows,num_cols);
	int offset=0;
	num_cols=0;
	for(auto Jn : Jn_set)
		{
			int twice_J,n;
			std::tie(twice_J,n)=Jn;
			int j_index=J_index_table[twice_J];
			int num_cols=blocks[j_index][n].cols();
			irrep_family_block.block(0,offset,num_rows,num_cols)=blocks[j_index][n];
			offset+=num_cols;
		}
}


void GetUnitaryTransformation(
  spncci::OperatorBlock& block, 
  spncci::OperatorBlock& transformation_matrix,
  const double zero_threshold
  )
{
  double norm=block.squaredNorm();
  int gamma_max=block.rows();
  int num_cols=block.cols();
  // Maximum probablity of a single irrep
  double max_probability=block.rowwise().squaredNorm().maxCoeff();
  std::cout<<"gamma_max:  "<<gamma_max<<"  norm:  "<<norm<<std::endl;
  std::cout<<"initial max probability "<<max_probability<<std::endl;
  std::cout<<block.rowwise().squaredNorm()<<std::endl;    
  //SVD decomposition

  // Eigen::JacobiSVD<spncci::OperatorBlock> svd(block,ComputeThinV);
  Eigen::JacobiSVD<spncci::OperatorBlock> svd(block,Eigen::ComputeFullU);
  // svd.setThreshold(Eigen::Default);
  // std::cout<<std::endl<<svd.singularValues()<<std::endl<<std::endl;;
  svd.setThreshold(zero_threshold);
  std::cout<< std::endl<<svd.singularValues() << std::endl;
  transformation_matrix=svd.matrixU().transpose();
  int rank=svd.rank();
  std::cout<<"rank "<<rank<<std::endl;

  // Temporary
  spncci::OperatorBlock transformed_block=transformation_matrix*block;

	max_probability=transformed_block.rowwise().squaredNorm().maxCoeff();
	std::cout<<"final max probability "<<max_probability<<std::endl<<std::endl;
	std::cout<<transformed_block.rowwise().squaredNorm()<<std::endl;
	std::cout<<"---------------------------"<<std::endl<<std::endl;
}

void  DefineIrrepFamilyTransformations(
  const std::vector<std::pair<int,int>>& Jn_set,
  const std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  std::map<int,std::map<int,int>>& J_index_lookup_table,
  std::map<int,spncci::OperatorBlock>& transformations,
  double zero_threshold
)
//Set of transformations for a given set of Jn pairs
{
  std::cout<<"defining rotation"<<std::endl;
  for(auto it=irrep_family_blocks.begin(); it!=irrep_family_blocks.end(); ++it)
  	{
  		int irrep_family_index=it->first;
  		const std::vector<spncci::OperatorBlocks>& blocks=it->second;
			spncci::OperatorBlock block;
			std::map<int,int>& J_index_table=J_index_lookup_table[irrep_family_index];
			spncci::RegroupBlocks(Jn_set, blocks,J_index_table, block);
  		
  		int gamma_max=block.rows();
      int irrep_dim=block.cols();
          
      if(gamma_max<2)
        continue;

			spncci::OperatorBlock& transformation_matrix=transformations[irrep_family_index];
			spncci::GetUnitaryTransformation(block, transformation_matrix, zero_threshold);
  	}

}

void WriteTransformationMatrices(  
	const std::map<int,spncci::OperatorBlock>& transformations,
  const std::string& filename
)
{
  std::ios_base::openmode mode_argument = std::ios_base::out | std::ios_base::binary;
  std::ofstream out_file;
  out_file.open(filename,mode_argument);

  if (!out_file)
    {
      std::cerr << "Could not open file '" << filename << "'!" << std::endl;
      return;
    }
	mcutils::WriteBinary<int>(out_file,lgi::binary_float_precision);
  for(auto it=transformations.begin(); it!=transformations.end(); ++it)
  	{
  		int irrep_family_index=it->first;

			const spncci::OperatorBlock& transformation_matrix=it->second;

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


void ReadTransformationMatrices(  
	const std::string& filename,
	std::map<int,spncci::OperatorBlock>& transformations
)
{
  std::ios_base::openmode mode_argument = std::ios_base::in | std::ios_base::binary;
  std::ifstream in_stream;
  in_stream.open(filename,mode_argument);

  if (!in_stream)
    {
      std::cerr << "Could not open file '" << filename << "'!" << std::endl;
      return;
    }

  // Check floating precision is correct
  int binary_float_precision;
  mcutils::ReadBinary<int>(in_stream,binary_float_precision);
  assert(binary_float_precision==lgi::binary_float_precision);

  int irrep_family_index, gamma_max;
  while(!in_stream.eof())
   	{
   		mcutils::ReadBinary<int>(in_stream,irrep_family_index);
   		mcutils::ReadBinary<int>(in_stream,gamma_max);

      // Read matrix.  Order is column major (Eigen default)
      if(lgi::binary_float_precision==4)
        {
          float buffer[gamma_max*gamma_max];
          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
          transformations[irrep_family_index]
          		=Eigen::Map<Eigen::MatrixXf>(buffer,gamma_max,gamma_max).cast<double>();
        }
      else if (lgi::binary_float_precision==8)
        {
          double buffer[gamma_max*gamma_max];
          in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
          transformations[irrep_family_index]
          	=Eigen::Map<Eigen::MatrixXd>(buffer,gamma_max,gamma_max);
        }
   	}
}

}