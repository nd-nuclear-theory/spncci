/****************************************************************
  transform_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

****************************************************************/
#include "spncci/transform_basis.h"

#include <fstream>
#include <omp.h>  

#include "lgi/lgi_unit_tensors.h"
#include "mcutils/eigen.h"
#include "mcutils/io.h"
#include "spncci/parameters.h"
#include "spncci/spncci_common.h"


namespace spncci
{
// void RegroupIntoIrrepFamilies(
//     const spncci::SpaceSpJ& spj_space,
//     int num_irrep_families,
//     int num_eigenvalues,
//     const std::vector<spncci::Matrix>& eigenvectors,
//     std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks
//   )
// 	// irrep_family_blocks: by J, by n, by irrep family
//   {
//     // std::cout<<"regrouping"<<std::endl;
//     //By J, by irrep_family_inde: dimensions of subspace (gamma_max, num_states)
//     std::vector<std::vector<std::pair<int,int>>> irrep_family_subspaces(spj_space.size());

//     // std::cout<<"geting dimensions"<<std::endl;
//     // iterate over j spaces and count number of states in each subspace and get gamma_max
//     for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
//       // for each J subspace
//       {
//         const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
//         auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
       	
//        	//resize container
//         irrep_family_subspacesJ.resize(num_irrep_families);
        
//         for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
//           {
//             // retrieve basis state information
//             StateSpJ spj_state(spj_subspace,spj_state_index);
//             int degeneracy = spj_state.degeneracy();
//             int irrep_family_index = spj_state.irrep_family_index();
//             int gamma_max=spj_state.gamma_max();
//             int upsilon_max=int(degeneracy/gamma_max);
//             auto& dimensions=irrep_family_subspacesJ[irrep_family_index];

//             // increment count of number of states by upsilon max 
//             irrep_family_subspacesJ[irrep_family_index]
//               =std::pair<int,int>(gamma_max,dimensions.second+upsilon_max);
//           }
//       }

//     // std::cout<<"initialize irrep family subspace blocks "<<std::endl;
//     irrep_family_blocks.resize(spj_space.size());
//     for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)    
//       {
//         const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
//         auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
        
//       	// resize by number of eigenvalues at each J
//       	irrep_family_blocks[spj_subspace_index].resize(num_eigenvalues);
      	
//       	// For each degenerate eigenstate of J
//       	for(int n=0; n<num_eigenvalues; ++n)
//       		{
//       			//Resize block
// 						spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_subspace_index][n];
//         		blocks.resize(irrep_family_subspacesJ.size());

//         		//Zero initilize each block if gamma_max > 1
// 		        for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
// 		          {
// 		            int gamma_max,num_states;
// 		            std::tie(gamma_max,num_states)=irrep_family_subspacesJ[irrep_family_index];
// 		            // std::cout<<irrep_family_index<<"  "<<gamma_max<<"  "<<num_states<<std::endl;
		            
//                 //  irrep may not branch to given J, especially for low Nmax
//                 if(gamma_max>0)
//   		  	        blocks[irrep_family_index]=Eigen::MatrixXd::Zero(num_states,gamma_max);
// 		          }  
// 		      }      
//       }

//     // std::cout<<"populating blocks"<<std::endl;
//     for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
//       {
//         auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_subspace_index];
//         const SubspaceSpJ& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
//         const spncci::Matrix& eigenvectors_J = eigenvectors[spj_subspace_index];
//         const int num_eigenvectors = eigenvectors_J.cols();

//         // std::cout<<num_eigenvectors<<" eigenvectors for "<<num_eigenvalues<<" eigenvalues"<<std::endl;
//         // std::cout<<eigenvectors_J<<std::endl;
//         // std::cout<<"-------------"<<std::endl<<std::endl;

//         for(int n=0; n<num_eigenvalues; ++n)
//         	{
// 		        spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_subspace_index][n];
		        
// 		        // initialize offsets to zero for each irrep family
// 		        std::vector<int> offsets(irrep_family_subspacesJ.size(),0);
		        
// 		        for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
// 		          {
// 		            // retrieve basis state information
// 		            StateSpJ spj_state(spj_subspace,spj_state_index);
// 		            int gamma_max=spj_state.gamma_max();

                
//                 // std::cout<<"n here "<<n<<std::endl;
//                 // skip irreps that don't contribute to given J space 
// 		            if(gamma_max==0)
// 		            	continue;

// 		            int degeneracy = spj_state.degeneracy();
// 		            int irrep_family_index = spj_state.irrep_family_index();
		            
// 		            int upsilon_max=degeneracy/gamma_max;
// 		            int eigen_offset=spj_state.offset(); //Starting position in eigenvector

// 		            int offset=offsets[irrep_family_index];
		            
//                 for(int gamma=1; gamma<=gamma_max; ++gamma)
// 		              {
// 		                //Taking the nth eigenvector for the given J value
//                     // std::cout<<"n "<<n<<std::endl;
// 		                // std::cout<<"gamma "<<gamma<<" of "<<gamma_max<<std::endl;
//                     // std::cout<<"index"<<spj_state_index<<" of "<<spj_subspace.size()<<std::endl;
//                     blocks[irrep_family_index].block(offset,gamma-1,upsilon_max,1)
// 		                  =eigenvectors_J.block(eigen_offset,n,upsilon_max,1);

// 		                // Increment offset in eigenvector
// 		                eigen_offset+=upsilon_max;
// 		              }
//                 // std::cout<<blocks[irrep_family_index]<<std::endl<<std::endl;;
// 		            // Increment offset in irrep family block 
// 		            offsets[irrep_family_index]+=upsilon_max;
		            
// 		           //  if(irrep_family_index==6)
// 		           //  {
// 		           //  std::cout<<blocks[irrep_family_index]<<std::endl<<std::endl;
// 		           //  std::cout<<"--------------------------------"<<std::endl<<std::endl;
// 		          	// }
// 		          }
//             // std::cout<<"--------------------------------"<<std::endl<<std::endl;
// 		      }
//       }
//   }

void RegroupIntoIrrepFamilies(
    const std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
    // const spncci::SpaceSpJ& spj_space,
    int num_irrep_families,
    int num_eigenvalues,
    const std::vector<spncci::Matrix>& eigenvectors,
    std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks
  )
  // irrep_family_blocks: by J, by n, by irrep family
  {
    // std::cout<<"regrouping"<<std::endl;
    //By J, by irrep_family_index: dimensions of subspace (gamma_max, num_states)
    std::vector<std::vector<std::pair<int,int>>> irrep_family_subspaces(spaces_spbasis.size());

    // std::cout<<"geting dimensions"<<std::endl;
    // iterate over j spaces and count number of states in each subspace and get gamma_max
    for(int spj_space_index=0; spj_space_index<spaces_spbasis.size(); ++spj_space_index)
      {
        const spncci::SpaceSpBasis& spj_space=spaces_spbasis[spj_space_index];
        // subspaces for a given J?
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_space_index];
        irrep_family_subspacesJ.resize(num_irrep_families);
        
        for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
          {
            const SubspaceSpBasis& spj_subspace = spj_space.GetSubspace(spj_subspace_index);           
            int gamma_max=spj_subspace.gamma_max();
            int irrep_family_index = spj_subspace.irrep_family_index();

            for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
              {
                // retrieve basis state information
                StateSpBasis spj_state(spj_subspace,spj_state_index);
                int kappa_max=spj_state.kappa_max();
                int upsilon_max=spj_state.upsilon_max();
                // int degeneracy = spj_state.degeneracy();
                
                int irrep_degeneracy=upsilon_max*kappa_max;

                auto& dimensions=irrep_family_subspacesJ[irrep_family_index];

                // increment count of number of states by upsilon max 
                irrep_family_subspacesJ[irrep_family_index]
                  =std::pair<int,int>(gamma_max,dimensions.second+irrep_degeneracy);
              }
          }
      }

    // std::cout<<"initialize irrep family subspace blocks "<<std::endl;
    irrep_family_blocks.resize(spaces_spbasis.size());
    
    for(int spj_space_index=0; spj_space_index<spaces_spbasis.size(); ++spj_space_index)
      {
        // std::cout<<"for spj space "<<spj_space_index<<std::endl;
        // const spncci::SpaceSpBasis& spj_space=spaces_spbasis[spj_space_index];
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_space_index];
        irrep_family_blocks[spj_space_index].resize(num_eigenvalues);

        for(int n=0; n<num_eigenvalues; ++n)
          {
            // std::cout<<"for eigenvalue "<<n<<std::endl;
            //Resize block[J]
            spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_space_index][n];
            blocks.resize(num_irrep_families);

            //Zero initilize each block if gamma_max > 1
            for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
              {
                int gamma_max,num_states;
                std::tie(gamma_max,num_states)=irrep_family_subspacesJ[irrep_family_index];
                
                // std::cout<<"for irrep family "<<irrep_family_index<<"  "<<gamma_max<<"  "<<num_states<<std::endl;
                //  irrep may not branch to given J, especially for low Nmax
                if(gamma_max>0)
                  blocks[irrep_family_index]=Eigen::MatrixXd::Zero(num_states,gamma_max);
              }        
          }
      }
    
    // std::cout<<"populating blocks"<<std::endl;
    for(int spj_space_index=0; spj_space_index<spaces_spbasis.size(); ++spj_space_index)
      {
        const spncci::SpaceSpBasis& spj_space=spaces_spbasis[spj_space_index];
        auto& irrep_family_subspacesJ=irrep_family_subspaces[spj_space_index];
        const spncci::Matrix& eigenvectors_J = eigenvectors[spj_space_index];
        const int num_eigenvectors = eigenvectors_J.cols();        
        


        std::vector<std::vector<int>> eigen_offsets;
        spncci::GetSpBasisOffsets(spj_space,eigen_offsets);

        // std::cout<<"-----------------"<<std::endl;
        // std::cout<<eigenvectors_J<<std::endl<<std::endl;

        for(int n=0; n<num_eigenvalues; ++n)
          {
            
            // In some low Nex cases, there will be fewer eigenvalues the num_eigenvalues.  In this case,
            // skip and regrouped block is be zero padded.
            if(n>=num_eigenvectors)
              continue;

            // initialize offsets to zero for each irrep family
            std::vector<int> offsets(irrep_family_subspacesJ.size(),0);

            // std::cout<<"for the "<<n<<"th eigenvalue "<<std::endl;
            spncci::OperatorBlocks& blocks=irrep_family_blocks[spj_space_index][n];

            for (int spj_subspace_index=0; spj_subspace_index<spj_space.size(); ++spj_subspace_index)
              {     
                const SubspaceSpBasis& spj_subspace = spj_space.GetSubspace(spj_subspace_index);
                int gamma_max=spj_subspace.gamma_max();
                int irrep_family_index = spj_subspace.irrep_family_index();

                // std::cout<<"  irrep family "<<irrep_family_index<<"  "<<gamma_max<<std::endl;

                // skip irreps that don't contribute to given J space 
                if(gamma_max==0)
                  continue;
                
                for (int spj_state_index=0; spj_state_index<spj_subspace.size(); ++spj_state_index)
                  {
                    // retrieve basis state information
                    StateSpBasis spj_state(spj_subspace,spj_state_index);
                    int kappa_max = spj_state.kappa_max();
                    int upsilon_max=spj_state.upsilon_max();
                    
                    int eigen_offset=eigen_offsets[spj_subspace_index][spj_state_index];
                    
                    // For kappa, for gamma, insert into matrix 
                    for(int kappa=1; kappa<=kappa_max; ++kappa)
                    {
                      // std::cout<<"offsets "<<spj_state_index<<"  "<<kappa<<"  "<<offsets[irrep_family_index]<<std::endl;
                      //Add in cols corresponding to gamma
                      for(int gamma=1; gamma<=gamma_max; ++gamma)      
                        {
                          int offset=offsets[irrep_family_index];
                          // Get upsilon_max sub-vector from eigenvector for given n and 
                          // add it to block for given irrep family
                          // std::cout<<offset<<"  "<<gamma-1<<"  "<<upsilon_max<<"  "<<blocks[irrep_family_index].rows()<<"  "
                          // <<blocks[irrep_family_index].cols()<<std::endl;
                          // std::cout<<eigenvectors_J.rows()<<"  "<<eigenvectors_J.cols()<<"  "<<eigen_offset<<"  "<<n<<std::endl;
                           blocks[irrep_family_index].block(offset,gamma-1,upsilon_max,1)
                            =eigenvectors_J.block(eigen_offset,n,upsilon_max,1);

                          // Increment offset in eigenvector
                          eigen_offset+=upsilon_max;
                        }
  
                      // Increment offset for next upsilon_max sub-vector
                      offsets[irrep_family_index]+=upsilon_max;
                    }


                  }// end spj_state_index 
                // std::cout<<" block "<<n<<"  "<<irrep_family_index<<std::endl;
                // std::cout<<blocks[irrep_family_index]<<std::endl<<std::endl;   
              } //end spj_subspace_index
          
          } //end n
      } //end spj_space_index
  }



  // void WriteIrrepFamilyBlocks(  
  //   const spncci::SpaceSpJ& spj_space,
  //   int num_irrep_families,
  //   int num_eigenvalues,
  //   const std::vector<int>& lgi_full_space_index_lookup,
  //   const std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  //   const std::string& filename
  // )
  // {
  //   std::ios_base::openmode mode_argument = std::ios_base::out | std::ios_base::binary;
  //   std::ofstream out_file;
  //   out_file.open(filename,mode_argument);

  //   std::cout<<"writing irrep family blocks to file"<<std::endl;

  //   if (!out_file)
  //     {
  //       std::cerr << "Could not open file '" << filename << "'!" << std::endl;
  //       return;
  //     }

  //   // Number of J values, num of eigenstates with given J eigenvalue 
  //   // and num irrep families in full space
  //   mcutils::WriteBinary<int>(out_file,lgi::binary_float_precision);
  //   mcutils::WriteBinary<int>(out_file,spj_space.size());
  //   mcutils::WriteBinary<int>(out_file,num_eigenvalues);
  //   mcutils::WriteBinary<int>(out_file,num_irrep_families);

  //   //for each irrep family
		// for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
		// 	{
		//     // Write irrep family index in full space 
  //   		int full_space_irrep_family_index=lgi_full_space_index_lookup[irrep_family_index];
  //   		mcutils::WriteBinary<int>(out_file,full_space_irrep_family_index);
				
		//     //for reach J eigenvalue
		//     for( int j_index=0; j_index<spj_space.size(); ++j_index)
		// 	    {
		// 				// get J
		// 				HalfInt J = spj_space.GetSubspace(j_index).J();
		      	
		//       	// get number of rows and columns
		// 				const spncci::OperatorBlock& block=irrep_family_blocks[j_index][0][irrep_family_index];
		// 				int rows=block.rows();
		// 				int cols=block.cols();
  //           // std::cout<<"rows "<<rows<<" cols "<<cols<<std::endl;
  //           // std::cout<<block<<std::endl<<std::endl;

		// 				// Write block information, col and rows are transposed
		// 				mcutils::WriteBinary<int>(out_file,TwiceValue(J));
		//     		mcutils::WriteBinary<int>(out_file,cols);
		//     		mcutils::WriteBinary<int>(out_file,rows);

		//     // 		// Only write if gamma_max (corresponding to rows) is >1
		// 				// if(cols==0)
		// 				// 	continue;

		//     		int size=rows*cols;

		//       	// For each n value
		//       	// std::cout<<"Transpose matrix and write to file"<<std::endl;
		//       	for(int n=0; n<num_eigenvalues; ++n)
		//       		{
		// 		        const spncci::OperatorBlock& block=irrep_family_blocks[j_index][n][irrep_family_index];

  //               // if(n==0 && j_index==0)
  //               // {
  //               //   // double max_probability=block.rowwise().squaredNorm().maxCoeff();
  //               //   // std::cout<<"final max probability "<<max_probability<<std::endl<<std::endl;
  //               //   // std::cout<<block.rowwise().squaredNorm()<<std::endl;
  //               //   // std::cout<<"---------------------------"<<std::endl<<std::endl;
  //               // }



		// 		        // write matrix.  Order is column major (Eigen default)
		// 		        if(lgi::binary_float_precision==4)
		// 		          {
		// 		            Eigen::MatrixXf buffer_matrix=block.transpose().cast<float>();
		// 		            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);
				            
		// 		          }  
				          
		// 		        else if (lgi::binary_float_precision==8)
		// 		          {
		// 		            Eigen::MatrixXd buffer_matrix=block.transpose();
		// 		            out_file.write(reinterpret_cast<char*>(buffer_matrix.data()),size*lgi::binary_float_precision);

		// 		          }
		// 		      }
		// 			}
		// 	}
  //   out_file.close();    
  // }


    void WriteIrrepFamilyBlocks(
    std::vector<HalfInt> J_values,  
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

    std::cout<<"writing irrep family blocks to file"<<std::endl;

    if (!out_file)
      {
        std::cerr << "Could not open file '" << filename << "'!" << std::endl;
        return;
      }

    // Number of J values, num of eigenstates with given J eigenvalue 
    // and num irrep families in full space
    mcutils::WriteBinary<int>(out_file,lgi::binary_float_precision);
    mcutils::WriteBinary<int>(out_file,J_values.size());
    mcutils::WriteBinary<int>(out_file,num_eigenvalues);
    mcutils::WriteBinary<int>(out_file,num_irrep_families);

    //for each irrep family
    for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
      {
        // Write irrep family index in full space 
        int full_space_irrep_family_index=lgi_full_space_index_lookup[irrep_family_index];
        mcutils::WriteBinary<int>(out_file,full_space_irrep_family_index);
        
        //for reach J eigenvalue
        for( int j_index=0; j_index<J_values.size(); ++j_index)
          {
            // get J
            HalfInt J = J_values[j_index];
            
            // get number of rows and columns
            const spncci::OperatorBlock& block=irrep_family_blocks[j_index][0][irrep_family_index];
            // mcutils::ChopMatrix(block,1e-10);

            int rows=block.rows();
            int cols=block.cols();

            // Write block information, col and rows are transposed
            mcutils::WriteBinary<int>(out_file,TwiceValue(J));
            mcutils::WriteBinary<int>(out_file,cols);
            mcutils::WriteBinary<int>(out_file,rows);


            int size=rows*cols;
            if(size==0)
              continue;

            // For each n value
            // std::cout<<"Transpose matrix and write to file"<<std::endl;
            for(int n=0; n<num_eigenvalues; ++n)
              {
                const spncci::OperatorBlock& block=irrep_family_blocks[j_index][n][irrep_family_index];
                // std::cout<<"irrep family "<<irrep_family_index<<" n "<<n<<std::endl;
                // std::cout<<block<<std::endl<<std::endl;
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
  // std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,//by irrep family, by J, by n
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
	
  // irrep_family_blocks.resize(num_irrep_families);

	//for each irrep family
	for(int i=0; i<num_irrep_families; ++i)
		{
			// std::cout<<"Read irrep_family index (corresponds to full space)"<<std::endl;
			int irrep_family_index;
  		mcutils::ReadBinary<int>(in_stream,irrep_family_index);
		
      // std::cout<<"irrep family index "<<irrep_family_index<<std::endl;

			std::vector<spncci::OperatorBlocks>& blocks=irrep_family_blocks[irrep_family_index];
			blocks.resize(num_J_values);
		  
		  for( int j_index=0; j_index<num_J_values; ++j_index)
		  	{
					int twice_J, rows, cols;
					mcutils::ReadBinary<int>(in_stream,twice_J);
				  mcutils::ReadBinary<int>(in_stream,rows);
					mcutils::ReadBinary<int>(in_stream,cols);

          // std::cout<<"irrep family index "<<irrep_family_index<<"  "
          // <<j_index<<"  "<<rows<<"  "<<cols<<"  "<<twice_J<<std::endl;

          // irrep did not contribute to J space 
					if(rows==0)
            continue;
					
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
				    // std::cout<<blocks[j_index][n]<<std::endl<<std::endl;
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
  // std::cout<<"Regrouping irrep family blocks "<<blocks[0].size()<<std::endl;
	
	// std::cout<<"Counting pass to get num columns of composite Jn block"<<std::endl;
	int total_num_cols=0;
  int total_num_rows=0;
	for(auto Jn : Jn_set)
		{
			int twice_J,n;
			std::tie(twice_J,n)=Jn;
			int j_index=J_index_table[twice_J];

      //check that irrep family block exists for given Jn
      if(blocks[j_index].size()==0)
        continue;

      // std::cout<<"incrementing total "<<j_index<<"  "<<n<<std::endl;
      // std::cout<<blocks[j_index][n]<<std::endl;
      // std::cout<<blocks.size()<<std::endl;
      // std::cout<<blocks[j_index].size()<<"  "<<n<<"  "<<blocks[j_index][n].cols()<<std::endl;

			total_num_cols+=blocks[j_index][n].cols();
      total_num_rows=std::max(total_num_rows,int(blocks[j_index][n].rows()));
		}

  if(total_num_rows==0)
    return;

	// std::cout<<"Create superblock containing all Jn blocks"<<std::endl;
 //  std::cout<<total_num_cols<<"  "<<total_num_rows<<std::endl;
	irrep_family_block=Eigen::MatrixXd::Zero(total_num_rows,total_num_cols);
	int offset=0;
	int num_cols=0;
	for(auto Jn : Jn_set)
		{
			int twice_J,n;
			std::tie(twice_J,n)=Jn;
			int j_index=J_index_table[twice_J];
			int num_cols=blocks[j_index][n].cols();
      
      // if irrep did not contribute to given Jn space, continue to next irrep
      if(num_cols==0)
        continue;
			
      irrep_family_block.block(0,offset,total_num_rows,num_cols)=blocks[j_index][n];
			offset+=num_cols;
		}
}


void GetUnitaryTransformation(
  spncci::OperatorBlock& block, 
  spncci::OperatorBlock& transformation_matrix,
  const std::pair<std::string,double>& truncation_mode
  )
{
  double norm=block.squaredNorm();
  int gamma_max=block.rows();
  int num_cols=block.cols();
  // Maximum probablity of a single irrep
  std::cout<<block<<std::endl;
  double max_probability=block.rowwise().squaredNorm().maxCoeff();
  std::cout<<"gamma_max:  "<<gamma_max<<"  norm:  "<<norm<<std::endl;
  std::cout<<"initial max probability "<<max_probability<<std::endl;
  std::cout<<block.rowwise().squaredNorm()<<std::endl;    
  //SVD decomposition

  // Eigen::JacobiSVD<spncci::OperatorBlock> svd(block,ComputeThinV);
  Eigen::JacobiSVD<spncci::OperatorBlock> svd(block,Eigen::ComputeFullU);
  svd.setThreshold(truncation_mode.second);
  spncci::OperatorBlock Umatrix=svd.matrixU().transpose();

  // if(Umatrix.determinant()<0)
  //   Umatrix.block(gamma_max-1,0,1,gamma_max)=-1*Umatrix.block(gamma_max-1,0,1,gamma_max);

  int rows;

  // If truncation mode is None, then use default zero threshold and keep
  // all irreps in irrep family (full unitary matrix) 
  if(truncation_mode.first=="None")
    rows=gamma_max;

  else if(truncation_mode.first=="Rank")
    rows=svd.rank();

  else if(truncation_mode.first=="Threshold")
    {
      rows=0;
      // spncci::OperatorBlock temp_transformation=svd.matrixU().transpose();
      spncci::OperatorBlock transformed_block=Umatrix*block;
      const auto& probabilities=transformed_block.rowwise().squaredNorm();

      for(int gamma=1; gamma<=gamma_max; ++gamma)
        {
          double probability=probabilities[gamma-1];
          if(probability>truncation_mode.second)
            rows++;
      }
    }
  else
    std::cout<<"Invalid truncation mode."<<std::endl
      <<"Valid truncation modes are 'None', 'Rank', and 'Threshold'"<<std::endl;

  std::cout<<"rows "<<rows<<std::endl;
  std::cout<< "rank "<<svd.rank()<<std::endl<<svd.singularValues() << std::endl;

  // truncating based on rank of matrix
  if(rows>0)
  {
    // transformation_matrix=Eigen::MatrixXd::Identity(gamma_max,gamma_max);
    std::cout<<"getting transformation matrix "<<std::endl;
    transformation_matrix=Umatrix.transpose().block(0,0,rows,gamma_max);
    std::cout<<"got transformation matrix"<<std::endl;
  }
  
  // Temporary
  if(rows>0)
    {
      spncci::OperatorBlock transformed_block=transformation_matrix*block;
      max_probability=transformed_block.rowwise().squaredNorm().maxCoeff();
      std::cout<<"final max probability "<<max_probability<<std::endl<<std::endl;
      std::cout<<transformed_block.rowwise().squaredNorm()<<std::endl;
      std::cout<<"---------------------------"<<std::endl<<std::endl;
    }
  if(rows==0)
    std::cout<<"no contributions from irrep with block "<<block<<std::endl;

}



// TODO: Switch from vector of blocks to map of blocks. 










void  DefineIrrepFamilyTransformations(
  const std::vector<std::pair<int,int>>& Jn_set,
  std::map<int,std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  // const std::vector<std::vector<spncci::OperatorBlocks>>& irrep_family_blocks,
  std::map<int,std::map<int,int>>& J_index_lookup_table,
  spncci::OperatorBlocks& transformations,
  const std::pair<std::string,double>& truncation_mode
)
//Set of transformations for a given set of Jn pairs
{
  std::cout<<"defining rotation"<<std::endl;
  int num_irrep_families=irrep_family_blocks.size();
  transformations.resize(num_irrep_families);

  for(auto it=irrep_family_blocks.begin(); it!=irrep_family_blocks.end(); ++it)
  // for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
  	{
      int irrep_family_index=it->first;
      std::cout<<"irrep family "<<irrep_family_index<<std::endl;
			const std::vector<spncci::OperatorBlocks>& blocks=it->second;
      // const std::vector<spncci::OperatorBlocks>& blocks=irrep_family_blocks[irrep_family_index];
      std::cout<<irrep_family_index<<"  "<<blocks.size()<<std::endl;
      if(blocks.size()==0)
        continue;
      
      spncci::OperatorBlock block;
			std::map<int,int>& J_index_table=J_index_lookup_table[irrep_family_index];
      
      std::cout<<"regrouping"<<std::endl;
			spncci::RegroupBlocks(Jn_set, blocks,J_index_table, block);
  		
  		int gamma_max=block.rows();
      int irrep_dim=block.cols();
      if(gamma_max==0)
        continue;

      // std::cout<<"get transformation index "<<std::endl; 
			spncci::OperatorBlock& transformation_matrix=transformations[irrep_family_index];

      std::cout<<"transforming"<<std::endl;
      std::cout<<num_irrep_families<<"  "<<irrep_family_index<<std::endl;
  		spncci::GetUnitaryTransformation(block,transformation_matrix,truncation_mode);
  	}

}

void WriteTransformationMatrices(
  const spncci::OperatorBlocks& transformations,  
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
	
  int num_irrep_families=transformations.size();
  mcutils::WriteBinary<int>(out_file,lgi::binary_float_precision);
  mcutils::WriteBinary<int>(out_file,num_irrep_families);

  for(int irrep_family_index=0; irrep_family_index<num_irrep_families; ++irrep_family_index)
  	{
			const spncci::OperatorBlock& transformation_matrix=transformations[irrep_family_index];



      // std::cout<<"for irrep "<<irrep_family_index<<std::endl;
      // std::cout<<transformation_matrix<<std::endl;

			int rows=transformation_matrix.rows();
      int cols=transformation_matrix.cols();

	    // Write irrep family index and rows (should be same as number of columns)
  		mcutils::WriteBinary<int>(out_file,irrep_family_index);
  		mcutils::WriteBinary<int>(out_file,rows);
      mcutils::WriteBinary<int>(out_file,cols);

      if(rows==0)
        continue;

  		int size=rows*cols;

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


void WriteTruncatedLGIs(
    const std::array<int,2>& nuclide,
    // int Nsmax, int Nmax, int truncation_file_num,
    const spncci::OperatorBlocks& transformations,
    const std::string& truncated_lgi_filename
  )
  {
    // std::cout<<"read lgi families"<<std::endl;
    
    // Reading in original lgi families list
    bool intrinsic=true; 
    std::string lgi_filename="lgi_families.dat";
    lgi::MultiplicityTaggedLGIVector lgi_families;
    HalfInt Nsigma0 = lgi::Nsigma0ForNuclide(nuclide,intrinsic);
    lgi::ReadLGISet(lgi_filename, Nsigma0,lgi_families);

    //  Getting new LGI family list after transformation and truncation 
    int num_irrep_families=transformations.size();
    lgi::MultiplicityTaggedLGIVector lgi_families_truncated;
    for(int irrep_family_index=0; irrep_family_index<num_irrep_families; irrep_family_index++)
      {
        const auto& lgi_family=lgi_families[irrep_family_index];
        // std::tie(Nex, sigma,Sp,Sn,S)=lgi_family.irrep.Key();

        int gamma_max=transformations[irrep_family_index].rows();

        if(gamma_max==0)
          continue;
        
        // MultiplicityTagged<lgi::LGI> new_lgi_family(lgi_family.irrep,gamma_max);
        lgi_families_truncated.emplace_back(lgi_family.irrep,gamma_max);
        std::cout<<lgi_family.Str()<<std::endl;
      }

    // write truncated lgi family labels to file  


    lgi::WriteLGILabels(lgi_families_truncated, truncated_lgi_filename);

  }


void ReadTransformationMatrices(  
  	const std::string& filename,
    spncci::OperatorBlocks& transformations
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

    int num_lgi_families;
    mcutils::ReadBinary<int>(in_stream,num_lgi_families);

    std::cout<<"num lgi families "<<num_lgi_families<<std::endl;
    // resize transformation container
    transformations.resize(num_lgi_families);

    int irrep_family_index, rows, cols;
    for(int family_num=0; family_num<num_lgi_families; ++family_num)
    // while(!in_stream.eof())
     	{
     		mcutils::ReadBinary<int>(in_stream,irrep_family_index);
     		mcutils::ReadBinary<int>(in_stream,rows);
        mcutils::ReadBinary<int>(in_stream,cols);

        // std::cout<<irrep_family_index<<"  "<<rows<<"  "<<cols<<std::endl;

        if(rows==0)
          continue;

        // Read matrix.  Order is column major (Eigen default)
        if(lgi::binary_float_precision==4)
          {
            float buffer[rows*cols];
            in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
            transformations[irrep_family_index]
            		=Eigen::Map<Eigen::MatrixXf>(buffer,rows,cols).cast<double>();
          }
        else if (lgi::binary_float_precision==8)
          {
            double buffer[rows*cols];
            in_stream.read(reinterpret_cast<char*>(&buffer),sizeof(buffer));
            transformations[irrep_family_index]
            	=Eigen::Map<Eigen::MatrixXd>(buffer,rows,cols);
          }
        // std::cout<<transformations[irrep_family_index]<<std::endl;
     	}
  }

void TransformSeeds(
  int bra_index,int ket_index,
  spncci::OperatorBlocks& transformations,
  basis::OperatorBlocks<double>& unit_tensor_seed_blocks
  )
  {
    spncci::OperatorBlock& bra_transformation=transformations[bra_index];
    spncci::OperatorBlock& ket_transformation=transformations[ket_index];


    // std::cout<<"----------------------------------------------------------"<<std::endl;
    // std::cout<<bra_index<<"  "<<ket_index<<std::endl;
    // std::cout<<bra_transformation<<std::endl<<std::endl;
    // std::cout<<ket_transformation<<std::endl<<std::endl;
    // // std::cout<<"          ------------------------------          "<<std::endl;

    for(int i=0; i<unit_tensor_seed_blocks.size(); ++i)
      { 
        spncci::OperatorBlock& block=unit_tensor_seed_blocks[i];
        spncci::OperatorBlock block2=unit_tensor_seed_blocks[i];
        // std::cout<<"block "<<i<<std::endl;
        // std::cout<<block<<std::endl<<std::endl;
        block=bra_transformation*block*ket_transformation.transpose();
        // std::cout<<block<<std::endl;
        // std::cout<<"          ------------------------------          "<<std::endl;
        // spncci::OperatorBlock temp=block-block2;
        // mcutils::ChopMatrix(temp, 1e-7);
        // std::cout<<temp<<std::endl;
      }
  }



}//end namespace