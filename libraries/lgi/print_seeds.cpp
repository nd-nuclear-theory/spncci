/****************************************************************
  print_seeds.cpp

  Anna E. McCoy 
  TRIUMF

  SPDX-License-Identifier: MIT

  10/11/19 (aem): Created.
  
  For a given pair of LGI, print out RMEs for each unit tensor 
  between LGI pair
****************************************************************/
#include "lgi/lgi_unit_tensors.h"
#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "u3shell/relative_operator.h"
int main(int argc, char **argv)
{
	if(argc<7)
		{
			std::cout<<"Syntax: Nmax N1v J0 T0 <operators_filename> <rme_filename>"<<std::endl;
			std::cout<<"	Nmax: Basis truncation parameters"<<std::endl;
			std::cout<<" 	N1v: N of valence space"<<std::endl;
			std::cout<<"	J0: operator angular momentum.  J0=-1 or all J0"<<std::endl;
			std::cout<<"	T0: operator isospin.  T0=-1 for all T0"<<std::endl;
			std::cout<<"	<operators_filename>: text file containing list of unit tensor operator labels"<<std::endl;
			std::cout<<"	<rme_filename>: binary file containing RMEs"<<std::endl;
			std::cout<<"		Note: Nmax, N1v, J0, T0 must match values used to generate rmes in rme file"<<std::endl;
			std::exit(EXIT_FAILURE);
		}
	int Nmax=std::stoi(argv[1]);
	int N1v=std::stoi(argv[2]);
	int J0=std::stoi(argv[3]);
	int T0=std::stoi(argv[4]);

	std::string lgi_unit_tensor_filename=argv[5];
	std::string seed_filename=argv[6];

	std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_unit_tensor_labels;
	u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nmax,N1v,relative_unit_tensor_labels,J0,T0,false);


	//Read in unit tensor labels 
	std::vector<u3shell::RelativeUnitTensorLabelsU3ST> lgi_unit_tensors;
  std::vector<int> rho0_values;

  bool files_found=lgi::ReadUnitTensorLabels(lgi_unit_tensor_filename,lgi_unit_tensors,rho0_values);
  if(not files_found)
  	{
  		std::cout<<"file "<<lgi_unit_tensor_filename<<" not found"<<std::endl;
  		std::exit(EXIT_FAILURE);
  	}

  // Reads in unit tensor seed blocks and stores them in a vector of blocks. 
  // Order corresponds to order of (unit_tensor,rho0) pairs in corresponding operator file. 
  basis::OperatorBlocks<double> unit_tensor_seed_blocks;
  files_found=lgi::ReadBlocks(seed_filename, lgi_unit_tensors.size(), unit_tensor_seed_blocks);
  
  if(not files_found)
  	{
  		std::cout<<"file "<<seed_filename<<" not found"<<std::endl;
  		std::exit(EXIT_FAILURE);
  	}

  u3shell::RelativeUnitTensorLabelsU3ST unit_tensor;
  int tensor_index=-1;
  for(int i=0; i<lgi_unit_tensors.size(); ++i)
  	{
  		int rho0=rho0_values[i];
  		const basis::OperatorBlock<double> block=unit_tensor_seed_blocks[i];
  		const u3shell::RelativeUnitTensorLabelsU3ST new_tensor=lgi_unit_tensors[i];

  		// If new tensor print tensor index
  		if(not (unit_tensor==new_tensor))
  			{
  				//Check that new tensor has rho0 starting with 1
  				assert(rho0==1);

  				unit_tensor=new_tensor;
  				// iterate through list of tensors to find tensor index 
  				// When found, print tensor information
  				for(int tensor_index=0; tensor_index<relative_unit_tensor_labels.size(); ++tensor_index)
  					{
  						if(relative_unit_tensor_labels[tensor_index]==unit_tensor)
  							{
  								std::cout<<"Tensor: "<<tensor_index<<std::endl;
  								std::cout<<unit_tensor.Str()<<std::endl;	
  								break;
  							}

  						if(tensor_index==relative_unit_tensor_labels.size())
  							{
  								std::cout<<"tensor not found "<<std::endl;
  								std::exit(EXIT_FAILURE);
  							}

  					}
  			}
  		std::cout<<"	rho0: "<<rho0<<std::endl;
  		std::cout<<mcutils::FormatMatrix(block, "+.6f", "	")<<std::endl<<std::endl;
  	}
}
