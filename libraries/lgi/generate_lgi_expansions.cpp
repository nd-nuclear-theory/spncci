/****************************************************************
  lgi_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/7/16 (aem,mac): Created.
  2/15/18 (aem): Update tests for Nsigma0ForNuclide and ReadLGISet
  10/24/21 (aem): Add test for lgi vector construction using lsu3shell
****************************************************************/
#include "lgi/lgi.h"

#include <cstdlib>

#include "lgi/lgi_gen.h"
#include "lgi/dimensions.h"
#include "utilities/nuclide.h"
#include "mcutils/eigen.h"


// operator_dir = ${SPNCCI_OPERATOR_DIR}/rununittensor01/

int main(int argc, char **argv)
{
	if(argc<4)
	{
		std::cout<<"Syntax: Z N Nsigma_max operator_dir"<<std::endl;
	}
  MPI_Init(&argc, &argv);
  int Z = std::stoi(argv[1]);
  int N = std::stoi(argv[2]);
  int Nsigma_max = std::stoi(argv[3]);
 	std::string operator_dir = argv[4];

  nuclide::NuclideType nuclide({Z,N});
  bool intrinsic = true;
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
  unsigned int N0 = nuclide::N0ForNuclide(nuclide);
  int N1v = nuclide::ValenceShellForNuclide(nuclide);
  
  //Generate LGI vector by finding possible cmf LGI by counting arguments 
  lgi::MultiplicityTaggedLGIVector lgi_vector = lgi::get_lgi_vector(nuclide,Nsigma0,Nsigma_max);
  
  if(true)  
    {
      for(const auto& lgi : lgi_vector)
        std::cout<<lgi.Str()<<std::endl;
    }

	MPI_Comm world_comm=MPI_COMM_WORLD;


  basis::OperatorBlocks<double> lgi_expansions
  	=lgi::generate_lgi_expansion(nuclide,Nsigma_max,lgi_vector,operator_dir,world_comm);

 //  std::cout<<"printing everything "<<std::endl;
	// for(int i=0; i<lgi_vector.size(); ++i)
	// 	{
	// 		const auto&[lgi,gamma_max] = lgi_vector[i];
	//     //Get labels
	//     const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();
	//     fmt::print("{} {} [{} {} {}]\n",Nex,sigma,Sp,Sn,S);
	//     const basis::OperatorBlock<double>& expansion = lgi_expansions[i];
	//     std::cout<<mcutils::FormatMatrix(expansion,".3f")<<std::endl;
	// 	}

	MPI_Finalize();

}
