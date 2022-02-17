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

#include "u3shell/relative_operator.h"

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


  std::string generator_dir = operator_dir+"/relative_generators";
  basis::OperatorBlocks<double> lgi_expansions
  	=lgi::generate_lgi_expansion(nuclide,Nsigma_max,lgi_vector,generator_dir,world_comm);


  int J0=-1;
  int T0=-1;
  std::vector<u3shell::RelativeUnitTensorLabelsU3ST> unit_tensor_labels;
  u3shell::GenerateRelativeUnitTensorLabelsU3ST(Nsigma_max,N1v,unit_tensor_labels,J0,T0,false);
  // 0(0,1)[1/2,1/2,0]
  // int i_bra=4;
  // int i_ket=2;
  for(int i_bra=0; i_bra<lgi_vector.size(); ++i_bra)
    for(int i_ket=0; i_ket<lgi_vector.size(); ++i_ket)
      { 
        // if (i_bra!=4 || i_ket!=2)
        //   continue;
        std::pair<MultiplicityTagged<lgi::LGI>,MultiplicityTagged<lgi::LGI>>
          lgi_pair(lgi_vector[i_bra],lgi_vector[i_ket]);
        
        std::cout<<i_bra<<"  "<<lgi_vector[i_bra].irrep.Str()<<std::endl;
        std::cout<<i_ket<<"  "<<lgi_vector[i_ket].irrep.Str()<<std::endl;
        const basis::OperatorBlock<double>& lgi_expansion_bra = lgi_expansions[i_bra];
        const basis::OperatorBlock<double>& lgi_expansion_ket = lgi_expansions[i_ket];

        std::vector<basis::OperatorBlocks<double>> lgi_unit_tensor_rmes
          =lgi::ComputeSeeds(
            nuclide,Nsigma_max,N1v,
            lgi_pair,
            lgi_expansion_bra,
            lgi_expansion_ket,
            unit_tensor_labels,
            operator_dir,
            world_comm
            );

        for(const auto& blocks : lgi_unit_tensor_rmes)
          for(const basis::OperatorBlock<double>& block : blocks)
            std::cout<<mcutils::FormatMatrix(block,"3.2f")<<std::endl<<std::endl;

      }





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
