/****************************************************************
  lsu3shel_interface.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/5/16 (aem,mac): Created.
****************************************************************/
#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell_io/lsu3shell_interface.h"
#include "spncci/sp_basis.h"
#include "spncci/lgi_unit_tensor.h"
#include "u3shell/moshinsky.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/unit_tensor_expansion.h"

int main(int argc, char **argv)
{
	// u3::U3CoefInit();

 //  int Nsigma_0=11;
	// std::string lsu3_filename
 //    ="/Users/annamccoy/projects/spncci/data/lgi_set/lsu3_basis.dat";
	// lsu3shell::LSU3Vector lsu3basis_vector;
	// // lsu3shell::ReadLSU3Vector(lsu3_filename, lsu3basis_vector);

	// // for(int i=0; i<lsu3basis_vector.size(); ++i)
	// // 	std::cout<<i<<"  "<<lsu3basis_vector[i].Str()<<std::endl;

	// spncci::LGIVectorType lgi_vector;
	// int Nsigma_begin;

	// std::string lgi_filename
 //    ="/Users/annamccoy/projects/spncci/data/lgi_set/lgi.dat";

	// int Nsigma_min=0, Nsigma_max=0; 
	// std::string brel_filename="None", nrel_filename="None";

	// std::string lgi_expansion_filename
	// 	="/Users/annamccoy/projects/spncci/data/lgi_set/lgi_lsu3.dat";

	// // lsu3shell::GenerateLSU3ShellExpansionLGI(
	// // 		Nsigma_0,Nsigma_min, Nsigma_max,
	// // 		lsu3_filename, brel_filename, nrel_filename,
	// // 		lgi_filename,lgi_expansion_filename
	// 	// );

	// // Eigen::MatrixXd ket;
	// // spncci::ReadLGISU3Expansion(lgi_expansion_filename, ket, "ket");
	// // Eigen::MatrixXd bra;
	// // spncci::ReadLGISU3Expansion(lgi_expansion_filename, bra, "bra");

	// // std::cout<<bra<<std::endl<<std::endl;
	// // std::cout<<ket<<std::endl;

	// int Nmin=0,Nmax=2;
	// std::vector<u3shell::RelativeUnitTensorLabelsU3ST> relative_tensor_labels;

	// u3shell::RelativeUnitTensorCoefficientsU3ST Brel_operator;
	// u3shell::BrelRelativeUnitTensorExpansion(Nmax,Nmax, Brel_operator);
	// // lsu3shell::GenerateLSU3ShellOperators(Nmax, Brel_operator, 10);

 //  u3shell::RelativeUnitTensorCoefficientsU3ST Nrel_operator;
 //  std::string filename="Nrel";
 //  u3shell::NrelRelativeUnitTensorExpansion(Nmin,Nmax, Nrel_operator,8);
 //  lsu3shell::GenerateLSU3ShellOperators(Nmax, Nrel_operator,filename);

  
 //  u3shell::TwoBodySpaceU3ST space(Nmax);
 //  u3shell::TwoBodyUnitTensorCoefficientsU3ST Nrel_twobody;
 //  u3shell::TransformRelativeTensorToTwobodyTensor(Nrel_operator,space,Nrel_twobody, "NAS");
 //  for(auto it=Nrel_twobody.begin(); it!=Nrel_twobody.end(); ++it)
 //  	{
 //  		std::cout<<it->first.Str()<<"  "<<it->second<<std::endl;
 //  	}
  
 //  // std::string filename_identity="identity";
 //  // u3shell::RelativeUnitTensorCoefficientsU3ST Id_operator;
 //  // IdentityRelativeUnitTensorExpansion(0,Nmax, Id_operator);
 //  // lsu3shell::GenerateLSU3ShellOperators(Nmax, Id_operator,filename_identity);

}//end main