/****************************************************************
  get_spncci_dimensions.cpp

  WARNING: Will not build since requires old version of
  lgi::GenerateLGIExpansion with stream arguments.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/25/17 (aem): Created.
****************************************************************/
#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_operator.h"
#include "lgi/lgi.h"
#include "lgi/lgi_solver.h"
#include "sp3rlib/u3coef.h"
#include "u3shell/unit_tensor_expansion.h"
#include "spncci/spncci_basis.h"

int main(int argc, char **argv)
{
u3::U3CoefInit();
  double zero_threshold=1e-6;
  // parameters for null space solver
  double threshold=1e-4;
  bool verbose=false;

  if(argc<8)
    std::cout<<"Syntax : Z N twice_Nsigma_0  Nsigma0_ex, Nmax  N1B  <basis>  <nrel>  <brel>"<<std::endl;
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int A=N+Z;
  int twice_Nsigma_0=std::stoi(argv[3]);
  HalfInt Nsigma_0=HalfInt(twice_Nsigma_0,2);
  int Nsigma0_ex=std::stoi(argv[4]);
  int Nmax=std::stoi(argv[5]);
  int N1B=std::stoi(argv[6]);
  std::string lsu3_filename = argv[7];
  std::string nrel_filename = argv[8];
  std::string brel_filename = argv[9];

  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;

  // Read in the basis
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);
  
  // Solving for lgi expansion from Brel+Ncm
  std::cout<<"Generating LGI expansion"<<std::endl;
  std::ifstream is_nrel(nrel_filename);
  std::ifstream is_brel(brel_filename);

  basis::MatrixVector lgi_expansion_matrix_vector;
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  lgi::GenerateLGIExpansion(A,Nsigma_0+Nsigma0_ex,basis_table,space, is_brel,
  	is_nrel,lgi_vector,lgi_expansion_matrix_vector);
  
  is_nrel.close();
  is_brel.close();

  std::cout<<"Finished Generating expansion "<<std::endl;
  
  // Writing lgi family labels to file
  std::ofstream os("LGI_file.dat");
  lgi::WriteLGILabels(lgi_vector,os);
  os.close();
  
  spncci::SpNCCISpace spncci_space;
  spncci::SigmaIrrepMap sigma_irrep_map;  // persistent container to store branchings
  spncci::NmaxTruncator truncator(Nsigma_0+Nsigma0_ex,Nmax);
  spncci::GenerateSpNCCISpace(lgi_vector,truncator,spncci_space,sigma_irrep_map);

  
  std::ofstream dim_os(fmt::format("dimensions_{:02d}_{:02d}_{:02d}_{:02d}.dat",Z,N,Nsigma0_ex,Nmax));
  //For each N and each J, calculate dimension
  for(HalfInt J=0; J<=4; ++J)
	  {
		  for(N=Nsigma0_ex; N<=Nmax; N+=2)
		  {
			  spncci::SpNCCISpace spncci_space;
			  spncci::SigmaIrrepMap sigma_irrep_map;  // persistent container to store branchings
			  spncci::NmaxTruncator truncator(Nsigma_0,N);
			  spncci::GenerateSpNCCISpace(lgi_vector,truncator,spncci_space,sigma_irrep_map);
		  	int dim=spncci::TotalDimensionU3LSJConstrained(spncci_space,J);
				dim_os << N <<" " <<std::setw(2)<<J<<" "<<std::setw(2)<<dim<<std::endl;
		  }
	  }
}
