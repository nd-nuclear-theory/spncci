/****************************************************************
  get_lgi_expansion.cpp

  Code generating expansion of lgi subspaces in terms of 
  lsu3shell su3 reduced basis states 

  Inputs:
    directory lsu3shell_rme containing files:
      * lsu3shell_basis.dat
      * Brel.rme
      * Nrel.rme

  Output: 
    "lgi_families.dat" : files containing lgi labels 
        Nex  lambda mu 2Sp 2Sn 2S gamma

    "lgi_expansions.dat" : file containing expansion of lgi's in lsu3shell basis
  Anna E. McCoy
  TRIUMF

  3/30/19 (aem): Recreated. 
****************************************************************/
#include <fstream>  
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "sp3rlib/u3coef.h"
#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"
#include "lgi/lgi_unit_tensors.h"
#include "spncci/spncci_basis.h"
#include "u3shell/relative_operator.h"

#include "mcutils/io.h"

// Testing function 
namespace lgi{}//namespace

int main(int argc, char **argv)
{
  if(argc<3)
  {
    std::cout<<"Syntax: Z N Nmax "<<std::endl;
    std::exit(EXIT_FAILURE);
  }

  // nuclide 
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);

  // Basis parameters
  int Nmax=std::stoi(argv[3]);
  

  // zero tolerance 
  lgi::zero_tolerance=1e-6;
  
  // output mode
  lgi::binary_format_code = 1;
  lgi::binary_float_precision=8;

  std::array<int,2> nuclide; // proton and neutron numbers
  nuclide[0]=Z;
  nuclide[1]=N; 

  bool intrinsic=true;

  // su3rme output files
  std::string su3rme_filename_base="lsu3shell_rme";
  std::string lsu3shell_basis_filename=su3rme_filename_base+"/lsu3shell_basis.dat"; // Will need to include path to file

  // Generate Nsigma0 from nuclei and type 
  HalfInt Nsigma0 = lgi::Nsigma0ForNuclide(nuclide,intrinsic);

  // Operator parameters
  std::string Brel_filename=su3rme_filename_base+"/Brel.rme";
  std::string Nrel_filename=su3rme_filename_base+"/Nrel.rme";

  // Unit tensor parameters
  int J0=-1;
  int T0=-1;

  int N1v=spncci::ValenceShellForNuclide(nuclide);

  //TODO
  // std::string relative_unit_tensor_filename_template = su3rme_filename_base + "/relative_unit_{:06d}.rme";

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis
  ////////////////////////////////////////////////////////////////
  std::cout << "Read lsu3shell basis..." << std::endl;
  // read lsu3shell basis (regroup into U3SPN subspaces)
  lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3ShellBasis(
      Nsigma0, lsu3shell_basis_filename,lsu3shell_basis_table,
      lsu3shell_basis_provenance,lsu3shell_space
    );

  ////////////////////////////////////////////////////////////////
  // solve for LGIs
  ////////////////////////////////////////////////////////////////
  std::cout << "Solve for LGIs..." << std::endl;
  lgi::MultiplicityTaggedLGIVector lgi_families;
  lsu3shell::OperatorBlocks lgi_expansions;
  std::vector<int> lsu3hsell_index_lookup_table;

  lgi::GetLGIExpansion(
      lsu3shell_space,lsu3shell_basis_table,
      Brel_filename,Nrel_filename,Z+N, Nsigma0,
      lgi_families, lgi_expansions,
      lsu3hsell_index_lookup_table
    );
  
  std::string lgi_filename="lgi_families.dat";
  lgi::WriteLGILabels(lgi_families, lgi_filename);

  // std::cout<<"write expansion to file "<<std::endl;
  std::string lgi_expansion_filename="lgi_expansions.dat";
  lgi::WriteExpansion(lgi_expansion_filename,lgi_expansions);
    
}
