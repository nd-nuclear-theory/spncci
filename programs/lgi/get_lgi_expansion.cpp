#include <fstream>  
#include "cppformat/format.h"
#include "sp3rlib/u3coef.h"
#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"


int main(int argc, char **argv)
{
  // Basis parameters
  int Nmax=2;
  int N=3;
  int Z=3;
  std::array<int,2> nuclide;
  nuclide[0]=Z;
  nuclide[1]=N; // proton and neutron numbers
  bool intrinsic=true;
  HalfInt Nsigma0 = lgi::Nsigma0ForNuclide(nuclide,intrinsic);
  std::string lsu3shell_basis_filename="lsu3shell_basis.dat"; // Will need to include path to file

  // Operator parameters
  std::string Brel_filename;
  std::string Nrel_filename;

  // Output parameters
  std::string lgi_filename;
  std::string expansion_filename;
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
  
  lgi::GetLGIExpansion(
      lsu3shell_space,lsu3shell_basis_table,
      Brel_filename,Nrel_filename,Z+N, Nsigma0,
      lgi_families, lgi_expansions
    );

  ////////////////////////////////////////////////////////////////
  // writing to file
  ////////////////////////////////////////////////////////////////

  lgi::WriteLGILabels(lgi_families,lgi_filename);


  lgi::WriteLGIExpansion(Z, N, Nmax,lgi_families,lgi_expansions,expansion_filename);



}
