/****************************************************************

****************************************************************/

#include <iostream>

#include "am/am.h"
#include "am/halfint.h"
#include "am/wigner_gsl.h"

#include "fmt/format.h"
#include "spncci/computation_control.h"

int main(int argc, char **argv)
{

  // process arguments
  if(argc<1+1)
    {
      std::cout << "Syntax: J_max" << std::endl;
      std::exit(EXIT_SUCCESS);
    }

// inputs must include Nmax, lgi filename, J_values
  ///////////////////////////////////////////////////////////////////////////////////////
  //// set up SpNCCI spaces
  ///////////////////////////////////////////////////////////////////////////////////////
  // Get LGI families
  // Must be LGI families up through Nmax+2
  std::string lgi_filename="lgi_families.dat";
  lgi::MultiplicityTaggedLGIVector lgi_families
  lgi::ReadLGISet(lgi_filename, Nsigma0,lgi_families);

  // Build SpNCCI basis for Nmax
  spncci::NmaxTruncator truncator0(Nsigma0,Nmax);
  spncci::SpNCCISpace spncci_space0;
  spncci::SigmaIrrepMap sigma_irrep_map0; //Container for spncci space irreps.  
  spncci::GenerateSpNCCISpace(lgi_families,truncator0,spncci_space0,sigma_irrep_map0,restrict_sp3r_to_u3_branching=false);
  baby_spncci_space0=spncci::BabySpNCCISpace(spncci_space0);

  spaces_spbasis0.resize(J_values.size());
  for(int j=0; j<J_values.size(); ++j)
    {
      const HalfInt& J=J_values[j];
      spaces_spbasis[j]=spncci::SpaceSpBasis(baby_spncci_space,J);
    }

  // Build SpNCCI basis for Nmax+2
  spncci::NmaxTruncator truncator2(Nsigma0,Nmax+2);
  spncci::SpNCCISpace spncci_space2;
  spncci::SigmaIrrepMap sigma_irrep_map2; //Container for spncci space irreps.  
  spncci::GenerateSpNCCISpace(lgi_families,truncator2,spncci_space2,sigma_irrep_map2,restrict_sp3r_to_u3_branching=false);
  baby_spncci_space2=spncci::BabySpNCCISpace(spncci_space2);

  spaces_spbasis2.resize(run_parameters.J_values.size());
  for(int j=0; j<run_parameters.J_values.size(); ++j)
    {
      const HalfInt& J=run_parameters.J_values[j];
      spaces_spbasis[j]=spncci::SpaceSpBasis(baby_spncci_space,J);
    }

// Read in eigenvector for Nmax and redistrubte in Nmax+2 container 
// Read eigenvector (may need loop over desired J values )
  std::string eigenvector_filename=...;
  spncci::OperatorBlock eigenvectors;
  spncci::ReadEigenvectors(filename,eigenvectors);

// spaces_spbasis for Nmax and Nmax+2
// Q matrix for Nmax+2
// Eigenvector for Nmax


}
