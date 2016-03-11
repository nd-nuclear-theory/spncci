/****************************************************************
  sp_basis_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  3/10/16 (aem,mac): Created.

****************************************************************/

#include "spncci/sp_basis.h"


int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // construction and string conversion
  ////////////////////////////////////////////////////////////////

  HalfInt Nsigma_0 = HalfInt(11,1);
  std::string filename = "libraries/spncci/lgi-3-3-2-fql-mini.dat";
  
  std::vector<spncci::LGI> basis;
  spncci::GenerateLGIVector(basis,filename,Nsigma_0);

  for (int i=0; i<basis.size(); ++i)
    std::cout << i << " " << basis[i].Str() << std::endl;

} //main
