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
  // LGI list
  ////////////////////////////////////////////////////////////////

  HalfInt Nsigma_0 = HalfInt(11,1);
  std::string filename = "libraries/spncci/lgi-3-3-2-fql-mini.dat";
  
  spncci::LGIVectorType lgi_vector;
  spncci::GenerateLGIVector(lgi_vector,filename,Nsigma_0);

  for (int i=0; i<lgi_vector.size(); ++i)
    std::cout << i << " " << lgi_vector[i].Str() << std::endl;

  ////////////////////////////////////////////////////////////////
  // build subspaces
  ////////////////////////////////////////////////////////////////

//  SigmaSpaceMapType sigma_space_map;
//  spncci::GenerateSp3RSpaces(lgi_vector,sigma_space_map);
//  for (auto it = sigma_space_map.begin(); it != sigma_space_map.end(); ++it)
//    {
//      const Sp3RSpace& space = (*it);
//      
//    }
//  std::cout << i << " " << lgi_vector[i].Str() << std::endl;





} //main
