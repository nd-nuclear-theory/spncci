/****************************************************************
  sp_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/sp_basis.h"

#include <cstdio>



namespace spncci
{

  void GenerateLGIVector(std::vector<LGI>& basis, const std::string& filename, const HalfInt& Nsigma_0)
  {

    FILE* lgi_file;
    
    lgi_file = std::fopen(filename.c_str(),"r");

    //   Nex 2Sp 2Sn 2S lambda mu count
  }


}  // namespace
