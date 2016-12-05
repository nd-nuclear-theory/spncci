/****************************************************************
  compute_unit_tensor_rmes.cpp

  compute relative unit tensor rmes in spncci basis using recurrance
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  12/3/16 (aem): Created.  Based on unit_tensor_test.cpp
****************************************************************/
#include <cstdio>
#include <sys/resource.h>

#include "lgi/lgi.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h" 
#include "spncci/unit_tensor.h"

// TODO::Move out of global space 
spncci::SpIrrepVector sp_irrep_vector;

std::unordered_map<u3::U3,vcs::MatrixCache, boost::hash<u3::U3> >  K_matrix_map;

int main(int argc, char **argv)
{
  u3::U3CoefInit();
  //unit tensor cache 
  u3::UCoefCache u_coef_cache;
	u3::g_u_cache_enabled = true;

  // parse arguments
  if (argc<2)
    {
      std::cout << "Syntax: lgi_file_name Nmax" << std::endl;
      std::exit(1);
    }
  std::string filename(argv[1]);
  int Nmax = std::stoi(argv[2]);

  // Generate vector of LGI's from input file 
  lgi::LGIVector lgi_vector;
  lgi::ReadLGISet(lgi_vector,filename);

  // // For generating the lgi_vector, using Li-6 as example;  
  
  // HalfInt Nsigma_0 = lgi_vector[0].firstHalfInt(11,1);
  // int N1b=2;


}
