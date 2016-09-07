/****************************************************************
  lsu3shell_basis_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/5/16 (aem,mac): Created.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "cppformat/format.h"

#include "spncci/sp_basis.h"
#include "spncci/lgi_unit_tensor.h"
#include "u3shell/moshinsky.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/unit_tensor_expansion.h"

#include "lsu3shell/lsu3shell_basis.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  int Nsigma_0=6;

  // Reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  std::string lsu3_filename="lsu3basis_table.dat";
  // lsu3_filename="/Users/annamccoy/projects/spncci/data/lgi_set/lsu3_basis.dat";
  lsu3shell::LSU3BasisTable lsu3basis_table;
  std::map<u3shell::U3SPN,int> subspace_dimensions;

  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, lsu3basis_table, subspace_dimensions);

  int dim, start_index;
  u3shell::U3SPN omegaSPN;
  for(auto lgi_group : lsu3basis_table)
    {
      std::tie(omegaSPN,dim,start_index)=lgi_group;
      std::cout<<omegaSPN.Str()<<"  "<<dim<<"  "<<start_index<<std::endl;
    }
  std::cout<<" "<<std::endl;
  int dim_tot;
  for(auto a : subspace_dimensions)
    {
      omegaSPN=a.first;
      dim_tot=a.second;
      std::cout<<omegaSPN.Str()<<"  "<<dim_tot<<"  "<<std::endl;
    }

  // Constructing the space from subspace_dimensions
  u3shell::SpaceU3SPN space(subspace_dimensions);



}//end main
