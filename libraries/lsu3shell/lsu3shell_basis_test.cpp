/****************************************************************
  lsu3shell_basis_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/5/16 (aem,mac): Created.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "cppformat/format.h"
#include "sp3rlib/u3coef.h"
#include "spncci/sp_basis.h"
// #include "spncci/lgi_unit_tensor.h"
#include "u3shell/moshinsky.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/unit_tensor_expansion.h"

#include "lsu3shell/lsu3shell_basis.h"

int main(int argc, char **argv)
{
  u3::U3CoefInit();

  // setup for test case
  int Nsigma_0=11;  // 11 for 6Li

  // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
  std::string lsu3_filename("lsu3basis_table.dat");
  lsu3shell::LSU3BasisTable basis_table;
  lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
  u3shell::SpaceU3SPN space;
  lsu3shell::ReadLSU3Basis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);

  // dump lsu3shell basis information
  for(const lsu3shell::LSU3BasisGroupData& group : basis_table)
    {
      std::cout
        << fmt::format(
            "{:20} dim {:6} start_index {:6}",
            group.omegaSPN.Str(), group.dim, group.start_index
          )
        <<std::endl;
    }
  std::cout<<" "<<std::endl;

  // dump U3SPN basis subspace info
  std::cout << "space" << " " << space.size() << std::endl;
  for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
    {
      const u3shell::SubspaceU3SPN& subspace = space.GetSubspace(subspace_index);
      std::cout
        << fmt::format("subspace {} labels {} dim {}",
                       subspace_index,
                       subspace.U3SPN().Str(),
                       subspace.size()
          )
        << std::endl;
      std::cout
        << fmt::format("provenance dim {}",basis_provenance[subspace_index].size())
        << std::endl;
      for (int state_index=0; state_index < subspace.size(); ++state_index)
        {
          const lsu3shell::LSU3BasisGroupLabels& basis_group_labels = basis_provenance[subspace_index][state_index];
          std::cout
            << fmt::format("  state {} Np {} Nn {}",state_index,basis_group_labels.Np,basis_group_labels.Nn)
            << std::endl;
        }

    }




}//end main
