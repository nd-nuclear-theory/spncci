/****************************************************************
  lsu3shell_basis_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  7/5/16 (aem,mac): Created.
****************************************************************/

#include <fstream>
#include <ostream>  

#include "fmt/format.h"
#include "sp3rlib/u3coef.h"
#include "moshinsky/moshinsky_xform.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/unit_tensor_expansion.h"

#include "lsu3shell/lsu3shell_basis.h"

#include "SU3ME/proton_neutron_ncsmSU3Basis.h"
#include "LSU3/ncsmSU3xSU2Basis.h"


void test_lsu3shell_basis_from_ncsmSU3xSU2BasisLSU3Tablar()
  {

    // setup for test case
    int Nsigma_0=11;  // 11 for 6Li

    // reading in basis table obtained using ncsmSU3xSU2BasisLSU3Tabular
    // std::string lsu3_filename("LSU3ShellBasis_table.dat");
    std::string lsu3_filename("basis_table.dat");

    lsu3shell::LSU3ShellBasisTable basis_table;
    lsu3shell::U3SPNBasisLSU3Labels basis_provenance;
    u3shell::SpaceU3SPN space;
    lsu3shell::ReadLSU3ShellBasis(Nsigma_0,lsu3_filename, basis_table, basis_provenance, space);

    // dump lsu3shell basis information
    for(const lsu3shell::LSU3ShellBasisGroupData& group : basis_table)
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
            const lsu3shell::LSU3ShellBasisGroupLabels& basis_group_labels = basis_provenance[subspace_index][state_index];
            std::cout
              << fmt::format("  state {} Np {} Nn {}",state_index,basis_group_labels.Np,basis_group_labels.Nn)
              << std::endl;
          }

      }
  }


void test_lsu3shell_basis_generated_within_spncci()
  {
    nuclide::NuclideType nuclide={3,3};
    bool intrinsic=false;
    HalfInt Nsigma0=nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
    int Nmax=1;

    std::map<u3shell::U3SPN, unsigned int> u3spn_dimensions
      =lsu3shell::generate_lsu3shell_basis_dimensions(nuclide,Nsigma0,Nmax);

    // std::map<u3shell::U3SPN, lsu3shell::Dimensions> u3spn_dimensions;

    std::map<u3shell::U3SPN, unsigned int> u3spn_cmf_dimensions 
      =lsu3shell::lsu3shell_cmf_basis_dimensions(Nsigma0,Nmax,u3spn_dimensions);

    const std::string input_filename="data/Z03-N03-Nmax14_u3s-dim.dat";
    std::cout<<input_filename<<std::endl;
    std::map<u3shell::U3SPN, lsu3shell::Dimensions> u3spn_dimensions_test;
    lsu3shell::read_lsu3shell_basis_dimensions(
        input_filename,Nsigma0,Nmax,u3spn_dimensions_test
      );



    for(const auto& pair : u3spn_dimensions)
      {
        if(pair.second !=u3spn_dimensions_test[pair.first].total)
          {
            std::cout<<pair.first.Str()<<"  "<<pair.second
            <<"  "<<u3spn_dimensions_test[pair.first].total<<std::endl;
            exit(EXIT_FAILURE);
          }
      }

    for(const auto& pair : u3spn_dimensions)
      {
        std::cout<<fmt::format("{}:  {:6d}  {:6d}",
          pair.first.Str(),pair.second,u3spn_cmf_dimensions[pair.first]
        )<<std::endl;
      }


  

  }



int main(int argc, char **argv)
{
  u3::U3CoefInit();
  // test_lsu3shell_basis_from_ncsmSU3xSU2BasisLSU3Tablar();
  test_lsu3shell_basis_generated_within_spncci();




}//end main
