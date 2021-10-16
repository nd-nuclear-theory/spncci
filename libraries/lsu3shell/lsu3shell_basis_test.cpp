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


void generate_lsu3shell_basis_dimensions(
  const lgi::NuclideType& nuclide, 
  const HalfInt& Nsigma0,
  const int& Nmax, 
  std::map<u3shell::U3SPN, Dimensions>& u3spn_dimensions
  )
{
  const auto&[Z,N]=nuclide;

  //SU3ME/proton_neutron_ncsmSU3Basis.h
  proton_neutron::ModelSpace lsu3shell_model_space(Z,N,Nmax);
  
  //LSU3/ncsmSU3xSU2Basis.cpp
  int idiag=0; int ndiag=1;
  lsu3::CncsmSU3xSU2Basis lsu3shell_basis(lsu3shell_model_space, idiag, ndiag);

  // Iterate over basis and regroup by Nex,lambda,mu,Sp,Sn,S
  //  loop over (ip, in) pairs
  for (int ipin_block = 0; ipin_block < basis.NumberOfBlocks(); ipin_block++) 
  {
    // If block is empty, continue
    if (!basis.NumberOfStatesInBlock(ipin_block)) {continue;}
    
    int ip = basis.getProtonIrrepId(ipin_block);
    int in = basis.getNeutronIrrepId(ipin_block);
    int Nex = basis.nhw_p(ip) + basis.nhw_n(in);

    int alpha_p_max = basis.getMult_p(ip);
    int alpha_n_max = basis.getMult_n(in);

    HalfInt Sp(irrep_p(basis.getProtonSU3xSU2(ip)).S2,2);
    HalfInt Sn(irrep_n(basis.getNeutronSU3xSU2(in)).S2,2);

    for (int iwpn = basis.blockBegin(ipin_block); iwpn < basis.blockEnd(ipin_block); ++iwpn) 
    {
      omega_pn = basis.getOmega_pn(ip, in, iwpn);
      HalfInt S(omega_pn.S2,2);
      u3::U3 omega(Nex+Nsigma0, {omega_pn.lm,omega_pn.mu})
      int dim=alpha_n_max*alpha_p_max*omega_pn.rho;
      u3shell::U3SPN omegaSPN({omega,S},Sp,Sn);
      u3spn_dimensions[omegaSPN]=Dimensions(dim,dim,dim);
    
    }
  }
}



int main(int argc, char **argv)
{
  u3::U3CoefInit();

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




}//end main
