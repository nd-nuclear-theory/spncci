/****************************************************************
  lgi_gen.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT
****************************************************************/
#include "u3ncsm/dimensions.h"

#include <fstream>
#include <iostream>

#include "LSU3/ncsmSU3xSU2Basis.h"
#include "SU3ME/proton_neutron_ncsmSU3Basis.h"

#include "sp3rlib/u3boson.h"

namespace lsu3shell
{

  std::map<u3shell::U3SPN,unsigned int>
  LSU3ShellBasisDimensions(
    const nuclide::NuclideType& nuclide,
    const HalfInt& Nsigma0,
    const int Nmax
  )
  // For a given Nucleus, determine number of U(3)SpSnS irreps in an Nmax truncated basis
  {
    const auto&[Z,N]=nuclide;
    std::map<u3shell::U3SPN,unsigned int> basis_dimensions;
    for(int Nex=0; Nex<=Nmax; ++Nex)
    {
      //SU3ME/proton_neutron_ncsmSU3Basis.h
      proton_neutron::ModelSpace lsu3shell_model_space(Z,N,Nex);
      
      //LSU3/ncsmSU3xSU2Basis.cpp
      int idiag=0; int ndiag=1;
      lsu3::CncsmSU3xSU2Basis lsu3shell_basis(lsu3shell_model_space, idiag, ndiag);

      // Iterate over basis and regroup by Nex,lambda,mu,Sp,Sn,S
      //  loop over (ip, in) pairs
      for (int ipin_block = 0; ipin_block < lsu3shell_basis.NumberOfBlocks(); ipin_block++) 
        {
          // If block is empty, continue
          if (!lsu3shell_basis.NumberOfStatesInBlock(ipin_block)) {continue;}
          
          unsigned int ip = lsu3shell_basis.getProtonIrrepId(ipin_block);
          unsigned int in = lsu3shell_basis.getNeutronIrrepId(ipin_block);
          unsigned int Nex = lsu3shell_basis.nhw_p(ip) + lsu3shell_basis.nhw_n(in);

          unsigned int alpha_p_max = lsu3shell_basis.getMult_p(ip);
          unsigned int alpha_n_max = lsu3shell_basis.getMult_n(in);

          HalfInt Sp(lsu3shell_basis.getProtonSU3xSU2(ip).S2,2);
          HalfInt Sn(lsu3shell_basis.getNeutronSU3xSU2(in).S2,2);

          for (int iwpn = lsu3shell_basis.blockBegin(ipin_block); iwpn < lsu3shell_basis.blockEnd(ipin_block); ++iwpn) 
            {
              const auto& omega_pn = lsu3shell_basis.getOmega_pn(ip, in, iwpn);
              HalfInt S(omega_pn.S2,2);
              u3::U3 omega(Nex+Nsigma0, {omega_pn.lm,omega_pn.mu});
              unsigned int dim=alpha_n_max*alpha_p_max*omega_pn.rho;
              u3shell::U3SPN omegaSPN({omega,S},Sp,Sn);
              basis_dimensions[omegaSPN]+=dim;
            }
        }
    }
    return basis_dimensions;
  }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Removes irreps with center of mass contamination based on counting arguments
  std::map<u3shell::U3SPN, unsigned int>
  EliminateCMContributions(
    const int Nmax,
    const HalfInt& Nsigma0,
    const std::map<u3shell::U3SPN, unsigned int>& u3spn_dimensions
    )
  {
    std::map<u3shell::U3SPN, unsigned int> u3spn_cmf_dimensions = u3spn_dimensions;
  
    // Iterate through the lsu3shell basis and remove CM contaminated states
    for(unsigned int Nex=0; Nex<=Nmax; ++Nex)
    {
      //Makes use of the fact that the irreps within the map are first ordered by Nex
      for(const auto& [irrep_prime,dimension] : u3spn_cmf_dimensions)
        {
          auto Nex_prime = static_cast<unsigned int>(irrep_prime.N()-Nsigma0);
          if(Nex_prime >= Nex)
            break;

          unsigned int Ncm = Nex-Nex_prime;

          const auto& [omega_prime,Sp,Sn,S] = irrep_prime.FlatKey();
          auto cm_contaminated_irreps =  u3::KroneckerProduct(omega_prime,u3::U3(Ncm,{Ncm,0u}));
          for(const auto& [irrep,rho_max] : cm_contaminated_irreps)
            {
              assert(rho_max==1);
              u3shell::U3SPN u3spn({{irrep,S},Sp,Sn});
              // std::cout<<u3spn.Str()<<"  "<<u3spn_cmf_dimensions[u3spn]<<" "<<dimension<<std::endl;
              assert(u3spn_cmf_dimensions[u3spn]>=dimension);
              u3spn_cmf_dimensions[u3spn]-=dimension;
            }
        }
    }
    return u3spn_cmf_dimensions;
  }


  std::map<u3shell::U3SPN, unsigned int>
  LSU3ShellCMFBasisDimensions(
    const HalfInt& Nsigma0,
    const int& Nmax, 
    const std::map<u3shell::U3SPN, unsigned int>& u3spn_dimensions
  )
  {
    return EliminateCMContributions(Nmax,Nsigma0,u3spn_dimensions);
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::map<u3shell::U3SPN, unsigned int>
  LSU3ShellCMFBasisDimensions(
    const nuclide::NuclideType& nuclide, 
    const HalfInt& Nsigma0,
    const int& Nmax
  )
  {
    // initial with all lsu3shell dimensions
    std::map<u3shell::U3SPN,unsigned int> u3spn_dimensions
      = LSU3ShellBasisDimensions(nuclide,Nsigma0,Nmax);

    return EliminateCMContributions(Nmax,Nsigma0,u3spn_dimensions);
  }


  void PrintLSU3ShellBasisInfo(
    const nuclide::NuclideType& nuclide,
    const int Nmax,
    const int Nstep
  )
  // For a given Nucleus, print basis labels including
  // Nex, x_proton S_proton x_neutron S_neutron rho x S
  {
    const auto&[Z,N]=nuclide;
    std::map<u3shell::U3SPN,unsigned int> basis_dimensions;
    for(int Nex=0; Nex<=Nmax; Nex+=Nstep)
    {
      //SU3ME/proton_neutron_ncsmSU3Basis.h
      proton_neutron::ModelSpace lsu3shell_model_space(Z,N,Nex);

      //LSU3/ncsmSU3xSU2Basis.cpp
      int idiag=0; int ndiag=1;
      lsu3::CncsmSU3xSU2Basis lsu3shell_basis(lsu3shell_model_space, idiag, ndiag);

      //  loop over (ip, in) pairs
      for (int ipin_block = 0; ipin_block < lsu3shell_basis.NumberOfBlocks(); ipin_block++)
        {
          // If block is empty, continue
          if (!lsu3shell_basis.NumberOfStatesInBlock(ipin_block)) {continue;}

          unsigned int ip = lsu3shell_basis.getProtonIrrepId(ipin_block);
          unsigned int in = lsu3shell_basis.getNeutronIrrepId(ipin_block);
          unsigned int Nex_proton = lsu3shell_basis.nhw_p(ip);
          unsigned int Nex_neutron = lsu3shell_basis.nhw_n(in);
          unsigned int Nex = Nex_proton+Nex_neutron;

          unsigned int alpha_p_max = lsu3shell_basis.getMult_p(ip);
          unsigned int alpha_n_max = lsu3shell_basis.getMult_n(in);

          auto irrep_proton = lsu3shell_basis.getProtonSU3xSU2(ip);
          auto irrep_neutron = lsu3shell_basis.getNeutronSU3xSU2(in);
          HalfInt S_proton(irrep_proton.S2,2);
          HalfInt S_neutron(irrep_neutron.S2,2);
          u3::SU3 x_proton(irrep_proton.lm,irrep_proton.mu);
          u3::SU3 x_neutron(irrep_neutron.lm,irrep_neutron.mu);

          for (int iwpn = lsu3shell_basis.blockBegin(ipin_block); iwpn < lsu3shell_basis.blockEnd(ipin_block); ++iwpn)
            {
              const auto& omega_pn = lsu3shell_basis.getOmega_pn(ip, in, iwpn);
              HalfInt S(omega_pn.S2,2);
              u3::SU3 x(omega_pn.lm,omega_pn.mu);
              unsigned int dim=alpha_n_max*alpha_p_max*omega_pn.rho;
              std::cout<<fmt::format("{}  {} {} {}  {} {} {}  {} {}  {}",
                Nex,
                Nex_proton,x_proton.Str(),S_proton,
                Nex_neutron,x_neutron.Str(),S_neutron,
                x.Str(),S,
                dim
              )<<std::endl;

            }
        }
    }

  }


}//lsu3shell
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace lgi
{

MultiplicityTaggedLGIVector GetLGIVector(
    const std::map<u3shell::U3SPN, unsigned int>& lsu3shell_cmf_dimensions,
    const HalfInt& Nsigma0,
    const unsigned int Nmax
  )
{
  // Initialize lgi_dimension with cmf U3SPN dimensions from lsu3shell basis
  std::map<u3shell::U3SPN, unsigned int> lgi_dimensions = lsu3shell_cmf_dimensions;

  // Iterate through basis and identify LGI dimension by substracting
  // U(3) irreps obtained by laddering from lower grade LGI.
  for(const auto& [lgi,dimension] : lgi_dimensions)
    {
      HalfInt Sp(lgi.Sp()),Sn(lgi.Sn()), S(lgi.S());
      int Nn_max = Nmax - int(lgi.N() - Nsigma0);

      // Eventually be replaced by single call.
      const auto& sigma=lgi.U3();
      bool subspace_labels_only = true;
      sp3r::Sp3RSpace sp3r_space(sigma,Nn_max,subspace_labels_only);

      for(const auto& subspace : sp3r_space)
        {
          if(subspace.U3()==sigma)
            continue;

          u3shell::U3SPN omegaSpSnS(subspace.U3(),Sp,Sn,S);
          lgi_dimensions[{subspace.U3(),Sp,Sn,S}] -= subspace.upsilon_max()*dimension;
        }
    }

  //Create LGI vector used in SpNCCI basis construction
  MultiplicityTaggedLGIVector lgi_vector;
  for(const auto& [lgi_u3spn,dim] : lgi_dimensions)
    {
      int Nsex=int(lgi_u3spn.N()-Nsigma0);
      if(Nsex%2==Nmax%2)
      {
        lgi::LGI lgi(lgi_u3spn,Nsex);
        lgi_vector.emplace_back(lgi,dim);

      }
    }

  return lgi_vector;
}

}// end namespace
