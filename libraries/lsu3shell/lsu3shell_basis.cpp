/****************************************************************
  lsu3shell_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lsu3shell/lsu3shell_basis.h"


#include <fstream>
#include <iostream>
#include <algorithm>

#include "cppformat/format.h"

namespace lsu3shell
{

  void 
  ReadLSU3Basis(
      HalfInt Nsigma_0, 
      const std::string& filename, 
      LSU3BasisTable& lsu3_basis_table, 
      std::map<u3shell::U3SPN,int>& subspace_dimensions
    )
  {
    std::ifstream basis_stream(filename.c_str());
    if(not basis_stream)
      std::cout<<fmt::format("File {} did not open",filename)<<std::endl;
    std::string line;
    while(std::getline(basis_stream,line))
      {
        if(not std::isdigit(line[0]))
          continue;
        std::istringstream line_stream(line);
        // alpha_n_max and alpha_p_max correspond multiplicities arising in the coupling of protons among shells and 
        // neutrons among shells. rho0_max is outer multplicity of coupling proton to neutron. 
        int Np, Nn, Nex,twice_Sp, twice_Sn, twice_S, lambda, mu, ip,in,alpha_p_max,alpha_n_max,rho0_max,lambda_p,mu_p,lambda_n,mu_n;
        line_stream
          >>ip>>alpha_p_max>>Np>>lambda_p>>mu_p>>twice_Sp 
          >>in>>alpha_n_max>>Nn>>lambda_n>>mu_n>>twice_Sn 
          >>rho0_max>>lambda >> mu>>twice_S;

        //conversions
        Nex=Nn+Np;
        HalfInt Sp=HalfInt(twice_Sp,2);
        HalfInt Sn=HalfInt(twice_Sn,2);
        HalfInt S=HalfInt(twice_S,2);
        u3::SU3 x(lambda,mu);
        // std::cout<<fmt::format("Nex {}  Nsigma_0 {}  x {}",Nex, Nsigma_0,x.Str())<<std::endl;
        u3::U3 omega(Nsigma_0+Nex,x);
        u3::U3S omegaS(omega,S);
        u3shell::U3SPN omegaSPN(omegaS,Sp,Sn);

        int start_index=subspace_dimensions[omegaSPN];
        int dim=alpha_n_max*alpha_p_max*rho0_max;
        subspace_dimensions[omegaSPN]+=dim;
        
        LSU3BasisGroup mult_group(omegaSPN,dim,start_index);
        lsu3_basis_table.push_back(mult_group);
      } 
    basis_stream.close();
  }

}// end namespace
