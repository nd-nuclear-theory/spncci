/****************************************************************
  lsu3shell_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lsu3shell/lsu3shell_basis.h"


#include <fstream>
#include <iostream>
#include <algorithm>

#include "fmt/format.h"

namespace lsu3shell
{

  void ReadLSU3ShellBasis(
      HalfInt Nsigma_0, 
      const std::string& filename, 
      LSU3ShellBasisTable& lsu3_basis_table,
      U3SPNBasisLSU3Labels& basis_provenance,
      u3shell::SpaceU3SPN& space
    )
  {

    // open basis file
    std::ifstream basis_stream(filename.c_str());
    if(not basis_stream)
      std::cout<<fmt::format("File {} did not open",filename)<<std::endl;

    // set up temporary storage
    std::map<u3shell::U3SPN,int> subspace_dimensions;
    std::map<u3shell::U3SPN,std::vector<LSU3ShellBasisGroupLabels>> subspace_provenances;

    // read basis table entries
    std::string line;
    while(std::getline(basis_stream,line))
      {
        // skip initial header line
        if(not std::isdigit(line[0]))
          continue;

        // parse line
        std::istringstream line_stream(line);
        // alpha_n_max and alpha_p_max correspond multiplicities arising in the coupling of protons among shells and 
        // neutrons among shells. rho0_max is outer multplicity of coupling proton to neutron. 
        int Np, Nn, Nex,twice_Sp, twice_Sn, twice_S, lambda, mu, ip,in,alpha_p_max,alpha_n_max,rho0_max,lambda_p,mu_p,lambda_n,mu_n;
        line_stream
          >>ip>>alpha_p_max>>Np>>lambda_p>>mu_p>>twice_Sp 
          >>in>>alpha_n_max>>Nn>>lambda_n>>mu_n>>twice_Sn 
          >>rho0_max>>lambda >> mu>>twice_S;

        // convert parameter data types
        Nex=Nn+Np;
        HalfInt Sp=HalfInt(twice_Sp,2);
        HalfInt Sn=HalfInt(twice_Sn,2);
        HalfInt S=HalfInt(twice_S,2);
        u3::SU3 xp(lambda_p,mu_p);
        u3::SU3 xn(lambda_n,mu_n);
        u3::SU3 x(lambda,mu);

        // perform sanity checks on input basis labels
        assert(am::AllowedTriangle(Sp,Sn,S));
        assert(u3::OuterMultiplicity(xp,xn,x)==rho0_max);

        // std::cout<<fmt::format("Nex {}  Nsigma_0 {}  x {}",Nex, Nsigma_0,x.Str())<<std::endl;
        u3::U3 omega(Nsigma_0+Nex,x);
        u3::U3S omegaS(omega,S);
        u3shell::U3SPN omegaSPN(omegaS,Sp,Sn);

        // computing indexing information
        int start_index=subspace_dimensions[omegaSPN];
        int dim=alpha_n_max*alpha_p_max*rho0_max;
        subspace_dimensions[omegaSPN]+=dim;
        
        // store entry
        LSU3ShellBasisGroupLabels lsu3_basis_group_labels(omegaSPN,ip,in,Np,Nn,Nex);
        LSU3ShellBasisGroupData mult_group(lsu3_basis_group_labels,dim,start_index);
        lsu3_basis_table.push_back(mult_group);

        // store provenence records
        std::vector<LSU3ShellBasisGroupLabels>& subspace_provenance_list = subspace_provenances[omegaSPN];
        for (int multiplicity_index=0; multiplicity_index < dim; ++multiplicity_index)
          subspace_provenance_list.push_back(lsu3_basis_group_labels);
      }

    // done with basis file
    basis_stream.close();

    // transfer subspace provenances into vector in canonical subspace order
    //
    // We can only do this now that we know all U3SPN labels, so we
    // know the final subspace indexing within the space.
    for (auto& omegaSPN_provenance : subspace_provenances)
      {

        // define aliases for key and value
        const u3shell::U3SPN& omegaSPN = omegaSPN_provenance.first;
        const std::vector<LSU3ShellBasisGroupLabels>& subspace_provenance_list = omegaSPN_provenance.second;
        
        // re-save provenance data for subspace
        basis_provenance.push_back(subspace_provenance_list);
      }

    // construct space
    space = u3shell::SpaceU3SPN(subspace_dimensions);

  }

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
}// end namespace
