/****************************************************************
  get_u3s_subspaces_su3pn.cpp

  DEPRECATED from the beginning (!): This is ad hoc adaptation of
  get_u3s_subspaces to include SUp(3) and SUn(3) in the binning key.

  Syntax:

    get_u3s_subspaces basis_filename dimensions_filename

      basis_filename: input LSU3shell basis filename, typically
      Z{}-N{}-Nmax{:02d}_basis.dat

      dimensions_filename: output U3S dimensions filename, typically
      u3s_subspace_labels_Z{}_N{}_Nmax{:02d}.dat

  Example:
 
    get_u3s_subspaces ~/results/mcaprio/lsu3shell/runmac0549/results/Z3-N3-Nmax06_basis.dat u3s_subspace_labels_Z3_N3_Nmax06.dat

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  03/21/21 (mac): Created, by patching get_u3s_subspaces key definition and
  output statement.

****************************************************************/
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>

// #include "spncci/computation_control.h"
#include "fmt/format.h"
#include "am/halfint.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
// #include "u3shell/u3spn_scheme.h"
// #include "spncci/results_output.h"

namespace lsu3shell
{
  void ReadLSU3ShellBasisAndAccumulate(
      const std::string& filename, 
      std::map<std::tuple<int,int,int,int,int,int,int,int,int,int,int,int>,int>& u3s_subspace_labels_set  // patched for pn (mac)
    )
  // Read basis states labels from LSU3Shell basis file and accumulate by N(lambda,mu)S labels
  // Based on ReadLSU3ShellBasis in lsu3shell/lsu3shell_basis.cpp
  {
    // open basis file
    std::ifstream basis_stream(filename.c_str());
    if(not basis_stream)
      std::cout<<fmt::format("File {} did not open",filename)<<std::endl;

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

        int dim=alpha_n_max*alpha_p_max*rho0_max;
        std::tuple<int,int,int,int,int,int,int,int,int,int,int,int>
          labels(Nex,lambda,mu,twice_Sp,twice_Sn,twice_S,Np,lambda_p,mu_p,Nn,lambda_n,mu_n);  // patched for pn (mac); put wp and wn last to preserve existing block structure in sort
        u3s_subspace_labels_set[labels]+=dim;
      }

    // done with basis file
    basis_stream.close();
  }

void WriteU3SLabels(std::string filename, const std::map<std::tuple<int,int,int,int,int,int,int,int,int,int,int,int>,int>& u3s_subspace_labels_set)  // patched for pn (mac)
// Write U3S labels and dimensions.
  {
    std::ofstream outfile;
    outfile.open(filename);
    for(auto it=u3s_subspace_labels_set.begin(); it!=u3s_subspace_labels_set.end(); ++it)
      {
        int Np,lambda_p,mu_p,Nn,lambda_n,mu_n;  // patched for pn (mac)
        int Nex,lambda,mu,twice_S,twice_Sp,twice_Sn;
        std::tie(
            Nex,lambda,mu,twice_Sp,twice_Sn,twice_S,
            Np,lambda_p,mu_p,Nn,lambda_n,mu_n  // patched for pn (mac)
          )=it->first;
        int dim=it->second;
        outfile<<fmt::format(
            "{:2d} {:3d} {:3d} {:2d} {:3d} {:3d} {:2d} {:3d} {:3d} {:3d} {:3d} {:3d} {:4d}",  // patched for pn (mac)
            Np,lambda_p,mu_p,Nn,lambda_n,mu_n,  // patched for pn (mac)
            Nex,lambda,mu,twice_Sp,twice_Sn,twice_S,dim
          )<<std::endl;
      }
    outfile.close();
  }

}//namespace

int main(int argc, char **argv)
{
  
  if(argc<2+1)
  {
    // std::cout<<"Syntax: nucleus_label filename"<<std::endl;
    exit(0);
  }

  std::string lsu3shell_filename=argv[1];
  std::string u3s_filename=argv[2];

    // SU(3) caching
  u3::U3CoefInit();

  
  std::map<std::tuple<int,int,int,int,int,int,int,int,int,int,int,int>,int> u3s_subspace_labels_set;  // patched for pn (mac)
  lsu3shell::ReadLSU3ShellBasisAndAccumulate(lsu3shell_filename, u3s_subspace_labels_set);
  
  lsu3shell::WriteU3SLabels(u3s_filename,u3s_subspace_labels_set);
  
  // std::string u3sl_filename=fmt::format("u3sl_subspace_labels_Z{}_N{}_Nmax{:02d}.dat",Z,N,Nmax);
  // lsu3shell::WriteU3SBranchedLabels(u3sl_filename,u3s_subspace_labels_set);

} 
