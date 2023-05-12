/****************************************************************
  get_u3s_subspaces.cpp

  Syntax:

    get_u3s_subspaces Z N Nmax dimensions_filename

      dimensions_filename: output U3S dimensions filename, typically
      u3s_subspace_labels_Z{}_N{}_Nmax{:02d}.dat

  Example:
 
    get_u3s_subspaces 3 3 6 u3s_subspace_labels_Z3_N3_Nmax06.dat

  Anna E. McCoy[2,3] and Mark A. Caprio[1]
  [1] University of Notre Dame
  [2] TRIUMF
  [3] Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  11/12/19 (aem): Created.
  08/29/20 (mac): Simplify command line arguments.
  03/09/22 (aem): Generate lsu3shell basis rather than read from file.

****************************************************************/
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>

#include "fmt/format.h"
#include "am/halfint.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "utilities/nuclide.h"
#include "u3ncsm/dimensions.h"

namespace lsu3shell
{

void WriteU3SLabels(
  const HalfInt& Nsigma0,
  std::string filename,
  const std::map<u3shell::U3SPN,unsigned int>& u3s_subspace_labels_set
  )
// Write U3SpSnS labels and dimensions.
  {
    std::ofstream outfile;
    outfile.open(filename);
    for(const auto& [subspace_labels,dim] : u3s_subspace_labels_set)
      {
        const auto& [omega,Sp,Sn,S] = subspace_labels.FlatKey();
        const auto& x = omega.SU3();
        outfile<<fmt::format("{:2d} {:3d} {:3d} {:3d} {:3d} {:3d} {:4d}",
          int(omega.N()-Nsigma0),x.lambda(),x.mu(),TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S),dim
        )<<std::endl;
      }
    outfile.close();
  }

void WriteU3SBranchedLabels(
  const HalfInt& Nsigma0,
  std::string filename,
  const std::map<u3shell::U3SPN,unsigned int>& u3s_subspace_labels_set
  )
  // Write U3LSpSnS labels and dimensions
  {

    std::ofstream outfile;
    outfile.open(filename);
    for(const auto& [subspace_labels,dim] : u3s_subspace_labels_set)
      {
        const auto& [omega,Sp,Sn,S] = subspace_labels.FlatKey();
        const auto& x = omega.SU3();
        int Nex = int(omega.N()-Nsigma0);

        MultiplicityTagged<unsigned int>::vector branched_states=u3::BranchingSO3(x);
        for(const auto& [L,kappa_max] : branched_states)
          {
            for(int kappa=1; kappa<=kappa_max; ++kappa)
              outfile<<fmt::format("{:2d} {:3d} {:3d} {:3d} {:3d} {:4d}",
                Nex,x.lambda(),x.mu(),kappa,L,TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S),dim)<<std::endl;
          }
      }
    outfile.close();
  }

}//namespace

int main(int argc, char **argv)
{
  
  if(argc<4+1)
  {
    std::cout<<"Syntax: Z N Nmax outputfile_name"<<std::endl;
    exit(0);
  }

  int Z = std::stoi(argv[1]);
  int N = std::stoi(argv[2]);
  int Nmax = std::stoi(argv[3]);
  nuclide::NuclideType nuclide({Z,N});
  std::string u3s_filename=argv[4];

  // // SU(3) caching
  // int max_lambda_plus_mu=39;
  // u3::U3CoefInit(max_lambda_plus_mu);

  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide);

  auto u3s_subspace_labels_set
    = lsu3shell::LSU3ShellBasisDimensions(nuclide,Nsigma0,Nmax);


  lsu3shell::WriteU3SLabels(Nsigma0,u3s_filename,u3s_subspace_labels_set);
  
  // std::string u3sl_filename=fmt::format("u3sl_subspace_labels_Z{}_N{}_Nmax{:02d}.dat",Z,N,Nmax);
  // lsu3shell::WriteU3SBranchedLabels(u3sl_filename,u3s_subspace_labels_set);

} 
