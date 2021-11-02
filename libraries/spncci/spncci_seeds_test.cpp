/****************************************************************
  spncci_seeds_test.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

 10/29/21 (aem): Created.
****************************************************************/
#include "spncci/spncci_seeds.h"

#include <cppitertools/itertools.hpp>
#include <string>
#include <mpi.h>

#include "fmt/format.h"
#include "lgi/lgi.h"
#include "utilities/nuclide.h"
#include "utilities/utilities.h"

int main(int argc, char** argv)
{

 if (argc<4)
  {
    std::cout<<"Usage: recurrence seeds_test <Z> <N> <Nsigma_max>"<<std::endl;
  }
  MPI_Init(&argc, &argv);
  

  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);
  int Nsigma_max=std::stoi(argv[3]);
  
  nuclide::NuclideType nuclide({3,3});
  bool intrinsic=true;
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
  std::string spncci_root_dir=get_spncci_project_root_dir();
  std::string operator_dir = fmt::format("{}/spncci/data/relative_operators/",spncci_root_dir);
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //Generate LGI vector by finding possible cmf LGI by counting arguments 
  lgi::MultiplicityTaggedLGIVector lgi_vector = lgi::get_lgi_vector(nuclide,Nsigma0,Nsigma_max);

  int eigensolver_max_iterations=1000;
  double eigensolver_tolerance=1e-8;
  MPI_Comm world_comm = MPI_COMM_WORLD;
  std::cout<<"about to generate lgi expansion"<<std::endl;
  spncci::seeds::generate_lgi_expansion(
    nuclide,Nsigma_max,operator_dir,
    eigensolver_max_iterations,eigensolver_tolerance,world_comm
  );
  

  MPI_Finalize();
}
