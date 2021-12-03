/****************************************************************
  lgi_gen.h

  Generate lgi expansion in lsu3shell basis 

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  11/1/21 (aem): Created.
****************************************************************/

#ifndef LGI_GEN_H_
#define LGI_GEN_H_

#include <vector>
#include <Eigen/Eigen>
#include "mpi.h"
#include "am/am.h"
#include "lgi/lgi.h"


namespace lgi
{

basis::OperatorBlocks<double> generate_lgi_expansion(
  const nuclide::NuclideType& nuclide,
  const int Nsigma_max,
  const lgi::MultiplicityTaggedLGIVector& lgi_vector,
  const std::string& operator_dir,
  const MPI_Comm world_comm=MPI_COMM_WORLD
  );

basis::OperatorBlocks<double> generate_lgi_expansion(
  const nuclide::NuclideType& nuclide,
  const int Nsigma_max,
  const std::string& operator_dir,
  const MPI_Comm world_comm=MPI_COMM_WORLD
  );
  // If lgi_vector not given, generate lgi_vector for full
  //  Nmax model sapce with Nmax=Nsigma_max.


  void ComputeSeeds(
    const nuclide::NuclideType& nuclide,
    const int Nsigma_max, const int N1v,
    const std::pair<MultiplicityTagged<lgi::LGI>,MultiplicityTagged<lgi::LGI>>& lgi_pair,
    const basis::OperatorBlock<double>& lgi_expansion_bra,
    const basis::OperatorBlock<double>& lgi_expansion_ket,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& unit_tensor_labels,
    const std::string& operator_dir,
    const MPI_Comm world_comm = MPI_COMM_WORLD,
    const bool& restrict_op_J0=false
    );


}
#endif
