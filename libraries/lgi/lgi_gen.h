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
#include "LSU3/ncsmSU3xSU2Basis.h"
#include "fmt/format.h"

namespace lgi
{
  inline std::string lgi_expansion_filename(const int Z, const int N, const lgi::LGI& lgi)
    {
      const auto&[Nex,sigma,Sp,Sn,S]=lgi.Key();
      std::string filename
        = fmt::format("lgi_expansion_Z{:02d}_N{:02d}_Nex{:02d}_lm{:02d}_mu{:02d}_2Sp{:02}_2Sn{:02}_2S{:02}.dat",
              Z,N,Nex,sigma.SU3().lambda(),sigma.SU3().mu(),TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S));
      return filename;
    }


basis::OperatorBlocks<double> generate_lgi_expansion(
  const nuclide::NuclideType& nuclide,
  const int Nsigma_max,
  const lgi::MultiplicityTaggedLGIVector& lgi_vector,
  const std::string& operator_dir,
  const MPI_Comm world_comm=MPI_COMM_WORLD
  );


 basis::OperatorBlock<double> generate_lgi_expansion(
    const nuclide::NuclideType& nuclide,
    const int Nsigma_max,
    const unsigned int N0,
    const CBaseSU3Irreps& baseSU3Irreps,
    const lgi::MultiplicityTaggedLGIVector& lgi_vector,
    const int lgi_index,
    const std::string& operator_dir,
    const MPI_Comm& individual_comm,
    double zero_threshold=1e-6
    );

basis::OperatorBlocks<double> generate_lgi_expansion(
  const nuclide::NuclideType& nuclide,
  const int Nsigma_max,
  const std::string& operator_dir,
  const MPI_Comm world_comm=MPI_COMM_WORLD
  );
  // If lgi_vector not given, generate lgi_vector for full
  //  Nmax model sapce with Nmax=Nsigma_max.


  std::vector<basis::OperatorBlocks<double>>
   ComputeSeeds(
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
