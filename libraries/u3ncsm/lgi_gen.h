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
#include "u3ncsm/u3ncsm_interface.h"

namespace lgi{

 // basis::OperatorBlock<double> generate_lgi_expansion(
 //    const nuclide::NuclideType& nuclide,
 //    const int Nsigma_max,
 //    const unsigned int N0,
 //    const CBaseSU3Irreps& baseSU3Irreps,
 //    const lgi::MultiplicityTaggedLGIVector& lgi_vector,
 //    const int lgi_index,
 //    const std::string& operator_dir,
 //    const MPI_Comm& individual_comm,
 //    double zero_threshold=1e-6
 //    );

basis::OperatorBlock<double> GenerateLGIExpansion(
    const nuclide::NuclideType& nuclide,
    const HalfInt& Nsigma0,
    const unsigned int N0,
    const CBaseSU3Irreps& baseSU3Irreps,
    const std::map<u3shell::U3SPN, unsigned int>& u3ncsm_basis_map,
    const MultiplicityTagged<lgi::LGI>& tagged_lgi,
    const std::string& operator_dir,
    const MPI_Comm& individual_comm,
    double zero_threshold = 1e-6
  );

/// Generate expansion of the lgi_index-th LGI in lgi_vector in terms of SU(3)SpSnS irreps
/// of a no-core shell model (NCSM) basis
///
/// Args:
///   nuclide : Array containing {Z,N}
///   Nsigma0 : Minimum Nsigma value in basis.
///   N0 : Minimum Pauli allowed number of quanta of a state in the NCSM basis.
///   baseSU3Irreps: Class storing allowed SU(3)xSU(2) irreps in each harmonic oscillator
///     shell for numbers of nucleons in each shell ranging from 1,...,A.
///     Defined in lsu3shell.
///   u3ncsm_basis_map : map keyed by U3SPN irreps in u3ncsm basis.
///   tagged_lgi : pair containing LGI and multiplicity  gamma_max
///   operator_dir: directory containing operator files Nrel.{PN,PPNN} and Brel.{PN,PPNN}
}
#endif
