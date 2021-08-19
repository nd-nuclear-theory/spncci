/****************************************************************
  parameters.h

  Code run parameters.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  6/17/17 (mac): Created, extracted from spncci.cpp.
  6/26/17 (mac): Remove dependence of parameters on command line arguments.
  10/24/17 (mac): Add eigensolver parameters and pass-through interaction
    parameters to input file.
****************************************************************/

#ifndef SPNCCI_SPNCCI_PARAMETERS_H_
#define SPNCCI_SPNCCI_PARAMETERS_H_

#include <array>
#include <string>
#include <vector>

#include "am/halfint.h"

namespace spncci
{
  //TODO extract out lsu3shell file information.  Only used in explicit construction.

  struct RunParameters
  // Structure to store input parameters for run.
  {

    // constructor
    RunParameters();

    // basis parameters
    int A;
    HalfInt Nsigma0;
    int Nsigmamax;
    int N1v;
    int Nmax;
    int gex;  // assumes single parity runs for now (to revisit later)

    // run mode
    bool count_only;  // basis counting run
    bool transform_lgi; //Apply basis transformation/truncation to lgi and seed rmes.
    // upstream information
    std::array<int,2> nuclide;  // (Z,N): proton and neutron numbers
    std::string interaction_name;
    bool use_coulomb;

    // filenames
    std::string lsu3shell_rme_directory;
    std::string lsu3shell_basis_filename;
    std::string Brel_filename;
    std::string Arel_filename;
    std::string Nrel_filename;
    std::string relative_unit_tensor_filename_template;

    // many-body problem
    std::string observable_directory;
    std::vector<std::string> observable_filenames;  // first observable is used as Hamiltonian
    std::vector<int> observable_J0_values;
    int num_observables;
    std::vector<HalfInt> J_values;
    std::vector<double> hw_values;

    // eigensolver
    int num_eigenvalues;
    int eigensolver_num_convergence;  // whatever exactly this is...
    int eigensolver_max_iterations;
    double eigensolver_tolerance;

  };

}  // namespace

#endif
