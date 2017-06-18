/****************************************************************
  parameters.h

  Code run parameters.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/17/17 (mac): Created, extracted from spncci.cpp.
****************************************************************/

#ifndef SPNCCI_SPNCCI_PARAMETERS_H_
#define SPNCCI_SPNCCI_PARAMETERS_H_

#include <string>

#include "am/halfint.h"

namespace spncci
{

  struct RunParameters
  // Structure to store input parameters for run.
  //
  // Data members:
  //   A (int): Atomic mass.
  //   ...
  {

    // constructor
    RunParameters(int argc, char **argv); 

    // basis parameters
    int A;
    HalfInt Nsigma0;
    int Nsigmamax;
    int N1v;
    int Nmax;

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
    std::vector<int> observable_Jvalues;
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
