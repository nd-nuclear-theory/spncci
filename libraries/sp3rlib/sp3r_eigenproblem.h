/****************************************************************
  sp3r_eigenproblem.h

  Control code for eigensolver

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/
#ifndef SP3R_EIGENPROBLEM_H_
#define SP3R_EIGENPROBLEM_H_

#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include "Spectra/SymEigsSolver.h"  // from spectra

#include "fmt/format.h"
#include "mcutils/eigen.h"

#include "basis/basis.h"
#include "basis/degenerate.h"
#include "basis/operator.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/config_parameters.h"


namespace sp3r
{

// Solve Eigenproblem
void SolveEigenproblem(
  basis::OperatorBlock<MatrixFloatType>& operator_matrix,
  int num_eigenvalues,
  int eigensolver_num_convergence,
  int eigensolver_max_iterations,
  double eigensolver_tolerance,
  sp3r::Vector& eigenvalues,
  sp3r::Matrix& eigenvectors,
  bool verbose=true
);
// Solve the Hamiltonian matrix (in a single J-space) for energy eigenvalues and vectors 
//
// Arguments:
//   hamiltonian_matrix (input): Hamiltonian matrix block for given J-subspace
//   num_eigenvalues (input): number of eigenvalues to solve for
//   eigensolver_num_convergence (input): Arnoldi eigensolver
//        parameter (whatever exactly this is...)
//   eigensolver_max_iterations (input): Arnoldi eigensolver maximum iterations
//   eigensolver_tolerance (input): Arnoldi tolerance parameters
//   eigenvalues (output): eigenvalues for J-subspace
//   eigenvectors (output): eigenvectors for J-subspace (as column vectors)

void WriteMatrixToFile(basis::OperatorBlock<MatrixFloatType>& operator_matrix,  double identifier);
 // Write Hamiltonian matrix to file 

}

#endif