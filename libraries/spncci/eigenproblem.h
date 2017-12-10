/****************************************************************
  eigenproblem.h

  Control code for eigensolver.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/6/17 (mac): Extracted from computation_control.
  12/10/17 (mac): Clean comments and remove J argument.

****************************************************************/

#ifndef SPNCCI_SPNCCI_EIGENPROBLEM_H_
#define SPNCCI_SPNCCI_EIGENPROBLEM_H_

#include "am/halfint.h"
#include "spncci/spncci_common.h"

namespace spncci
{
 void 
  SolveHamiltonian(
      const spncci::OperatorBlock& hamiltonian_matrix,
      int num_eigenvalues,
      int eigensolver_num_convergence,
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      spncci::Vector& eigenvalues,
      spncci::Matrix& eigenvectors
    );
 // Solve the Hamiltonian matrix (in a single J-space) for energy eigenvalues and vectors 
 //
 // Arguments:
 //   hamiltonian_matrix (input): Hamiltonian matrix block for given J-subspace
 //   num_eigenvalues (input): number of eigenvalues to solve for
 //   eigensolver_num_convergence (input): Arnoldi eigensolver
 //     parameter (whatever exactly this is...)
 //   eigensolver_max_iterations (input): Arnoldi eigensolver maximum iterations
 //   eigensolver_tolerance (input): Arnoldi tolerance parameters
 //   eigenvalues (output): eigenvalues for J-subspace
 //   eigenvectors (output): eigenvectors for J-subspace (as column vectors)

}  // namespace

#endif
