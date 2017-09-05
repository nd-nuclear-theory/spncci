/****************************************************************
  eigenproblem.h

  Control code for eigensolver.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/6/17 (mac): Extracted from computation_control.

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
      const HalfInt& J,
      int num_eigenvalues,
      int eigensolver_num_convergence,  // whatever exactly this is...
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      spncci::Vector& eigenvalues,  // eigenvalues for J-subspace
      spncci::Matrix& eigenvectors  // eigenvectors for J-subspace
    );
  // Solve the Hamiltonian matrix (in a single J-space) for energy eigenvalues and vectors 


}  // namespace

#endif
