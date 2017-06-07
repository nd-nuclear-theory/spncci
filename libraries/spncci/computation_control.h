/****************************************************************
  computation_control.h

  High-level control code for computations in SpNCCI calculation
  programs.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).

****************************************************************/

#ifndef SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_
#define SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_

#include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "spncci/spncci_basis.h"
#include "spncci/unit_tensor.h"
#include "sp3rlib/u3.h"
#include "u3shell/relative_operator.h"
// #include "u3shell/upcoupling.h"
#include "u3shell/u3spn_scheme.h"
#include "spncci/branching_u3s.h"
#include "spncci/branching_u3lsj.h"

namespace spncci
{

  // convenience typedef for use in iteration over J sectors
  typedef std::pair<HalfInt,HalfInt> JPair;

  void 
  ConstructBranchedObservables(
    const spncci::SpaceU3S& space_u3s,
    const std::vector<std::vector<spncci::SectorLabelsU3S>>& observable_sectors_u3s,
    const std::vector<basis::MatrixVector>& observable_matrices_u3s,
    std::map<HalfInt,spncci::SpaceLS>& spaces_lsj,
    int num_observables,
    const std::vector<HalfInt>& J_values,
    const std::vector<int>& observable_Jvalues,
    std::vector<std::map<spncci::JPair,spncci::MatrixType>>& observable_matrices
    );
  // Construct J branched observable matrices 


  void 
  SolveHamiltonian(
      const spncci::MatrixType& hamiltonian_matrix,
      const HalfInt& J,
      int num_eigenvalues,
      int eigensolver_num_convergence,  // whatever exactly this is...
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      std::map<HalfInt,Eigen::VectorXd>& eigenvalues,  // map: J -> eigenvalues
      std::map<HalfInt,spncci::MatrixType>& eigenvectors  // map: J -> eigenvectors
    );
  // Solve the hamiltonian matrix for energy eigenvalues and vectors 


}  // namespace

#endif
