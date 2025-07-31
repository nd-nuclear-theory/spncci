/****************************************************************
  sp3r_eigenproblem.cpp

  Control code for eigensolver

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/

#include "sp3rlib/sp3r_eigenproblem.h"


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
)
{
  int hamiltonian_dimension = operator_matrix.rows();
  int actual_num_eigenvalues = std::min(num_eigenvalues,hamiltonian_dimension);
  if (hamiltonian_dimension==0)
    {
      std::cout << "  Skipping space of dimension zero" << std::endl;
      return;
    }
  bool revert_to_full_solver = !(
      (num_eigenvalues<hamiltonian_dimension)
      &&(eigensolver_num_convergence<=hamiltonian_dimension)
    );

  if (revert_to_full_solver)
  // use Eigen::SelfAdjointEigenSolver
  {
    if(verbose)
      std::cout<<"  Using solver: Eigen::SelfAdjointEigenSolver"<<std::endl;

    // define eigensolver and compute
    Eigen::SelfAdjointEigenSolver<basis::OperatorBlock<sp3r::MatrixFloatType>> eigensolver(operator_matrix);

    // Verify status. Returns "Success" if computation was succesful, "NoConvergence" otherwise.
    int eigensolver_status = eigensolver.info();
    if(verbose)
    {
      std::cout
      << fmt::format("  Eigensolver reports: status {}",eigensolver_status)
      << std::endl;
    }
    assert(eigensolver_status==Eigen::Success);

    // Save eigenvalues and eigenvectors
    eigenvalues = eigensolver.eigenvalues().block(0,0,actual_num_eigenvalues,1);
    eigenvectors = eigensolver.eigenvectors().block(0,0,operator_matrix.rows(),actual_num_eigenvalues);
  }
  else
  // use Spectra::SymEigsSolver
  {
    if(verbose)
      std::cout<<"  Using solver: Spectra::SymEigsSolver"<<std::endl;

    // define eigensolver and compute
    Spectra::DenseSymMatProd<sp3r::MatrixFloatType> matvec(operator_matrix);
    Spectra::SymEigsSolver<Spectra::DenseSymMatProd<sp3r::MatrixFloatType> >
      eigensolver(
          matvec,
          num_eigenvalues,
          eigensolver_num_convergence
        );
    eigensolver.init();
    int converged_eigenvectors = eigensolver.compute(
        Spectra::SortRule::SmallestAlge,
        eigensolver_max_iterations,
        eigensolver_tolerance
      );

    // verify status
    //
    // From Spectra documentation:
    //
    //   enum  	Spectra::COMPUTATION_INFO {
    //     Spectra::SUCCESSFUL = 0,
    //     Spectra::NOT_COMPUTED,
    //     Spectra::NOT_CONVERGING,
    //     Spectra::NUMERICAL_ISSUE
    //   }

    auto eigensolver_status = eigensolver.info();
    int eigensolver_num_iterations = eigensolver.num_iterations();
    if(verbose)
    {
      std::cout
        << fmt::format("  Eigensolver reports: eigenvectors {} status {} num_iterations {}",converged_eigenvectors,eigensolver_status,eigensolver_num_iterations)
        << std::endl;
    }

    assert(converged_eigenvectors=eigensolver.eigenvalues().size());  // should this always be true?
    assert(converged_eigenvectors=eigensolver.eigenvectors().cols());  // should this always be true?
    assert(converged_eigenvectors==num_eigenvalues);  // require that expected number of eigenvectors

    // save eigenresults
    eigenvalues = eigensolver.eigenvalues();
    eigenvectors = eigensolver.eigenvectors();
  }
    if(verbose)
    {
      // diagnostic output: eigenvalues
      std::cout << fmt::format("  Eigenvalues:") << std::endl
                << mcutils::FormatMatrix(eigenvalues.transpose(),"8.5f","    ")
                << std::endl;
    }
  // check eigenvector norms
  sp3r::Vector eigenvector_norms(eigenvectors.cols());
  for (int eigenvector_index=0; eigenvector_index<actual_num_eigenvalues; ++eigenvector_index)
    {
      eigenvector_norms(eigenvector_index) = eigenvectors.col(eigenvector_index).norm();
      const sp3r::MatrixFloatType norm_tolerance=1e-8;
      assert(fabs(eigenvector_norms(eigenvector_index)-1)<norm_tolerance);
    }

  // normalize eigenvectors -- redundant with Spectra eigensolver
  for (int eigenvector_index=0; eigenvector_index<actual_num_eigenvalues; ++eigenvector_index)
    eigenvectors.col(eigenvector_index).normalize();
}


void WriteMatrixToFile(basis::OperatorBlock<MatrixFloatType>& operator_matrix, double identifier)
{
  std::string filename=fmt::format("matrix{:0.1f}.out",identifier);
  std::ofstream stream(filename);
  int rows=matrix.rows();
  stream<<rows<<std::endl;
  stream<<operator_matrix<<std::endl;
}

}