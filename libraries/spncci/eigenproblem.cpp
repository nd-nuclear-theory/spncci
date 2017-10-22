/****************************************************************
  computation_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/eigenproblem.h"

#include <iostream>

#include "SymEigsSolver.h"  // from spectra

#include "cppformat/format.h"
#include "mcutils/eigen.h"

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
    )
  {    

    // set up aliases
    std::cout << fmt::format("  Diagonalizing: J={}",J) << std::endl;

    // handle low-dimensional exceptions
    //
    // Determine appropriate solver:
    //
    // -- iterative solver requires number of eigenvalues to be
    // strictly less than matrix dimension.
    //
    // -- iterative solver also requires the ncv parameter to be less
    // than or equal to the matrix dimension ("ncv must satisfy nev <
    // ncv <= n, n is the size of matrix")

    int hamiltonian_dimension = hamiltonian_matrix.rows();
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

        std::cout << "  Using solver: Eigen::SelfAdjointEigenSolver" << std::endl;

        // define eigensolver and compute
        Eigen::SelfAdjointEigenSolver<spncci::OperatorBlock> eigensolver(hamiltonian_matrix);

        // verify status
        //
        // From Eigen documentation:
        //
        //   ComputationInfo info 	( 		) 	const
        //   	inline
        //   
        //   Reports whether previous computation was successful.
        //   
        //   Returns
        //       Success if computation was succesful, NoConvergence otherwise. 


        int eigensolver_status = eigensolver.info();
        std::cout
          << fmt::format("  Eigensolver reports: status {}",eigensolver_status)
          << std::endl;
        assert(eigensolver_status==Eigen::Success);

        // save eigenresults
        eigenvalues = eigensolver.eigenvalues();
        eigenvectors = eigensolver.eigenvectors();
      }
    else
      // use Spectra::SymEigsSolver
      {

        std::cout << "  Using solver: Spectra::SymEigsSolver" << std::endl;

        // define eigensolver and compute
        Spectra::DenseSymMatProd<spncci::MatrixFloatType> matvec(hamiltonian_matrix);
        Spectra::SymEigsSolver<spncci::MatrixFloatType,Spectra::SMALLEST_ALGE,Spectra::DenseSymMatProd<spncci::MatrixFloatType>>
          eigensolver(
              &matvec,
              num_eigenvalues,
              eigensolver_num_convergence
            );
        eigensolver.init();
        int converged_eigenvectors = eigensolver.compute(
            eigensolver_max_iterations,
            eigensolver_tolerance,
            Spectra::SMALLEST_ALGE  // sorting rule
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

        int eigensolver_status = eigensolver.info();
        int eigensolver_num_iterations = eigensolver.num_iterations();
        std::cout
          << fmt::format("  Eigensolver reports: eigenvectors {} status {} num_iterations {}",converged_eigenvectors,eigensolver_status,eigensolver_num_iterations)
          << std::endl;
        assert(converged_eigenvectors=eigensolver.eigenvalues().size());  // should this always be true?
        assert(converged_eigenvectors=eigensolver.eigenvectors().cols());  // should this always be true?
        assert(converged_eigenvectors==num_eigenvalues);  // require that expected number of eigenvectors

        // save eigenresults
        eigenvalues = eigensolver.eigenvalues();
        eigenvectors = eigensolver.eigenvectors();
      }

    // diagnostic output: eigenvalues
    std::cout << fmt::format("  Eigenvalues (J={}):",J) << std::endl
              << mcutils::FormatMatrix(eigenvalues.transpose(),"8.5f","    ")
              << std::endl;
    
    // check eigenvector norms
    spncci::Vector eigenvector_norms(eigenvectors.cols());
    for (int eigenvector_index=0; eigenvector_index<actual_num_eigenvalues; ++eigenvector_index)
      {
        eigenvector_norms(eigenvector_index) = eigenvectors.col(eigenvector_index).norm();
        const spncci::MatrixFloatType norm_tolerance=1e-8;
        assert(fabs(eigenvector_norms(eigenvector_index)-1)<norm_tolerance);
      }

    // normalize eigenvectors -- redundant with Spectra eigensolver
    for (int eigenvector_index=0; eigenvector_index<actual_num_eigenvalues; ++eigenvector_index)
      eigenvectors.col(eigenvector_index).normalize();
  }


}  // namespace
