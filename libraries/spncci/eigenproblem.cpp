/****************************************************************
  eigenproblem.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT
****************************************************************/

#include "spncci/eigenproblem.h"

#include <iostream>
#include <fstream>
#include "Spectra/SymEigsSolver.h"  // from spectra

#include "fmt/format.h"
#include "mcutils/eigen.h"

namespace spncci
{
  void
  SolveEigenproblem(
      const spncci::OperatorBlock& matrix,
      int num_eigenvalues,
      int eigensolver_num_convergence,
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      spncci::Vector& eigenvalues,
      spncci::Matrix& eigenvectors,
      bool verbose
    )
  {

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

    int hamiltonian_dimension = matrix.rows();
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
          std::cout << "  Using solver: Eigen::SelfAdjointEigenSolver" << std::endl;

        // define eigensolver and compute
        Eigen::SelfAdjointEigenSolver<spncci::OperatorBlock> eigensolver(matrix);

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
        if(verbose)
          {
            std::cout
            << fmt::format("  Eigensolver reports: status {}",eigensolver_status)
            << std::endl;
          }
        assert(eigensolver_status==Eigen::Success);

        // save eigenresults
        eigenvalues = eigensolver.eigenvalues().block(0,0,actual_num_eigenvalues,1);
        eigenvectors = eigensolver.eigenvectors().block(0,0,matrix.rows(),actual_num_eigenvalues);
      }
    else
      // use Spectra::SymEigsSolver
      {
        if(verbose)
          std::cout << "  Using solver: Spectra::SymEigsSolver" << std::endl;

        // define eigensolver and compute
        Spectra::DenseSymMatProd<spncci::MatrixFloatType> matvec(matrix);
        Spectra::SymEigsSolver<Spectra::DenseSymMatProd<spncci::MatrixFloatType> >
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
        // TODO: Mark, what is going on here.  Should converged eigenvectors equal number of eigenvalues?

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


 void WriteMatrixToFile(spncci::OperatorBlock& matrix, double hw)
    {
      std::string filename=fmt::format("matrix{:0.1f}.out",hw);
      std::ofstream stream(filename);
      int rows=matrix.rows();
      stream<<rows<<std::endl;
      stream<<matrix<<std::endl;

      // for(int row=0; row<rows; ++row)
      //   for(int col=0; col<rows; ++col)
      //     {
      //       stream<<matrix(row,col);
      //     }

      // stream.close();
    }


}  // namespace
