/****************************************************************
  spectra_dense_timing.cpp

  Carry out timing trials for dense matrix diagonalization with
  Spectra solver.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/14/17 (mac): Create, based on example from spectra/README.md.
    + Update coding conventions (for includes and namespaces).
****************************************************************/

#include <iostream>

#include "eigen3/Eigen/Dense"
#include "SymEigsSolver.h"  // from spectra

#include "mcutils/profiling.h"


int main()
{

  std::vector<int> dimension_list({10,100,200,500,1000,2000,5000,10000});
  int num_eigenvalues = 3;
  int num_convergence = 2*num_eigenvalues;  // docs for SymEigsSolver say to take "ncv>=2*nev"
  int max_iterations = 100*num_eigenvalues;
  double tolerance = 1e-6;  // if too tight (default 1e-10) causes slow evaluation and convergence failure status

  typedef Eigen::MatrixXf MatrixType;
  typedef float FloatType;

  for (int dimension : dimension_list)
    {

      std::cout << "Dimension: " << dimension << std::endl;
      std::cout << std::endl;

      // define test matrix -- all ones
      MatrixType M = MatrixType::Constant(dimension,dimension,-1.);

      // define test matrix -- random symmetric
      //MatrixType A = MatrixType::Random(dimension,dimension);
      //MatrixType M = A + A.transpose();

      // Spectra test
      if (false)
        {
          std::cout << "Spectra::SymEigsSolver" << std::endl;

          // start timing
          Timer test_time;
          test_time.Start();

          // define eigensolver and compute
          Spectra::DenseSymMatProd<FloatType> matvec(M);
          Spectra::SymEigsSolver<FloatType,Spectra::SMALLEST_ALGE,Spectra::DenseSymMatProd<FloatType>> eigensolver(&matvec,num_eigenvalues,num_convergence);
          eigensolver.init();
          int nconv = eigensolver.compute(max_iterations,tolerance,Spectra::SMALLEST_ALGE);  // int maxit=1000, Scalar tol=1e-10, int sort_rule=LARGEST_ALGE

          // end timing
          test_time.Stop();
          std::cout << "(Elapsed time: " << test_time.ElapsedTime() << ")" << std::endl;

          // retrieve results
          std::cout << "Status: " << eigensolver.info() << std::endl;
          Eigen::VectorXf evalues;
          if(eigensolver.info() == Spectra::SUCCESSFUL)
            evalues = eigensolver.eigenvalues();
          std::cout << "Eigenvalues:\n" << evalues << std::endl;

          std::cout << std::endl;
          std::cout << std::endl;
        }

      // Eigen test
      if (true)
        {
          std::cout << "Eigen::SelfAdjointEigenSolver" << std::endl;

          // start timing
          Timer test_time;
          test_time.Start();

          // define eigensolver and compute
          Eigen::SelfAdjointEigenSolver<MatrixType> eigensolver(M);

          // end timing
          test_time.Stop();
          std::cout << "(Elapsed time: " << test_time.ElapsedTime() << ")" << std::endl;


          // retrieve results
          std::cout << "Status: " << eigensolver.info() << std::endl;
          Eigen::VectorXf evalues;
          if(eigensolver.info() == Eigen::Success)
            evalues = eigensolver.eigenvalues();
          evalues.resize(20);  // truncate
          std::cout << "Eigenvalues:\n" << evalues << std::endl;

          std::cout << std::endl;
          std::cout << std::endl;

        }
    }

  return EXIT_SUCCESS;
}
