/****************************************************************
  null_solver.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "lgi/null_solver.h"

#include <iostream>

#include "fmt/format.h"


namespace lgi
{

  void FindNullSpaceSVD(Eigen::MatrixXd& A, Eigen::MatrixXd& null_vectors, double threshold, bool verbose)
  {
    assert(A.rows()>=A.cols());

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,Eigen::ComputeFullV);  // for tall A, full and thin V should be identical

    Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType w = svd.singularValues();
    Eigen::MatrixXd V = svd.matrixV();

    if (verbose)
      {
        std::cout << std::endl
                  << "A: " << std::endl
                  << mcutils::FormatMatrix(A,"e") << std::endl;
        std::cout << std::endl
                  << "w values: " << std::endl
                  << mcutils::FormatMatrix(w,"e") << std::endl;
        std::cout << std::endl
                  << "V: " << std::endl
                  << mcutils::FormatMatrix(null_vectors,"e") << std::endl;
      }

    // extract nullity
    svd.setThreshold(threshold);
    int rank = svd.rank();
    int dimension = A.cols();
    int nullity = dimension - rank;

    // extract null vectors
    assert(V.rows()==V.cols());
    null_vectors = V.block(0,rank,dimension,dimension-rank);

    if (verbose)
      {
        std::cout << std::endl
                  << "null_vectors: " << std::endl
                  << mcutils::FormatMatrix(null_vectors,"e") << std::endl;
      }

  }

}
