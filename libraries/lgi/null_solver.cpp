/****************************************************************
  null_solver.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/

#include "lgi/null_solver.h"

#include <iostream>

#include "fmt/format.h"
#include "mcutils/eigen.h"

namespace lgi
{

  void FindNullSpaceSVD(Eigen::MatrixXd& A, Eigen::MatrixXd& null_vectors, double threshold, bool verbose)
  {
    assert(A.rows()>=A.cols());

    // Since the eigen SVD decomposition has trouble when matrix is all zeros, we include 
    // a check to see if all matrix elements are below threshold.  If so, then the null vectors are the identy matrix
    bool zero_matrix=mcutils::IsZero(A,threshold);

    if(zero_matrix)
      {

        null_vectors=Eigen::MatrixXd::Identity(A.cols(),A.cols());
        if(verbose)
          std::cout<<"Matrix was zero matrix.  Null vectors given by "<<std::endl<<null_vectors<<std::endl;
        return;
      }

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
                  << mcutils::FormatMatrix(V,"e") << std::endl;
      }

    // extract nullity
    svd.setThreshold(threshold);
    int rank = svd.rank();
    int dimension = A.cols();
    int nullity = dimension - rank;
    if(verbose)
      std::cout<<std::endl<<"rank: "<<rank<<"  dimension: "<<dimension<<"  nullity: "<<nullity<<std::endl;

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
