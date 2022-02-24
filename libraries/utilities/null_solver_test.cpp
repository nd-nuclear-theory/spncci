/****************************************************************
  null_solver_test.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  1/19/17 (aem) : Created.
****************************************************************/

#include "lgi/null_solver.h"

#include "fmt/format.h"

int main(int argc, char **argv)
{

  {
    std::cout << std::endl
              << "Test 1..." << std::endl;

    // generate test matrix
    Eigen::MatrixXd UTest(4,3);
    UTest << 1,0,0, 0,1,0, 0,0,1, 0,0,0;
    Eigen::MatrixXd WTest(3,3);
    WTest << 1,0,0, 0,0,0, 0,0,0;
    Eigen::MatrixXd VTest = Eigen::MatrixXd::Identity(3,3);
    Eigen::MatrixXd A = UTest * WTest * VTest.transpose();

    // solve!
    Eigen::MatrixXd null_vectors;
    double threshold = 1e-12;
    lgi::FindNullSpaceSVD(A,null_vectors,threshold,true);
  }

  {
    std::cout << std::endl
              << "Test 2..." << std::endl;
  
    // generate test matrix
    Eigen::MatrixXd UTest = Eigen::MatrixXd::Random(4,3);  // oops, not column orthogonal
    Eigen::MatrixXd WTest(3,3);
    WTest << 2,0,0, 0,1,0, 0,0,0;
    Eigen::MatrixXd VTest(3,3);
    VTest << 3,1,1, -1,2,1, -0.5,-2,-3.5;  // oops, need to normalize
    Eigen::MatrixXd A = UTest * WTest * VTest.transpose();

    // solve!
    Eigen::MatrixXd null_vectors;
    double threshold = 1e-12;
    lgi::FindNullSpaceSVD(A,null_vectors,threshold,true);
  }

  {
    std::cout << std::endl
              << "Test 3..." << std::endl;
  
    // generate test matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(1,1);

    // solve!
    Eigen::MatrixXd null_vectors;
    double threshold = 1e-12;
    lgi::FindNullSpaceSVD(A,null_vectors,threshold,true);
  }

  {
    std::cout << std::endl
              << "Test 3.5..." << std::endl;
  
    // generate test matrix
    // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2,2);
    Eigen::MatrixXd A(2,2);
    A << 1e-13, 1e-14, 1e-15,1e-16;

    // solve!
    Eigen::MatrixXd null_vectors;
    double threshold = 1e-12;
    lgi::FindNullSpaceSVD(A,null_vectors,threshold,true);
  }


  {
    std::cout << std::endl
              << "Test 4..." << std::endl;
  
    // generate test matrix
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(2,2);

    // solve!
    Eigen::MatrixXd null_vectors;
    double threshold = 1e-12;
    lgi::FindNullSpaceSVD(A,null_vectors,threshold,true);
    std::cout << fmt::format("{}x{}",null_vectors.rows(),null_vectors.cols()) << std::endl;
  }

}
