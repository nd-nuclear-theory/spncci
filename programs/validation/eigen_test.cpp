/****************************************************************
  eigen_test.cpp

  Simple testbed program for Eigen.

  Compilation and linkage:
    g++ -o eigen_test eigen_test.cpp -I${HOME}/local/opt/eigen-3.0.3/include

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  SPDX-License-Identifier: MIT

  3/11/16 (aem,mac): Created.
****************************************************************/



#include <Eigen/Eigen>
#include <Eigen/Eigenvalues>
#include <iostream>


void EigenSolverTest()
// http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
//
// Add namespace qualifiers.
{
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4f> es;
  Eigen::Matrix4f X = Eigen::Matrix4f::Random(4,4);
  Eigen::Matrix4f A = X + X.transpose();
  es.compute(A);
  std::cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << std::endl;
  es.compute(A + Eigen::Matrix4f::Identity(4,4)); // re-use es to compute eigenvalues of A+I
  std::cout << "The eigenvalues of A+I are: " << es.eigenvalues().transpose() << std::endl;
}

int main(int argc, char **argv)
{
  EigenSolverTest();
}

