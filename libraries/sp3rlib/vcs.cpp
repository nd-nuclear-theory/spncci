/****************************************************************
  vcs.cpp

  Anna E. McCoy
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/
#include "sp3rlib/vcs.h"

#include <numeric>
#include "fmt/format.h"
#include <Eigen/Eigenvalues>
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3boson.h"
#include "mcutils/eigen.h"
#include "cppitertools/itertools.hpp"
#include "mcutils/eigen.h"

namespace vcs{

  // Temporary wrapper to construct K matrix map until we can completley eliminate 
  // K_matrix_map from spncci code and only use K_matrices stored as a part of the basis
  // indexing (Sp3RSpace).
  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map)
    {
      for(const auto& subspace : irrep)
        {
          K_matrix_map[subspace.omega()]=subspace.K_matrix();
        }
    }


  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map, vcs::MatrixCache& Kinv_matrix_map)
    {
      for(const auto& subspace : irrep)
        {
          K_matrix_map[subspace.omega()]=subspace.K_matrix();
          Kinv_matrix_map[subspace.omega()]=subspace.Kinv_matrix();
        }

    }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Uses U3Boson basis construction.
// Note: for numerical stability, the Smatrix recurrence used here differs from the recurrence relation given in,
// e.g.,
//    D. J. Rowe, A. E. McCoy and M. A. Caprio. Phys. Script. 91 (2016) 033003.
//    D. J. Rowe, J. Math Phys. 25 (1984) 2662.
//    D. J. Rowe, B. G. Wybourne and P.H. Butler. J. Phys. A 18 (1985) 939.
//
// by a factor of 8, i.e., S_Rowe = 16S_sp3rlib.  Consequently, after taking the squareroot of S to get K,
// we multiply K by a factor of 2^Nn.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unordered_map<u3::U3,basis::OperatorBlock<double>>
GenerateSmatrices(const u3::U3& sigma, const u3boson::U3BosonSpace& space)
  {
    std::unordered_map<u3::U3,basis::OperatorBlock<double>> SMatrices;
    SMatrices[sigma] = basis::OperatorBlock<double>::Identity(1,1);

    // Skip bra_index=0 (omega=sigma) since that's already accounted for above.
    for(int bra_index=1; bra_index<space.size(); ++bra_index)
      {
        const auto& bra_subspace = space.GetSubspace(bra_index);
        const auto& omega_bra=bra_subspace.omega();
        assert(omega_bra!=sigma);
        int Nn = int(omega_bra.N()-sigma.N());
        int bra_dimension = bra_subspace.dimension();
        auto omega_ket_vector = u3::KroneckerProduct(omega_bra,{-2,{0,2}});
        SMatrices[omega_bra]=basis::OperatorBlock<double>::Zero(bra_dimension,bra_dimension);
        for(const auto& [omega_ket,dummy] : omega_ket_vector)
          {
            int ket_index = space.LookUpSubspaceIndex(omega_ket);
            if(ket_index==-1)
              continue;

            const auto& ket_subspace = space.GetSubspace(ket_index);
            int ket_dimension = ket_subspace.dimension();

            basis::OperatorBlock<double> Coef1(bra_dimension,ket_dimension);
            basis::OperatorBlock<double> Coef2(ket_dimension,bra_dimension);
            for(int bra_state_index=0; bra_state_index<bra_subspace.size(); ++bra_state_index)
              for(int ket_state_index=0; ket_state_index<ket_subspace.size(); ++ket_state_index)
                {
                  const auto& bra_state = bra_subspace.GetState(bra_state_index);
                  const auto& ket_state = ket_subspace.GetState(ket_state_index);
                  const u3::U3& n_bra = bra_state.n();
                  const u3::U3& n_ket = ket_state.n();
                  const int rho_max_bra = bra_state.rho_max();
                  const int rho_max_ket = ket_state.rho_max();
                  double DeltaOmega = vcs::Omega(n_bra,omega_bra)-vcs::Omega(n_ket,omega_ket);

                  for(int rho_bra=1; rho_bra<=rho_max_bra; ++rho_bra)
                    for(int rho_ket=1; rho_ket<=rho_max_ket; ++rho_ket)
                    {
                      int row = bra_subspace.GetStateOffset(bra_state_index,rho_bra);
                      int col = ket_subspace.GetStateOffset(ket_state_index,rho_ket);
                      double boson_rme = u3boson::U3BosonCreationRME(sigma,n_bra,rho_bra,omega_bra,sigma,n_ket,rho_ket,omega_ket);
                      Coef1(row,col) = DeltaOmega*boson_rme;
                      Coef2(col,row) = boson_rme;
                    }
                }

            SMatrices[omega_bra]+=Coef1*SMatrices[omega_ket]*Coef2;
          }

          SMatrices[omega_bra]=SMatrices[omega_bra]/(8*Nn);
      }
    return SMatrices;
  }

vcs::KmatrixMap
GenerateKmatrices(
  const u3::U3& sigma,
  const u3boson::U3BosonSpace& space,
  const double zero_threshold
)
  {
    // vcs::U3BosonSpace space(sigma,Nn_max);
    auto SMatrices = GenerateSmatrices(sigma,space);

    KmatrixMap KMatrix_map;

    if(sp3r::ModifySp3RBranching(sigma))
      {
        for(const auto& [omega,SMatrix] : SMatrices)
          {
            if(SMatrix.rows()==1)
              {
                double k_squared = SMatrix(0,0);

                if(fabs(k_squared)>zero_threshold)
                  {
                    // Since 16 was factored out of the definition of the Smatrix,
                    // we need to add it back it sqrt(16)^(Nn/2)=2^Nn
                    double factor = pow(2,int(omega.N()-sigma.N()));

                    double k=std::sqrt(k_squared);
                    KMatrix_map[omega][0]
                      =basis::OperatorBlock<double>::Identity(1,1)*k*factor;
                    KMatrix_map[omega][1]
                      =basis::OperatorBlock<double>::Identity(1,1)/k/factor;
                  }
                continue;
              }

            // Get Eigenvalues and eigenvectors
            Eigen::SelfAdjointEigenSolver<basis::OperatorBlock<double>> 
              eigen_system(SMatrix);

            basis::OperatorBlock<double> eigenvectors=eigen_system.eigenvectors();
            const basis::OperatorBlock<double>& eigenvalues=eigen_system.eigenvalues();
            assert(Eigen::ComputationInfo::Success == eigen_system.info());

            // determinant of eigenvectors comming out of eigen solver can have determinant +-1.
            // For K matrix calculation, need eigenvectors to have determinant 1.
            eigenvectors*=eigenvectors.determinant();
            
            std::vector<double>eigenvalue_vector;
            for(int i=0; i<eigenvalues.rows(); ++i)
              {
                // If sigma is a unitary Sp(3,R) irrep, then the eigenvalues of S
                // should always be >=0.
                if(eigenvalues(i,0)<0)
                  {
                    // If negative eigenvalue is ``non-zero", exit with failure.
                    if(fabs(eigenvalues(i,0))>zero_threshold)
                      {
                        fmt::print("negative eigenvalue for {}: {}\n",omega,eigenvalues(i,0));
                        exit (EXIT_FAILURE);
                      }
                  }

                eigenvalue_vector.push_back(eigenvalues(i,0));
              }
            // Get list of non-zero eigenvalues and the index of the corresponding eigenvector. 
            auto nonzero_eigs 
              = iter::filter(
                  [zero_threshold](const auto& eig){return eig.second > zero_threshold;},
                  iter::enumerate(eigenvalue_vector)
                );

            // upsilon_max corresponds to the number of non-zero eigenvalues of S
            int upsilon_max = std::distance(nonzero_eigs.begin(),nonzero_eigs.end());
            if(upsilon_max==0)
              continue;

            // Define K as UD and Kinv=DinvU^dagger.  Where D is a matix of dimension
            // u3boson_dimension x upsilon_max with matrix element (i,j) given by
            // sqrt(eig(j)) where i is the ith non-zero eigenvalue corresponding to the jth
            // eigenvector.
            //
            // K is a matrix of dimension u3boson_dimension x upsilon_max and
            // Kinv is a matrix of dimension upsilon_max x u3boson_dimension.
            // |sigma,upsilon,omega> = sum[n,rho] (K_inv)_{upsilon,[n,rho]}|sigma,n,rho,omega>
            auto& K = KMatrix_map[omega][0];
            auto& Kinv = KMatrix_map[omega][1];
            basis::OperatorBlock<double> D = basis::OperatorBlock<double>::Zero(SMatrix.cols(),upsilon_max);
            basis::OperatorBlock<double> Dinv = basis::OperatorBlock<double>::Zero(upsilon_max,SMatrix.cols());

            int i=0;
            for(const auto& [j,k_squared] : nonzero_eigs)
              {
                double k = std::sqrt(k_squared);
                D(j,i)=k;
                Dinv(i,j)=1/k;
                i++;
              }

            K=eigenvectors*D;
            Kinv=Dinv*eigenvectors.transpose();

            // Since 16 was factored out of the definition of the Smatrix,
            // we now need to add it back it sqrt(16)^(Nn/2)=2^Nn
            double factor = pow(2,int(omega.N()-sigma.N()));
            K*=factor;
            Kinv*=1./factor;

            assert(mcutils::IsZero(Kinv*K-basis::OperatorBlock<double>::Identity(upsilon_max,upsilon_max),1e-6));
            assert(mcutils::IsZero(K.transpose()*Kinv.transpose()-basis::OperatorBlock<double>::Identity(upsilon_max,upsilon_max),1e-6));
          }
      }
    else
      {
        for(const auto& [omega,SMatrix] : SMatrices)
          {
            // From pulling out a factor of 16 from the definition of the Smatrix
            // We now need to add it back it sqrt(16)^(Nn/2)=2^Nn
            double factor =pow(2,int(omega.N()-sigma.N())); 
            
            //calculate K matrix
            Eigen::SelfAdjointEigenSolver<basis::OperatorBlock<double>> eigen_system(SMatrix);
            KMatrix_map[omega][0]=eigen_system.operatorSqrt().cast<double>();
            KMatrix_map[omega][1]=KMatrix_map[omega][0].inverse();
            KMatrix_map[omega][1]*=1./factor;
            KMatrix_map[omega][0]*=factor;

            // /////////////////////////////////////////////////////////////////////////////
            // // For Testing
            // /////////////////////////////////////////////////////////////////////////////
            // const basis::OperatorBlock<double>& eigenvectors=eigen_system.eigenvectors();
            // const basis::OperatorBlock<double>& eigenvalues=eigen_system.eigenvalues();
            // assert(Eigen::ComputationInfo::Success == eigen_system.info());

            // basis::OperatorBlock<double> D=basis::OperatorBlock<double>::Zero(SMatrix.rows(),SMatrix.cols());
            // for(int i=0; i<eigenvalues.rows(); ++i)
            //   D(i,i)=std::sqrt(eigenvalues(i,0));

            // basis::OperatorBlock<double> K_test= eigenvectors*D*eigenvectors.transpose()*factor;
            // // if(!mcutils::IsZero(K_test-KMatrix_map[omega][0],1e-12))
            // //   std::cout<<K_test<<std::endl<<"---"<<std::endl<<KMatrix_map[omega][0]<<std::endl;

            // assert(mcutils::IsZero(K_test-KMatrix_map[omega][0],1e-6));
            /////////////////////////////////////////////////////////////////////////////
          }

      }
    return KMatrix_map;
  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}  //  namespace
