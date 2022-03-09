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

namespace vcs{

  void GenerateSMatrices(const sp3r::Sp3RSpace& irrep, vcs::SMatrixCache& S_matrix_map, bool sp3r_u3_branch_restricted)
  {
    u3::U3 sigma=irrep.sigma();

    for(int i=0; i<irrep.size(); i++)
      {

        // Generate S_matrix = K_matrix^2
        sp3r::U3Subspace u3_subspace_p=irrep.GetSubspace(i);
        u3::U3 omega_p=u3_subspace_p.labels();

        int dimension_p=u3_subspace_p.size();
        vcs::SMatrixType S_matrix_p=vcs::SMatrixType::Zero(dimension_p,dimension_p);
        int upsilon_max=0;
        if (sigma==omega_p)
          {
            S_matrix_p(0,0)=1.0;
            upsilon_max=1;
          }
        // general case
        else
          {
            MultiplicityTagged<u3::U3>::vector omega_set=KroneckerProduct(omega_p, u3::U3(0,0,-2));
            // sum over omega
            for (int w=0; w<omega_set.size(); w++)
              {
                u3::U3 omega(omega_set[w].irrep);
                if (not irrep.ContainsSubspace(omega))
                  continue;

                if(not S_matrix_map.count(omega))
                  continue;

                const sp3r::U3Subspace& u3_subspace=irrep.LookUpSubspace(omega);//GetSubspace().LookUpSubspaceIndex(omega);

                int dimension=u3_subspace.size();
                // OPTCHECK: Try doing mat-mat-mat by hand
                vcs::SMatrixType coef1_matrix(dimension_p,dimension);
                vcs::SMatrixType coef2_matrix(dimension,dimension_p);

                // iterate over n1p,rho1p and n2p rho2p
                for (int i1=0; i1<dimension_p; i1++)
                  {

                    MultiplicityTagged<u3::U3> n1p_rho1p=u3_subspace_p.GetStateLabels(i1);
                    u3::U3 n1p(n1p_rho1p.irrep);
                    for (int i2=0; i2<dimension_p; i2++)
                      {
                        MultiplicityTagged<u3::U3> n2p_rho2p=u3_subspace_p.GetStateLabels(i2);
                        u3::U3 n2p(n2p_rho2p.irrep);

                        // Filling out coef1_matrix
                        for (int j1=0; j1<dimension; j1++)
                          {
                            MultiplicityTagged<u3::U3> n1_rho1=u3_subspace.GetStateLabels(j1);
                            u3::U3 n1(n1_rho1.irrep);
                            double long coef1;
                            if (u3::OuterMultiplicity(n1.SU3(), u3::SU3(2,0), n1p.SU3())>0)

                              coef1=(
                                     2./int(n1p.N())
                                     *(Omega(n1p, omega_p)-Omega(n1,omega))
                                     *U3BosonCreationRME(sigma,n1p_rho1p,omega_p,sigma, n1_rho1,omega)
                                     );

                            else
                              coef1=0.0;

                            coef1_matrix(i1,j1)=coef1;

                          }

                        // Filling out coef2_matrix
                        for (int j2=0; j2<dimension; j2++)
                          {
                            MultiplicityTagged<u3::U3> n2_rho2=u3_subspace.GetStateLabels(j2);
                            u3::U3 n2(n2_rho2.irrep);
                            if (u3::OuterMultiplicity(n2.SU3(),u3::SU3(2,0), n2p.SU3())>0)
                              coef2_matrix(j2,i2)=U3BosonCreationRME(sigma, n2p_rho2p, omega_p, sigma, n2_rho2, omega);
                            else
                              coef2_matrix(j2,i2)=0.0;
                          }
                      }
                  }
                double tolerance=1e-4;
                mcutils::ChopMatrix(coef1_matrix, tolerance);
                mcutils::ChopMatrix(coef2_matrix, tolerance);
                S_matrix_p+=coef1_matrix*S_matrix_map[omega]*coef2_matrix;
              }
          }//end else

        if(sp3r_u3_branch_restricted)
          {
            // Get Eigenvalues and eigenvectors
            Eigen::SelfAdjointEigenSolver<vcs::SMatrixType> eigen_system(S_matrix_p);
            const vcs::SMatrixType& eigenvectors=eigen_system.eigenvectors();
            const vcs::SMatrixType& eigenvalues=eigen_system.eigenvalues();

            // sqrt(sum(matrix elements)^2)
            double sum=0;
            for(int i=0; i<eigenvalues.size(); ++i)
              sum+=eigenvalues(i);

            double norm_factor=sum/eigenvalues.size();

            vcs::SMatrixType eigenvalues_matrix=vcs::SMatrixType::Zero(eigenvalues.size(),eigenvalues.size());
            if(norm_factor>1e-2)
              for(int i=0; i<eigenvalues.size(); ++i)
                {
                  double k2=eigenvalues(i);
                  if(fabs(k2/norm_factor)>1e-6)
                    eigenvalues_matrix(i,i)=k2;
                }
            vcs::SMatrixType S_matrix_p2=eigenvectors*eigenvalues_matrix*eigenvectors.transpose();

            // mcutils::ChopMatrix(S_matrix_p, 1e-4);
            // if(not mcutils::IsZero(S_matrix_p2-S_matrix_p,1e-6))
              // std::cout<<"S matrix diff"<<omega_p.Str()<<std::endl<<S_matrix_p2<<std::endl<<S_matrix_p<<std::endl;

            mcutils::ChopMatrix(S_matrix_p, 1e-6);
            S_matrix_map[omega_p]=S_matrix_p2;
          }
        else
          S_matrix_map[omega_p]=S_matrix_p;

      } // end for (i)
  }

  // K matrix computed by taking the built-in Eigen operator square-root of S=KK^dagger
  // K=UDU^dagger where U is the eigenvectors of S and D is the diagonal matrix with
  // diagonal values given by the square root of the eigenvalues of S.
  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map)
  {
    vcs::SMatrixCache S_matrix_map;
    bool sp3r_u3_branch_restricted=false;
    vcs::GenerateSMatrices(irrep,S_matrix_map,false);
    for(auto it=S_matrix_map.begin(); it!=S_matrix_map.end(); ++it)
      {
        //calculate K matrix
        Eigen::SelfAdjointEigenSolver<vcs::SMatrixType> eigen_system(it->second);
        K_matrix_map[it->first]=eigen_system.operatorSqrt().cast<double>();
      }
  }


  // K matrix obtained by solving for eigenvalues Lambda and eigenvectors U of KK^dagger
  // as descripted in D. J. Rowe, A. E. McCoy and M. A. Caprio, Phys. Scripta 91 (2016) 0330003.
  // K(i,j)=Sqrt(lambda_i)U(i,j)
  // Kinv(j,i)=Sqrt(lambda_i)^(-1)U(i,j).transpose
  // Note K compute here differs from K computed in function above and is not symmetric
  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, vcs::MatrixCache& K_matrix_map, vcs::MatrixCache& Kinv_matrix_map)
  {
    vcs::SMatrixCache S_matrix_map;
    bool sp3r_u3_branch_restricted=true;
    vcs::GenerateSMatrices(irrep,S_matrix_map,sp3r_u3_branch_restricted);
    for(auto it=S_matrix_map.begin(); it!=S_matrix_map.end(); ++it)
      {
        // Get Eigenvalues and eigenvectors
        Eigen::SelfAdjointEigenSolver<vcs::SMatrixType> eigen_system(it->second);
        const vcs::SMatrixType& eigenvectors=eigen_system.eigenvectors();
        const vcs::SMatrixType& eigenvalues=eigen_system.eigenvalues();

        // sqrt(sum(matrix elements)^2)
        double sum=0;
        for(int i=0; i<eigenvalues.size(); ++i)
          sum+=fabs(eigenvalues(i));

        double norm_factor=sum/eigenvalues.size();

        if(fabs(norm_factor)<1e-2)
          continue;

        // Loop through eigenvalues and identify which eigenvalues are non-zero
        std::vector<int> non_zero_eigen_positions;
        for(int i=0; i<eigenvalues.size(); ++i)
        {
          // std::cout<<eigenvalues(i)<<"  "<<norm_factor<<"  "<<eigenvalues(i)/norm_factor<<std::endl;
          if(fabs(eigenvalues(i)/norm_factor)>1e-6)
          {
            non_zero_eigen_positions.push_back(i);
          }
        }
        if(non_zero_eigen_positions.size()==0)
          continue;

        // Initialize K and Kinv
        int rows=non_zero_eigen_positions.size();
        int cols=eigenvalues.size();

        vcs::SMatrixType K(rows,cols);
        vcs::SMatrixType Kinv(cols,rows);


        // std::cout<<"Eigenvalues "<<non_zero_eigen_positions.size()<<std::endl<<eigenvalues<<std::endl;
        // Construct K and Kinv from non-zero eigenvalues and corresponding eigenvectors
        for(int i=0; i<non_zero_eigen_positions.size(); ++i)
          {
            int index=non_zero_eigen_positions[i];
            double k=sqrt(eigenvalues(index));
            K.row(i)=k*eigenvectors.row(index);
            Kinv.col(i)=1/k*eigenvectors.row(index).transpose();
          }

        // Eigen::MatrixXd K=K1;
        // Eigen::MatrixXd Kinv=K1inv;

        K_matrix_map[it->first]=K.cast<double>();
        Kinv_matrix_map[it->first]=Kinv.cast<double>();

      }
  }


// using KmatrixMap = std::unordered_map<u3::U3, std::array<basis::OperatorBlock<double>, 2>>;
using SmatrixMap = std::unordered_map<u3::U3, basis::OperatorBlock<double>>;

SmatrixMap
GenerateSmatrices(
  const u3::U3& sigma,
  const std::map<u3::U3, MultiplicityTagged<u3::U3>::vector>& u3_subspaces
  )
  {
    SmatrixMap Smatrix_map;
    Smatrix_map[sigma]=basis::OperatorBlock<double>::Identity(1,1);
    
    for(const auto& [omegap,state_vector_p] : u3_subspaces)
      {
        if(omegap==sigma)
          continue;

        auto Nnp = int(omegap.N()-sigma.N());
        auto omega_set = u3::KroneckerProduct(omegap,{-2,{0,2}});
        
        std::vector<MultiplicityTagged<u3::U3>> boson_states_p;
        for(const auto&[n,rho_max] : state_vector_p)
          for(int rho : iter::range(1,rho_max+1))
            boson_states_p.push_back({n,rho});

        // Initialize S matrix 
        int num_states_p = boson_states_p.size();
        basis::OperatorBlock<double>& S_target = Smatrix_map[omegap];
        S_target = Eigen::MatrixXd::Zero(num_states_p,num_states_p);

        for(auto& [omega,dummy] : omega_set)
          {
            if (not u3_subspaces.count(omega))
              continue;
            
            const auto& S_source = Smatrix_map.at(omega);
            const auto& state_vector = u3_subspaces.at(omega);

            std::vector<MultiplicityTagged<u3::U3>> boson_states;
            for(const auto&[n,rho_max] : state_vector)
              for(int rho : iter::range(1,rho_max+1))
                boson_states.push_back({n,rho});
            
            int num_states = boson_states.size();            
            Eigen::MatrixXd coef1_matrix(num_states_p,num_states);
            Eigen::MatrixXd coef2_matrix(num_states,num_states_p);
            for(int index1p=0; index1p<num_states_p; index1p++)
              for(int index1=0; index1<num_states; index1++)
                {
                  const auto& [n1p,rho1p] = boson_states_p[index1p];
                  const auto& [n1,rho1] = boson_states[index1];
                       
                  coef1_matrix(index1p,index1)
                    =(Omega(n1p, omegap)-Omega(n1,omega))
                      *U3BosonCreationRME(sigma,n1p,rho1p,omegap,sigma,n1,rho1,omega);
                  coef2_matrix(index1,index1p)
                    =U3BosonCreationRME(sigma,n1p,rho1p,omegap,sigma,n1,rho1,omega);
                }
            S_target+=coef1_matrix*S_source*coef2_matrix/(8*Nnp);
          }
      }
  
    return Smatrix_map;
  }


KmatrixMap
 GenerateKMatrices(
    const u3::U3& sigma,
    const std::map<u3::U3, MultiplicityTagged<u3::U3>::vector>& u3_subspaces
  )
 {
  
  KmatrixMap Kmatrix_map;
  Kmatrix_map[sigma][0]=basis::OperatorBlock<double>::Identity(1,1);
  Kmatrix_map[sigma][1]=basis::OperatorBlock<double>::Identity(1,1);
  
  auto Smatrix_map = GenerateSmatrices(sigma,u3_subspaces);
  // std::cout<<" Smatrices computed"<<std::endl;
  if (sp3r::ModifySp3RBranching(sigma))
    {
      fmt::print("Code for hadling A<6 irreps has not implemented in this function.\n");
      std::exit(EXIT_FAILURE);
    }
  else
    {
      for(const auto& [omega,Smatrix] :  Smatrix_map)
        {
          // From pulling out a factor of 16 from the definition of the Smatrix
          // We now need to add it back it sqrt(16)^(Nn/2)=2^Nn
          double factor =pow(2,int(omega.N()-sigma.N())); 
          
          //calculate K matrix
          Eigen::SelfAdjointEigenSolver<basis::OperatorBlock<double>> eigen_system(Smatrix);
          Kmatrix_map[omega][0]=factor*eigen_system.operatorSqrt().cast<double>();
          Kmatrix_map[omega][1]=Kmatrix_map[omega][0].inverse()/factor;
        }
    }

  return Kmatrix_map;
 }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Based on U3Boson basis construction
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unordered_map<u3::U3,basis::OperatorBlock<double>>
GetSMatrics(const u3::U3& sigma, const vcs::U3BosonSpace& space)
  {
    // vcs::U3BosonSpace space(sigma,Nn_max);
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

                  //testing
                  if(DeltaOmega<0)
                    continue;

                  for(int rho_bra=1; rho_bra<=rho_max_bra; ++rho_bra)
                    for(int rho_ket=1; rho_ket<=rho_max_ket; ++rho_ket)
                    {
                      int row = bra_subspace.GetStateOffset(bra_state_index,rho_bra);
                      int col = ket_subspace.GetStateOffset(ket_state_index,rho_ket);
                      double boson_rme = vcs::U3BosonCreationRME(sigma,n_bra,rho_bra,omega_bra,sigma,n_ket,rho_ket,omega_ket);
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
GetKMatrices(
  const u3::U3& sigma,
  const vcs::U3BosonSpace& space,
  const double zero_threshold
)
  {
    // vcs::U3BosonSpace space(sigma,Nn_max);
    auto SMatrices = GetSMatrics(sigma,space);

    KmatrixMap KMatrix_map;

    if(sp3r::ModifySp3RBranching(sigma))
      {
        for(const auto& [omega,SMatrix] : SMatrices)
          {
            if(SMatrix.rows()==1)
              {
                double k_squared = SMatrix(0,0);
                if(fabs(k_squared)>zero_threshold && k_squared>0)
                  {
                    double k=std::sqrt(k_squared);
                    KMatrix_map[omega][0]=basis::OperatorBlock<double>::Identity(1,1)*k;
                    KMatrix_map[omega][0]=basis::OperatorBlock<double>::Identity(1,1)/k;
                  }
                continue;
              }

            // Get Eigenvalues and eigenvectors
            Eigen::SelfAdjointEigenSolver<basis::OperatorBlock<double>> 
              eigen_system(SMatrix);

            const basis::OperatorBlock<double>& eigenvectors=eigen_system.eigenvectors();
            const basis::OperatorBlock<double>& eigenvalues=eigen_system.eigenvalues();
            assert(Eigen::ComputationInfo::Success == eigen_system.info());
            
            std::vector<double>eigenvalue_vector;
            for(int i=0; i<eigenvalues.rows(); ++i)
              {
                // Negative eigenvalues can appear if the irrep is non-unitary.  For now, we 
                // just eliminate those states from the basis TODO: determine proper
                // branching rule.
                if(eigenvalues(i,0)<0)
                  {
                    fmt::print("negative eigenvalue for {}: {}\n",omega,eigenvalues(i,0));
                    // eigenvalue_vector.push_back(0);
                  }
                else {eigenvalue_vector.push_back(eigenvalues(i,0));}
              }
            // Get list of non-zero eigenvalues and the index of the corresponding eigenvector. 
            auto nonzero_eigs 
              = iter::filter([zero_threshold](const auto& eig){return eig.second > zero_threshold;}, iter::enumerate(eigenvalue_vector));

            // upsilon_max corresponds to the number of non-zero eigenvalues of S
            int upsilon_max = std::distance(nonzero_eigs.begin(),nonzero_eigs.end());
            if(upsilon_max==0)
              continue;
            // K matrix defined such that 
            // |sigma,upsilon,omega> 
            //    = sum[n,rho] (K_inv)_{upsilon,[n,rho]}|sigma,n,rho,omega>
            auto& K = KMatrix_map[omega][0];
            auto& K_inv = KMatrix_map[omega][1];
            K.resize(SMatrix.cols(),upsilon_max);
            K_inv.resize(upsilon_max,SMatrix.cols());
            int i=0;
            for(const auto& [j,k_squared] : nonzero_eigs)
              {
                double k = std::sqrt(k_squared);
                K.col(i) = k*eigenvectors.col(j);
                K_inv.row(i) = 1/k*eigenvectors.col(j).transpose(); 
                i++;
              }

            // From pulling out a factor of 16 from the definition of the Smatrix
            // We now need to add it back it sqrt(16)^(Nn/2)=2^Nn
            double factor =pow(2,int(omega.N()-sigma.N()));
            K=K*factor;
            K_inv=K_inv/factor;
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
            KMatrix_map[omega][0]=factor*eigen_system.operatorSqrt().cast<double>();
            KMatrix_map[omega][1]=KMatrix_map[omega][0].inverse()/factor;
          }

      }
    return KMatrix_map;
  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}  //  namespace
