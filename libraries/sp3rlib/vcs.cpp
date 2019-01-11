/****************************************************************
  vcs.cpp                       

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "sp3rlib/vcs.h"

#include <omp.h>
#include "fmt/format.h"
#include <eigen3/Eigen/Eigenvalues>  
#include "sp3rlib/u3coef.h"   
#include "mcutils/eigen.h"

namespace vcs{

  double BosonCreationRME(const u3::U3& np, const u3::U3& n)
  //  SU(3) Reduced matrix element of a^\dagger boson creation operator
  {
    double rme=0;
    int n1=int(n.f1());
    int n2=int(n.f2());
    int n3=int(n.f3());

    if((np.f1()==(n.f1()+2))&&(np.f2()==n.f2())&&(np.f3()==n.f3()))
      rme=sqrt(
               (n1+4)*(n1-n2+2)*(n1-n3+3)
               /(2.*(n1-n2+3)*(n1-n3+4))
               );

    if ((np.f1()==n.f1())&&(np.f2()==(n.f2()+2))&&(np.f3()==n.f3()))
      rme=sqrt(
               (n2+3)*(n1-n2)*(n2-n3+2)
               /(2.*(n1-n2-1)*(n2-n3+3))
               );
    if ((np.f1()==n.f1())&&(np.f2()==n.f2())&&(np.f3()==(n.f3()+2)))
      rme=sqrt(
               (n3+2)*(n2-n3)*(n1-n3+1)/(2.*(n1-n3)*(n2-n3-1))
               );
    return rme;

  }

  double U3BosonCreationRME(
          const u3::U3& sigmap, const MultiplicityTagged<u3::U3>np_rhop,  const u3::U3& omegap,
          const u3::U3& sigma, const MultiplicityTagged<u3::U3> n_rho, const u3::U3& omega
          )
  {
    int rho0_max=u3::OuterMultiplicity(omega.SU3(),u3::SU3(2,0),omegap.SU3());
    int rhon_max=u3::OuterMultiplicity(n_rho.irrep.SU3(),u3::SU3(2,0),np_rhop.irrep.SU3());
    if ((sigmap==sigma)&&(rho0_max>0)&&(rhon_max>0))
      {  
        
        double rme=ParitySign(u3::ConjugationGrade(omegap)+u3::ConjugationGrade(omega))
                    *u3::U(u3::SU3(2,0),n_rho.irrep.SU3(),omegap.SU3(),sigma.SU3(),np_rhop.irrep.SU3(),1,np_rhop.tag,omega.SU3(),n_rho.tag,1)
                    // *u3::U(sigma.SU3(), n_rho.irrep.SU3(), omegap.SU3(), u3::SU3(2,0), omega.SU3(),n_rho.tag,1,np_rhop.irrep.SU3(),1,np_rhop.tag)
                    *BosonCreationRME(np_rhop.irrep,n_rho.irrep);
        return rme;
      }
    else
      return 0.0;
  }


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

        if(fabs(norm_factor<1e-2))
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




}  //  namespace 
