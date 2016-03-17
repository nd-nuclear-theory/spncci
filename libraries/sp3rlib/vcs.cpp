/****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy on 3/9/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py.

****************************************************************/

#include <cmath>  
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues>  
#include <string>

#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/sp3r.h"



namespace vcs{
  double BosonCreationRME(const u3::U3& np, const u3::U3& n)
	//  SU(3) Reduced matrix element of a^\dagger boson creation operator
  	{
  		double rme=0;
  		int n1=int(n.f1);
  		int n2=int(n.f2);
  		int n3=int(n.f3);

  		if((np.f1==(n.f1+2))&&(np.f2==n.f2)&&(np.f3==n.f3))
        rme=sqrt(
  				(n1+4)*(n1-n2+2)*(n1-n3+3)
  				/(2.*(n1-n2+3)*(n1-n3+4))
  				);

  		if ((np.f1==n.f1)&&(np.f2==(n.f2+2))&&(np.f3==n.f3))
  			rme=sqrt(
  				(n2+3)*(n1-n2)*(n2-n3+2)
  				/(2.*(n1-n2-1)*(n2-n3+3))
  				);
  		if ((np.f1==n.f1)&&(np.f2==n.f2)&&(np.f3==(n.f3+2)))
  			rme=sqrt(
  				(n3+2)*(n2-n3)*(n1-n3+1)/(2.*(n1-n3)*(n2-n3-1))
  				);
  		return rme;

  	}


  void GenerateKMatrices(const sp3r::Sp3RSpace& irrep, std::map<u3::U3,Eigen::MatrixXd>& K_matrix_map)
  {
    u3::U3 sigma=irrep.sigma();
    std::map<u3::U3,Eigen::MatrixXd> S_matrix_map;
    for(int i=0; i<irrep.size(); i++ )
    {
      
      // Generate S_matrix = K_matrix^2
      sp3r::U3Subspace u3_subspace_p=irrep.GetSubspace(i);
      u3::U3 omega_p=u3_subspace_p.GetSubspaceLabels();

      int dimension_p=u3_subspace_p.size();
      Eigen::MatrixXd S_matrix_p(dimension_p,dimension_p);
      Eigen::MatrixXd K_matrix_p(dimension_p,dimension_p);
      if (sigma==omega_p)
        {
          S_matrix_p(0,0)=1.0;
          K_matrix_p(0,0)=1.0;
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

              const sp3r::U3Subspace& u3_subspace=irrep.LookUpSubspace(omega);//GetSubspace().LookUpSubspaceIndex(omega);

              int dimension=u3_subspace.size();

              // OPTCHECK: Try doing mat-mat-mat by hand 
              Eigen::MatrixXd coef1_matrix(dimension_p,dimension);
              Eigen::MatrixXd coef2_matrix(dimension,dimension_p);
              
              // iterate over n1p,rho1p and n1p rho2p
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
                          double coef1=(
                            (Omega(n1p, omega_p)-Omega(n1,omega))
                            *u3::U(sigma.SU3(),n1.SU3(),omega_p.SU3(),u3::SU3(2,0),omega.SU3(),n1_rho1.tag,1,n1p.SU3(),1,n1p_rho1p.tag)
                            *BosonCreationRME(n1p,n1)
                            );
                          coef1_matrix(i1,j1)=coef1;
                        }

                      // Filling out coef1_matrix
                      for (int j2; j2<dimension; j2++)
                        {
                          MultiplicityTagged<u3::U3> n2_rho2=u3_subspace.GetStateLabels(j2);
                          u3::U3 n2(n2_rho2.irrep);
                          double coef2=(
                            (Omega(n2p, omega_p)-Omega(n2,omega))
                            *u3::U(sigma.SU3(),n2.SU3(),omega_p.SU3(),u3::SU3(2,0),omega.SU3(),n2_rho2.tag,1,n2p.SU3(),1,n2p_rho2p.tag)
                            *BosonCreationRME(n2p,n2)
                            );
                          coef2_matrix(i2,j2)=coef2;
                        }


                    }                  
                }
              S_matrix_p+=coef1_matrix*S_matrix_map[omega]*coef2_matrix;


            }
          //calculate K matrix 
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_system(S_matrix_p);
          K_matrix_p=eigen_system.operatorSqrt();
        }
      S_matrix_map[omega_p]=S_matrix_p;
      K_matrix_map[omega_p]=K_matrix_p;
    

      // else calculae S_matrix recuresively 
      //    generate list of possible omega
      //    for each omega, construct coefficient matrix 
      //    multiply coefficient matrix with S_matrix(omega)
      
    } // end for (i)  
  }





  // double SMatrix(const u3::U3& s, const u3::U3& wp, MultiplicityTagged<u3::U3>& nr1p, MultiplicityTagged<u3::U3>& nr2p)
  // {
  //   double smatrix=0.0;
  //   if (s==wp)
  //     return 1.0;

  //   u3::U3 n1p=nr1p.irrep;
  //   u3::U3 n2p=nr2p.irrep;
  //   std::vector<u3::U3> n1set;
  //   std::vector<u3::U3> n2set;
    
  //   if (n1p.f1-2>=n1p.f2)
  //     n1set.push_back(u3::U3(n1p.f1-2,n1p.f2,n1p.f3));
  //   if (n1p.f2-2>=n1p.f3)
  //     n1set.push_back(u3::U3(n1p.f1,n1p.f2-2,n1p.f3));
  //   if (n1p.f3-2>=0)
  //     n1set.push_back(u3::U3(n1p.f1,n1p.f2,n1p.f3-2));
    
  //   if (n2p.f1-2>=n2p.f2)
  //     n2set.push_back(u3::U3(n2p.f1-2,n2p.f2,n2p.f3));
  //   if (n2p.f2-2>=n2p.f3)
  //     n2set.push_back(u3::U3(n2p.f1,n2p.f2-2,n2p.f3));
  //   if (n2p.f3-2>=0)
  //     n2set.push_back(u3::U3(n2p.f1,n2p.f2,n2p.f3-2));

  //   MultiplicityTagged<u3::U3>::vector wset=KroneckerProduct(wp, u3::U3(0,0,-2));
  //   for (int i=0; i<wset.size(); i++)
  //     {
  //       u3::U3 w=wset[i].irrep;
  //       for(int j=0; j<n1set.size(); j++)
  //         {
  //           u3::U3 n1=n1set[j];
  //           int r1_max=OuterMultiplicity(s.SU3(),n1.SU3(),w.SU3());
  //           for(int r1=1; r1<=r1_max; r1++)
  //             {
  //               MultiplicityTagged<u3::U3> nr1(n1,r1);
  //               for(int k=0; k<n2set.size(); k++)
  //                 {
  //                   u3::U3 n2=n2set[k];
  //                   int r2_max=OuterMultiplicity(s.SU3(),n2.SU3(),w.SU3());
  //                   for (int r2=1; r2<=r2_max; r2++)
  //                     {
  //                       MultiplicityTagged<u3::U3> nr2(n2,r2);
  //                       smatrix+=(
  //                         (Omega(n1p, wp)-Omega(n1,w))
  //                         *u3::U(s.SU3(),n1.SU3(),wp.SU3(),u3::SU3(2,0),w.SU3(),r1,1,n1p.SU3(),1,nr1p.tag)
  //                         *BosonCreationRME(n1p,n1)
  //                         *u3::U(s.SU3(),n2.SU3(),wp.SU3(),u3::SU3(2,0),w.SU3(),r2,1,n2p.SU3(),1,nr2p.tag)
  //                         *BosonCreationRME(n2p,n2)
  //                         *SMatrix(s,w,nr1,nr2)
  //                         );
  //                     }


  //                 }

  //             }
  //         }

  //     }
  //   return 2.*smatrix/int(n1p.N());
  // }



}  //  namespace 
