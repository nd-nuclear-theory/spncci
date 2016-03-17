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

#include "sp3rlib/u3.h"
#include "sp3rlib/vcs.h"
#include "sp3rlib/sp3r.h"


int main(int argc, char **argv)
{
	u3::U3CoefInit();
	////////////////////////////////////////////////////////
	// Omega test
	////////////////////////////////////////////////////////
    u3::U3 w(HalfInt(9,2),HalfInt(1,2),HalfInt(1,2));
    u3::U3 nn(4,0,0);
    //  4.375.

	std::cout<< vcs::Omega(nn, w) << std::endl;

	////////////////////////////////////////////////////////
	// BosonCreationRME test
	////////////////////////////////////////////////////////
	// Checked against formula  
	// np.f1=n.f1+2
	u3::U3 np(6,4,2);
	u3::U3 n(4,4,2);
	std::cout << vcs::BosonCreationRME(np,n) << std::endl;
	// np.f2=n.f2+2
	np=u3::U3(6,4,2);
	n=u3::U3(6,2,2);
	std::cout << vcs::BosonCreationRME(np,n) << std::endl;
	// np.f3=n.f3+2
	np=u3::U3(6,4,2);
	n=u3::U3(6,4,0);
	std::cout << vcs::BosonCreationRME(np,n) << std::endl;

	////////////////////////////////////////////////////////
	//  Smatrix 
	////////////////////////////////////////////////////////
	u3::U3 s(20,13,10);
	sp3r::Sp3RSpace irrep(s,4);


	std::map<u3::U3,Eigen::MatrixXd> K_matrix_map;

  	vcs::GenerateKMatrices(irrep,K_matrix_map);

	// for(int i=0; i<irrep.size(); i++ )
 //    {
      
 //      // Generate S_matrix = K_matrix^2
 //      sp3r::U3Subspace u3_subspace_p=irrep.GetSubspace(i);
 //      u3::U3 omega_p=u3_subspace_p.GetSubspaceLabels();

 //      int dimension_p=u3_subspace_p.size();
 //      std::cout<<omega_p.Str()<<dimension_p<<std::endl;
 //      std::cout<<K_matrix_map[omega_p]<<std::endl;
 //  	}




  	std::map<u3::U3,Eigen::MatrixXd>::const_iterator it;
  	// std::cout<<K_matrix_map.begin()<<std::endl;
    for (auto it=K_matrix_map.begin(); it !=K_matrix_map.end(); ++it)
  		{
  			std::cout <<it->first.Str()<<it->second<<std::endl;
  		}



} // main 


