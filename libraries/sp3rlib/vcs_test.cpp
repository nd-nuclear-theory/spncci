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

  	//std::map<u3::U3,Eigen::MatrixXd>::const_interator it=K_matrix_map.begin();
  	//while (int it!=K_matrix_map.end())
  	//	{
  	//		std::cout<< it->first.Str() <<it->second<<std::endl;
  	//	}

      // s=U3State([20,13,10])
      // w=U3State([22,15,10])
      // n1=U3State([2,2,0])
      // r1=1
      // n2=U3State([4,0,0])
      // r2=1
      // smatrix(s,w,n1,r1,n1,r1) returns 
      //   1014.0
      // Smatrix(s,w,n1,r1,n2,r2) returns 
      //   -28.9827534924


} // main 


