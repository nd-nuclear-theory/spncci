/****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy on 3/9/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py.

****************************************************************/
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"


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


	vcs::MatrixCache K_matrix_map;
	vcs::GenerateKMatrices(irrep,K_matrix_map);
  	// vcs::GenerateKMatrices(irrep,K_matrix_map);

	// for(int i=0; i<irrep.size(); i++ )
 //    {
      
 //      // Generate S_matrix = K_matrix^2
 //      sp3r::U3Subspace u3_subspace_p=irrep.GetSubspace(i);
 //      u3::U3 omega_p=u3_subspace_p.labels();

 //      int dimension_p=u3_subspace_p.size();
 //      std::cout<<omega_p.Str()<<dimension_p<<std::endl;
 //      std::cout<<K_matrix_map[omega_p]<<std::endl;
 //  	}




	vcs::MatrixCache::const_iterator it;
	// std::cout<<K_matrix_map.begin()<<std::endl;
  for (auto it=K_matrix_map.begin(); it !=K_matrix_map.end(); ++it)
		{
			std::cout <<it->first.Str()<<std::endl<<it->second<<std::endl;
		}


	std::cout<<std::endl<<"checking restricted irreps"<<std::endl;
	HalfInt Nsigma(10,2);
	bool restrict_sp3r_to_u3_branching=true;
	u3::U3 sigma(Nsigma, u3::SU3(0,1));
	sp3r::Sp3RSpace irrep_restricted(sigma,2,restrict_sp3r_to_u3_branching);
	std::cout<<irrep_restricted.DebugStr()<<std::endl;

	vcs::MatrixCache K_matrix_map_restricted;

 	vcs::GenerateKMatrices(irrep_restricted,K_matrix_map_restricted);
  for (auto it=K_matrix_map_restricted.begin(); it !=K_matrix_map_restricted.end(); ++it)
		{
			std::cout <<it->first.Str()<<std::endl<<it->second<<std::endl;
			std::cout<<"inverse "<<std::endl;
			std::cout<<it->second.inverse()<<std::endl;
		}




} // main 


