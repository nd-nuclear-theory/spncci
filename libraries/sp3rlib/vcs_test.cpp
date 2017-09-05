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

	{
		std::cout<<std::endl<<"checking restricted irreps"<<std::endl;
		HalfInt Nsigma(10,2);
		bool restrict_sp3r_to_u3_branching=true;
		u3::U3 sigma(Nsigma, u3::SU3(2,0));

		sp3r::Sp3RSpace irrep_restricted;
		sp3r::ConstructRestrictedSp3RSpace(sigma,4,irrep_restricted);
		sp3r::Sp3RSpace irrep(sigma,4);

		std::cout<<irrep.DebugStr()<<std::endl;
		std::cout<<irrep_restricted.DebugStr()<<std::endl;

		vcs::MatrixCache K_matrix_map_restricted;
	 	vcs::GenerateKMatrices(irrep_restricted,K_matrix_map_restricted);

		vcs::MatrixCache K_matrix_map;
	 	vcs::GenerateKMatrices(irrep,K_matrix_map);


		for (auto it=K_matrix_map.begin(); it !=K_matrix_map.end(); ++it)
			{
				std::cout <<it->first.Str()<<std::endl<<it->second<<std::endl;
				if(K_matrix_map_restricted.count(it->first))
					std::cout<<"Restricted K "<<std::endl<<K_matrix_map_restricted[it->first]<<std::endl;

				std::cout<<"inverse "<<std::endl;
				std::cout<<it->second.inverse()<<std::endl;
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_system(it->second);
	      std::cout<<"eigenvalues"<<std::endl;
	      const Eigen::VectorXd& eigen_values=eigen_system.eigenvalues();
	      std::cout<<eigen_values<<std::endl;

			}

		u3::U3 w1(7,u3::SU3(2,1));
		u3::U3 w2(7,u3::SU3(4,0));
		u3::U3 n(2,u3::SU3(2,0));
		
		u3::U3 wp(9,u3::SU3(2,2));
		u3::U3 np(2,u3::SU3(2,0));
		std::cout<<w1.Str()<<"  "<<vcs::Omega(np,wp)-vcs::Omega(n,w1)<<std::endl;
		std::cout<<w2.Str()<<"  "<<vcs::Omega(np,wp)-vcs::Omega(n,w2)<<std::endl;


		// for (auto it=K_matrix_map_restricted.begin(); it !=K_matrix_map_restricted.end(); ++it)
		// 	{
		// 		std::cout <<it->first.Str()<<std::endl<<it->second<<std::endl;
		// 		std::cout<<"inverse "<<std::endl;
		// 		std::cout<<it->second.inverse()<<std::endl;
		// 		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_system(it->second);
	 //      std::cout<<"eigenvalues"<<std::endl;
	 //      const Eigen::VectorXd& eigen_values=eigen_system.eigenvalues();
	 //      std::cout<<eigen_values<<std::endl;
		// 	}
	}
} // main 


