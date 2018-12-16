/****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  Created by Anna E. McCoy on 3/9/16.   
 
  3/9/16 (aem): Created based on prototype vcs.py.
	10/4/17 (aem) : Add tests for A<6 Kmatrices
****************************************************************/
#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"
#include "mcutils/eigen.h"

void CreateU3BosonMatrix(
	const u3::U3& omegap, const u3::U3& omega,
	const sp3r::Sp3RSpace& irrep,
	Eigen::MatrixXd& boson_matrix
	)
  {
    const u3::U3& sigma(irrep.sigma());
    const sp3r::U3Subspace& omegap_subspace=irrep.LookUpSubspace(omegap);
    const sp3r::U3Subspace& omega_subspace=irrep.LookUpSubspace(omega);
    
	  boson_matrix.resize(omegap_subspace.size(),omega_subspace.size());
    for(int bra_index=0; bra_index<omegap_subspace.size(); ++bra_index)
    	for(int ket_index=0; ket_index<omega_subspace.size(); ++ket_index)
    		{
    			const MultiplicityTagged<u3::U3>& np_rhop(omegap_subspace.GetStateLabels(bra_index));
    			const MultiplicityTagged<u3::U3>& n_rho(omega_subspace.GetStateLabels(ket_index));
    			boson_matrix(bra_index,ket_index)=vcs::U3BosonCreationRME(sigma,np_rhop,omegap,sigma,n_rho,omega);
    		}    
 	}

void CreateOmegaMatrix(const u3::U3& omegap, const u3::U3& omega,
	const sp3r::Sp3RSpace& irrep,
	Eigen::MatrixXd& Omega_matrix
)
	{
    const sp3r::U3Subspace& omegap_subspace=irrep.LookUpSubspace(omegap);
    const sp3r::U3Subspace& omega_subspace=irrep.LookUpSubspace(omega);
    
	  Omega_matrix.resize(omegap_subspace.size(),omega_subspace.size());
    for(int bra_index=0; bra_index<omegap_subspace.size(); ++bra_index)
    	for(int ket_index=0; ket_index<omega_subspace.size(); ++ket_index)
    		{
    			const u3::U3& np(omegap_subspace.GetStateLabels(bra_index).irrep);
    			const u3::U3& n(omega_subspace.GetStateLabels(ket_index).irrep);
    			Omega_matrix(bra_index,ket_index)=vcs::Omega(np,omegap)-vcs::Omega(n,omega);
    		}

	}




int main(int argc, char **argv)
{
	u3::U3CoefInit();

if(false)
{
	Eigen::MatrixXd A(2,2);
	A(0,0)=1;
	A(0,1)=2;
	A(1,0)=4;
	A(1,1)=6;
	std::cout<<"Matrix A\n"<<A<<std::endl;
	Eigen::MatrixXd B=A.col(0)*A.row(1);
	Eigen::MatrixXd C=A.col(1)*A.row(1);
	std::cout<<B<<std::endl<<std::endl<<C<<std::endl;
	double a=(B.array()*C.array()).sum();
	std::cout<<a<<std::endl;
}	
if(false)
{
	Eigen::MatrixXd K(2,2);
	K(0,0)=1;
	K(0,1)=1;
	K(1,0)=1;
	K(1,1)=-1;

	std::cout<<"inverse "<<std::endl;
	Eigen::MatrixXd id=Eigen::MatrixXd::Identity(K.rows(),K.rows());
	Eigen::MatrixXd inv=K.colPivHouseholderQr().solve(id);
	// Eigen::MatrixXd inv2=mcutils::RightInverse(K);
	Eigen::MatrixXd inv3=K.inverse();
	// std::cout<<inv<<std::endl<<std::endl<<inv2<<std::endl<<std::endl<<inv3<<std::endl<<std::endl;
	std::cout<<K*inv<<std::endl;
}
if(false)
{
	Eigen::MatrixXd KK(2,3);
	KK(0,0)=1;
	KK(0,1)=2;
	KK(0,2)=3;
	KK(1,0)=3;
	KK(1,1)=5;
	KK(1,2)=4;

	std::cout<<KK.transpose()*KK<<std::endl;


	std::cout<<"inverse "<<std::endl;
	Eigen::MatrixXd id=Eigen::MatrixXd::Identity(KK.rows(),KK.rows());
	Eigen::MatrixXd inv=KK.colPivHouseholderQr().solve(id);
	// Eigen::MatrixXd inv2=mcutils::RightInverse(KK);
	// std::cout<<inv<<std::endl<<inv2<<std::endl;
}
if(false)
{

	Eigen::MatrixXd K(2,2);
	K(0,0)=1.55708;
	K(0,1)=-1.0545;
	K(1,0)=-1.0545;
	K(1,1)=1.37135;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_system(K);
  const Eigen::MatrixXd& eigenvectors=eigen_system.eigenvectors();
  const Eigen::MatrixXd& eigenvalues=eigen_system.eigenvalues();

    double s=eigenvalues(0);
    double k=sqrt(eigenvalues(0));
    std::cout<<"eigenvalues"<<eigenvalues<<std::endl;
    std::cout<<"eigenvectors"<<std::endl;
    std::cout<<eigenvectors<<std::endl;
    std::cout<<eigenvectors.transpose()*K*eigenvectors<<std::endl<<std::endl;
    
    std::cout<<eigen_system.operatorSqrt()<<std::endl;
		Eigen::MatrixXd Sqrt=Eigen::MatrixXd::Zero(eigenvectors.rows(),eigenvectors.cols());
		Eigen::MatrixXd Inverse=Eigen::MatrixXd::Zero(eigenvectors.rows(),eigenvectors.cols());
    for(int i=0; i<eigenvalues.size(); i++)
      {
        double k=sqrt(eigenvalues(i));
        Sqrt+=eigenvectors.row(i).transpose()*k*eigenvectors.row(i);
        // std::cout<<k*eigenvectors.row(i)<<std::endl;
        Inverse+=eigenvectors.row(i).transpose()*1/k*eigenvectors.row(i);
        // Kinv.col(i)=1/k*eigenvectors.col(i).transpose();
      }

    std::cout<<Sqrt<<std::endl<<std::endl;;
    std::cout<<Inverse<<std::endl;
    std::cout<<K.inverse()<<std::endl;
    

		// std::cout<<eigen_system.operatorSqrt()*eigen_system.operatorSqrt().transpose()<<std::endl;
}


if(false)
{
	u3::U3 s(20,13,10);
	sp3r::Sp3RSpace irrep(s,4);

	vcs::MatrixCache K_matrix_map,K_matrix_map2, Kinv_matrix_map;
	vcs::GenerateKMatrices(irrep,K_matrix_map);
	vcs::GenerateKMatrices(irrep,K_matrix_map2, Kinv_matrix_map);

	for(auto it=K_matrix_map.begin(); it!=K_matrix_map.end(); ++it)
		{
			const Eigen::MatrixXd& K=it->second;
			const u3::U3& omega=it->first;
			
			Eigen::MatrixXd Kinv=Kinv_matrix_map[it->first];
			std::cout<<"K matrix for "<<it->first.Str()<<std::endl;
			std::cout<<K<<std::endl<<K_matrix_map2[it->first]<<std::endl<<std::endl;
			std::cout<<K.inverse()<<std::endl<<Kinv<<std::endl<<std::endl;
		}
}

if(true)
{
	u3::U3 s(11,1,1);
	sp3r::Sp3RSpace irrep(s,20);
	// sp3r::Sp3RSpace irrep(s,25);
	vcs::MatrixCache K_matrix_map, Kinv_matrix_map;
	vcs::GenerateKMatrices(irrep,K_matrix_map, Kinv_matrix_map);

	for(auto it=K_matrix_map.begin(); it!=K_matrix_map.end(); ++it)
		{
			const Eigen::MatrixXd& K=it->second;
			const u3::U3& omega=it->first;
			
			if(not K_matrix_map.count(it->first))
				continue;
			Eigen::MatrixXd Kinv=Kinv_matrix_map.at(it->first);
			std::cout<<"-----------"<<std::endl;
			std::cout<<"K matrix for "<<it->first.Str()<<std::endl;
			
			std::cout<<K_matrix_map.at(it->first)<<std::endl;
			std::cout<<Kinv<<std::endl<<std::endl;
			std::cout<<K*Kinv<<std::endl;
			std::cout<<"-----------"<<std::endl;
		}
}



	// u3::U3 sigma(16,u3::SU3(4,0));
	// sp3r::Sp3RSpace irrep(sigma,4);

	// for(int i=0; i<irrep.size(); ++i)
	// 	for(int j=0; j<irrep.size(); ++j)
	// 	{
	// 		Eigen::MatrixXd boson_matrix;
	// 		const u3::U3& omegap(irrep.GetSubspace(j).labels());
	// 		const u3::U3& omega(irrep.GetSubspace(i).labels());
	// 		CreateU3BosonMatrix(omegap,omega,irrep,boson_matrix);
	// 		std::cout<<omegap.Str()<<"  "<<omega.Str()<<std::endl;
	// 		std::cout<<boson_matrix<<std::endl<<std::endl;
	// 	}

// vcs::MatrixCache K_matrix_map;
// vcs::MatrixCache K_matrix_map2;
// // GenerateKMatrices2(irrep,K_matrix_map2);
// vcs::GenerateKMatrices(irrep, K_matrix_map);

// for(auto it=K_matrix_map.begin(); it!=K_matrix_map.end(); ++it)
// 	{
// 		std::cout<<it->first.Str()<<std::endl;
// 		// std::cout<<it->second<<std::end;
// 		// std::cout<<K_matrix_map2[it->first]<<std::end;
// 	}

if(false)
{
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
 //      sp3r::U3Subspace u3_subspacep=irrep.GetSubspace(i);
 //      u3::U3 omegap=u3_subspacep.labels();

 //      int dimensionp=u3_subspacep.size();
 //      std::cout<<omegap.Str()<<dimensionp<<std::endl;
 //      std::cout<<K_matrix_map[omegap]<<std::endl;
 //  	}

}
	// HalfInt Nsigma(10,2);
	// bool restrict_sp3r_to_u3_branching=true;
	// u3::U3 sigma(Nsigma, u3::SU3(2,0));

	// sp3r::Sp3RSpace irrep(sigma,10);

	// vcs::MatrixCache K_matrix_map;
 // 	vcs::GenerateKMatrices(irrep,K_matrix_map);

	// for (auto it=K_matrix_map.begin(); it !=K_matrix_map.end(); ++it)
	// 	{
	// 		std::cout <<it->first.Str()<<std::endl<<it->second<<std::endl;
	// 		// Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_system(it->second);
 //   //    std::cout<<"eigenvalues"<<std::endl;
 //   //    const Eigen::VectorXd& eigen_values=eigen_system.eigenvalues();
 //      Eigen::MatrixXd Kmatrix=it->second;
 //      // std::cout<<eigen_values<<std::endl;
 //      // int num_eigs=0;
 //      // for(int i=0; i<eigen_values.size(); ++i)
 //      // 	if(eigen_values(i)>1e-6)
 //      // 		{
 //      // 			num_eigs++;
 //      // 		}
 //      // if(not num_eigs)
 //      // 	{
 //      // 		std::cout<<"remove "<<it->first.Str()<<std::endl;
 //      // 		continue;
 //      // 	}

 //      // Eigen::MatrixXd Kmatrix2=Kmatrix.block(0,0,Kmatrix.rows(),num_eigs);
 //      // std::cout<<Kmatrix2<<std::endl;
	// 		std::cout<<"inverse "<<std::endl;
	// 		Eigen::MatrixXd identity=Eigen::MatrixXd::Identity(Kmatrix.cols(),Kmatrix.cols());
	// 		Eigen::MatrixXd inverse=Kmatrix.transpose().colPivHouseholderQr().solve(identity).transpose();
	// 		std::cout<<inverse<<std::endl<<std::endl;
	// 		// std::cout<<it->second.inverse()<<std::endl;
 //     }
			

		// u3::U3 w1(7,u3::SU3(2,1));
		// u3::U3 w2(7,u3::SU3(4,0));
		// u3::U3 n(2,u3::SU3(2,0));
		
		// u3::U3 wp(9,u3::SU3(2,2));
		// u3::U3 np(2,u3::SU3(2,0));
		// std::cout<<w1.Str()<<"  "<<vcs::Omega(np,wp)-vcs::Omega(n,w1)<<std::endl;
		// std::cout<<w2.Str()<<"  "<<vcs::Omega(np,wp)-vcs::Omega(n,w2)<<std::endl;


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


} // main 


