/****************************************************************
  vcs.h                       

  Define vector coherent state methods for Sp(3,R).

  Anna E. McCoy
  University of Notre Dame

  SPDX-License-Identifier: MIT 
 
   3/9/16 (aem) : Created.
  10/4/17 (aem) : Add tests for A<6 Kmatrices
****************************************************************/
#include "sp3rlib/u3coef.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/vcs.h"
#include "mcutils/eigen.h"
#include "cppitertools/itertools.hpp"

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
namespace vcs
{
	using Matrix =  basis::OperatorBlock<double>;
}

int main(int argc, char **argv)
{
	u3::U3CoefInit(39);

/////////////////////////////////////////////////////////////////////////////
// Test Omega function
std::map<std::tuple<u3::U3,u3::U3>,double> 
Omega_validation_table ={
		{{u3::U3(0,0,0),u3::U3(5,4,3)},27},
		{{u3::U3(2,0,0),u3::U3(5,5,4)},30},
		{{u3::U3(2,0,0),u3::U3(6,5,3)},34},
		{{u3::U3(2,0,0),u3::U3(6,4,4)},32},
		{{u3::U3(2,0,0),u3::U3(7,4,3)},37},
		{{u3::U3(2,2,0),u3::U3(6,6,4)},39},
		{{u3::U3(2,2,0),u3::U3(6,5,5)},37},
		{{u3::U3(2,2,0),u3::U3(7,6,3)},44},
		{{u3::U3(2,2,0),u3::U3(7,5,4)},41},
		{{u3::U3(4,0,0),u3::U3(7,5,4)},38},
		{{u3::U3(4,0,0),u3::U3(8,5,3)},44},
		{{u3::U3(4,0,0),u3::U3(8,4,4)},42},
		{{u3::U3(4,0,0),u3::U3(9,4,3)},49},
		{{u3::U3(0,0,0),u3::U3(4.5_hi,2.5_hi,1.5_hi)},17.375},
		{{u3::U3(2,0,0),u3::U3(4.5_hi,4.5_hi,1.5_hi)},20.375},
		{{u3::U3(2,0,0),u3::U3(4.5_hi,3.5_hi,2.5_hi)},17.375},
		{{u3::U3(2,0,0),u3::U3(5.5_hi,3.5_hi,1.5_hi)},22.375},
		{{u3::U3(2,0,0),u3::U3(5.5_hi,2.5_hi,2.5_hi)},20.375},
		{{u3::U3(2,0,0),u3::U3(6.5_hi,2.5_hi,1.5_hi)},26.375},
		{{u3::U3(2,2,0),u3::U3(4.5_hi,4.5_hi,3.5_hi)},20.375},
		{{u3::U3(2,2,0),u3::U3(5.5_hi,4.5_hi,2.5_hi)},24.375},
		{{u3::U3(4,0,0),u3::U3(5.5_hi,4.5_hi,2.5_hi)},21.375},
		{{u3::U3(2,2,0),u3::U3(5.5_hi,3.5_hi,3.5_hi)},22.375},
		{{u3::U3(2,2,0),u3::U3(6.5_hi,4.5_hi,1.5_hi)},30.375},
		{{u3::U3(4,0,0),u3::U3(6.5_hi,4.5_hi,1.5_hi)},27.375},
		{{u3::U3(2,2,0),u3::U3(6.5_hi,3.5_hi,2.5_hi)},27.375},
		{{u3::U3(4,0,0),u3::U3(6.5_hi,3.5_hi,2.5_hi)},24.375},
		{{u3::U3(4,0,0),u3::U3(7.5_hi,3.5_hi,1.5_hi)},31.375},
		{{u3::U3(4,0,0),u3::U3(7.5_hi,2.5_hi,2.5_hi)},29.375},
		{{u3::U3(4,0,0),u3::U3(8.5_hi,2.5_hi,1.5_hi)},37.375}
	};

for(const auto& [u3_pair, value] : Omega_validation_table)
	{
		const auto& [n,omega] = u3_pair;
		assert(Omega_validation_table[u3_pair] == vcs::Omega(n,omega));
	}

/////////////////////////////////////////////////////////////////////////////
std::map<u3::U3,std::map<u3::U3,vcs::Matrix>> Kmatrix_validation_table = {
{
	u3::U3(14,{3,1}),
	{
		{u3::U3(14,{3,1}),vcs::Matrix{{ 1.00000000}}},
		{u3::U3(16,{3,2}),vcs::Matrix{{ 3.00000000}}},
		{u3::U3(16,{4,0}),vcs::Matrix{{ 2.64575131}}},
		{u3::U3(18,{2,2}),vcs::Matrix{{ 6.80213819,-0.25347458},{-0.25347458,  4.96008239}}}, 
		{u3::U3(20,{1,2}),vcs::Matrix{{12.72792206,-0.00000000},{-0.00000000, 12.72792206}}}, 
		{u3::U3(20,{3,4}),vcs::Matrix{{33.82805696,-0.68011935},{-0.68011935, 25.58783769}}},
		{u3::U3(20,{2,3}),vcs::Matrix{{19.62822980, 2.16141450,-0.24674332},{ 2.16141450, 19.15191779, -0.72961120},{-0.24674332,-0.72961120, 13.68965979}}}, 
		{u3::U3(20,{3,1}),vcs::Matrix{{20.33536264,-1.99295318, 0.03411572},{-1.99295318, 16.44081762, -0.63570756},{ 0.03411572,-0.63570756, 15.37589807}}} 
	}
},
{
	u3::U3(HalfInt(29,2),{2,1}),
	{
		{u3::U3(HalfInt(29,2),{2,1}),vcs::Matrix{{ 1.00000000}}},
		{u3::U3(HalfInt(33,2),{0,3}),vcs::Matrix{{ 2.64575131}}},
		{u3::U3(HalfInt(37,2),{1,2}),vcs::Matrix{{ 7.39228476,-0.14419682},{-0.14419682,5.68441207}}},
		{u3::U3(HalfInt(37,2),{3,1}),vcs::Matrix{{ 8.20961975,-0.51846909},{-0.51846909,7.00460728}}},
		{u3::U3(HalfInt(41,2),{2,1}),vcs::Matrix{{23.07368495,-1.24931673, 0.21040579},{-1.24931673,18.44600494,-0.42907995},{0.21040579,-0.42907995,17.93799375}}}
	}
}
};

/// Compare two functions for calculating Kmatrices
if(true)
{
	std::vector<u3::U3> sigma_vector={u3::U3(7,4,3), u3::U3(HalfInt(13,2),HalfInt(9,2),HalfInt(7,2))};
	for(const auto& [sigma,test_values] : Kmatrix_validation_table)
	{
		sp3r::Sp3RSpace sp3r_irrep(sigma,40);
		
		std::map<u3::U3, MultiplicityTagged<u3::U3>::vector> u3_subspace_map;
		for(const auto& u3_subspace : sp3r_irrep)
			{
				std::vector<u3::U3> n_vector;
				const auto& omega = u3_subspace.labels();
				for(auto i : iter::range(0,int(u3_subspace.size()),1))
					{
						const auto&[n,dummy] = u3_subspace.GetStateLabels(i);
						n_vector.push_back(n);
					}

				auto unique_n = iter::unique_everseen(n_vector);
				for(const auto& n : unique_n)
					u3_subspace_map[omega].push_back({n,u3::OuterMultiplicity(sigma,n,omega)});
			}

		vcs::MatrixCache K1_matrix_map;
		vcs::GenerateKMatrices(sp3r_irrep,K1_matrix_map);
		vcs::KmatrixMap K2_matrix_map = vcs::GenerateKMatrices(sigma,u3_subspace_map);

		for(const auto& [omega,KK] : K1_matrix_map)
			if(omega.N()-sigma.N()>36)
			{
				// std::cout<<omega.Str()<<std::endl;
				// std::cout<<KK<<std::endl<<std::endl;
				// vcs::Matrix Kinv = KK.inverse();
				// std::cout<<Kinv<<std::endl<<std::endl;
				// std::cout<<KK*Kinv<<std::endl;
				double factor = pow(2,double(omega.N()-sigma.N()));
				vcs::Matrix KKf = KK/factor;
				std::cout<<"factor "<<factor<<std::endl;
				std::cout<<KK.inverse()*factor<<std::endl<<KKf.inverse()<<std::endl<<std::endl;//<<KKf*KKf.inverse()<<std::endl<<std::endl;


			}


		for(const auto& [omega,K] : test_values)
		// for(const auto& [omega,K1] : K_matrix_map1)
			{
				const auto& K1=K1_matrix_map.at(omega);
				const auto& K2=K2_matrix_map.at(omega)[0];
				assert(mcutils::IsZero(K1-K,1e-6));
				assert(mcutils::IsZero(K2-K,1e-6));
			}
	}
}






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
	sp3r::Sp3RSpace irrep(s,10);

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

if(false)
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


