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
		sp3r::Sp3RSpace sp3r_irrep(sigma,20);
		
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

		for(const auto& [omega,K] : test_values)
		// for(const auto& [omega,K1] : K_matrix_map1)
			{
				const auto& K1=K1_matrix_map.at(omega);
				const auto& K2=K2_matrix_map.at(omega)[0];
				// std::cout<<K2<<std::endl<<K<<std::endl;
				assert(mcutils::IsZero(K1-K,1e-6));
				assert(mcutils::IsZero(K2-K,1e-6));
			}
	}
}

int a = 1 << 10;
std::cout<<a<<" or "<<pow(2,10)<<std::endl;


} // main 


