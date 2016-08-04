/****************************************************************
  sp_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "lsu3shell_interface.h"


#include <fstream>
#include <iostream>

#include "cppformat/format.h"

#include "u3shell/two_body_operator.h"
#include "u3shell/moshinsky.h"
#include "spncci/sp_basis.h"


namespace lsu3shell
{
	void ReadLSU3Vector(const std::string& filename, LSU3Vector& lsu3basis_vector)
	{
		std::ifstream basis_stream(filename.c_str());
		std::string line;
		while(std::getline(basis_stream,line))
			{
				std::istringstream line_stream(line);

				int Nex, twice_Sp, twice_Sn, twice_S, lambda, mu, index;
				line_stream >> index >> Nex >>twice_Sp >> twice_Sn >>twice_S >> lambda >> mu;

				//conversions
				HalfInt Sp=HalfInt(twice_Sp,2);
				HalfInt Sn=HalfInt(twice_Sn,2);
				HalfInt S=HalfInt(twice_S,2);
				u3::SU3 x(lambda,mu);
				lsu3basis_vector.push_back(lsu3shell::LSU3Irrep(Nex,x,Sp,Sn,S));
			}	
		basis_stream.close();
	}

	void GenerateLSU3ShellOperators(int Nmax, const u3shell::RelativeUnitTensorCoefficientsU3ST& relative_tensor_expansion, int operator_index)
	{
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
	  // declare coefficient containers
    u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
    u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
    //moshinsky transform and accumulate coefficients
    u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients
	  	=u3shell::TransformRelativeTensorToTwobodyTensor(relative_tensor_expansion,twobody_space);
    // convert to biquads
    u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
    // convert to pn scheme
    u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn);
    
    std::string operator_stream_filename = fmt::format("operator{:06d}.recoupler",operator_index);
    std::ofstream operator_stream(operator_stream_filename);
    WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
    operator_stream.close();
	}

	void GenerateLSU3ShellOperators(int Nmax, const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_tensor_labels)
		{
    u3shell::TwoBodySpaceU3ST  twobody_space(Nmax);
	  for(int i=0; i<relative_tensor_labels.size(); ++i)
		  {
		  	// declare coefficient containers
	      u3shell::TwoBodyUnitTensorCoefficientsU3ST two_body_unit_tensor_coefficients;
	      u3shell::TwoBodyUnitTensorCoefficientsU3ST biquad_coefficients;
	      u3shell::TwoBodyUnitTensorCoefficientsU3SPN biquad_coefficients_pn;
	      // Moshinsky transform unit tensors to two-body operators 
	  		MoshinskyTransformUnitTensor(relative_tensor_labels[i], 1.0, twobody_space,two_body_unit_tensor_coefficients);
	      // convert to biquads
	      u3shell::TransformTwoBodyUnitTensorToBiquad(two_body_unit_tensor_coefficients,biquad_coefficients);
	      // convert biquads to pn scheme
	      u3shell::TransformBiquadToPNScheme(biquad_coefficients,biquad_coefficients_pn);
	      // std::cout<<fmt::format("{} operator{:06d}",relative_unit_tensor_labels.Str(),i)<<std::endl;
	      std::string operator_stream_filename = fmt::format("operator{:06d}.recoupler",i);
	      std::ofstream operator_stream(operator_stream_filename);
	      WriteTwoBodyOperatorRecoupler(operator_stream,biquad_coefficients_pn);
	      operator_stream.close();
		  }
		}

	void ReadLSU3ShellRMEs(
				std::ifstream& is,
				int dimp, int dim, const u3::SU3& x0,
				const lsu3shell::LSU3Vector& lsu3shell_vector_bra,
				const lsu3shell::LSU3Vector& lsu3shell_vector_ket,
				std::vector<Eigen::MatrixXd>& matrix_vector
			)
	{
		u3::SU3 x, xp;
		for(int i=0; i<dimp; ++i)
			for(int j=0; j<dim; ++j)
				{
		  		xp=lsu3shell_vector_bra[i].x();
		  		x=lsu3shell_vector_ket[j].x();
		  		int rho0_max=u3::OuterMultiplicity(x,x0,xp);
		  		for(int rho0=1; rho0<=rho0_max; ++rho0)
		  			{
		  				if (rho0<matrix_vector.size())
		  					matrix_vector.push_back(Eigen::MatrixXd::Zero(dimp,dim));	
		  				is>>matrix_vector[rho0](i,j);
	  				}
	  		}
	}

	void 
	WriteLSU3ShellExpansionLGI(
			const std::map<lsu3shell::LSU3Irrep,int>& lgi_lsu3shell_map, 
			const Eigen::MatrixXd& kernel,
			std::ofstream& os
		)
	{
		// int row_dim=kernel.rows();
		// for(auto it=lgi_lsu3shell_map.begin(); it!= lgi_lsu3shell_map.end(); ++it)
		// 	{
		// 		int Nex;
		// 		u3::SU3 x;
		// 		HalfInt Sp, Sn, S;
		// 		std::tie(Nex,x,Sp,Sn,S)=it->first;
		// 		int column=it->second;
		// 		int num_nonzero=0;
		// 		for(int row=0; row<row_dim; ++row)
		// 				if(fabs(kernel(row,column))>10e-8)
		// 					++num_nonzero;
		// 		for(int row=0; row<row_dim; ++row)
		// 			{

		// 			}
		// 	}
		//Take map of LGI labels and column in L SU3shell basis
		// iterate over map
		//write out non-zero matrix elements
		//lsu3shell::LSU3Irrep(Nex,x,Sp,Sn,S)
	}

	void GenerateLSU3ShellExpansionLGI(
		int Nsigma_min,
		int Nsigma_max, 
		std::string basis_file, 
		std::string brel_filename, 
		std::string nrel_filename,
		std::string lgi_filename
	)
	// Construct Brel and Ncm matrix in lsu3shell basis and 
	// solve for null space 
	// write expansions to file 
	// std::string brel_file="operator.000000";
	// std::string brel_file="operator.000001";
	{
		LSU3Vector lsu3basis_vector;
		lsu3shell::ReadLSU3Vector(basis_file, lsu3basis_vector);
		int Nex=0;
		u3::SU3 x; 
		HalfInt Sp, Sn, S;
		int index=0;
		while(Nex==0)
			{
				std::tie(Nex, x, Sp, Sn, S)=lsu3basis_vector[index].Key();
			}
		spncci::LGIVectorType lgi_vector;
		//if starting at Nsigma=0, use LSU3 irreps to fill out LGIVector for Nsigma=0;
		if(Nsigma_min==0)
			{
				Nsigma_min=1;
				// fill in vector
			}
		// Read in already generated LGI's from file	
		else
				spncci::GenerateLGIVector(lgi_vector,lgi_filename, Nsigma_min);

		// For each Nsigma
		for(int Nsigma=Nsigma_min; Nsigma<=Nsigma_max; ++Nsigma)
		{
			// std::ofstream& os 

			////////////////////////////////////////////////////////////////////////
			// Setting up the matrix 
			////////////////////////////////////////////////////////////////////////
			int dim=lsu3basis_vector.size();

			int dim2,dim0;
			std::ifstream is_brel(brel_filename.c_str());
			is_brel>>dim0>>dim2;
			// if((dim0+dim2)!=dim)
			// 	std::cout<<"brel dimension mismatch"<<std::endl;
			Eigen::MatrixXd BNcmMatrix=Eigen::MatrixXd::Zero(dim0+dim2,dim2);
			////////////////////////////////////////////////////////////////////////
			//Read in matrix elements of Brel
			////////////////////////////////////////////////////////////////////////
			double rme;
			for(int i=0; i<dim0; ++i)
				for(int j=0; j<dim2; ++j)
					{
						is_brel>>rme;
						BNcmMatrix(i,j)=rme;
					}
			is_brel.close();
			////////////////////////////////////////////////////////////////////////
			//Read in matrix elements of Brel
			////////////////////////////////////////////////////////////////////////
			int dim01,dim02;
			std::ifstream is_nrel(nrel_filename.c_str());
			is_nrel>>dim01>>dim02;

			//Sanity check
			if((dim01!=dim02)||(dim01!=dim0))
				std::cout<<"dimension miss match between Brel and Nrel"<<std::endl;

			for(int i=0; i<dim2; ++i)
				for(int j=0; j<dim2; ++j)
					{
						is_nrel>>rme;
						//Ncm=Ntotal-Nre;;
						BNcmMatrix(dim0+i,j)=2-rme;
					}
			////////////////////////////////////////////////////////////////////////
			//Solve for null space 
			////////////////////////////////////////////////////////////////////////
			Eigen::MatrixXd BNcm_kernal=BNcmMatrix.fullPivLu().kernel();

			//Get list of LGI labels and append to lgi structure 

			////////////////////////////////////////////////////////////////////////
			//write to file  
			////////////////////////////////////////////////////////////////////////
		}
	}

	void WriteControlFile()
	{}
}// end namespace