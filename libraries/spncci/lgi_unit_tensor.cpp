/****************************************************************
  lgi_unit_tensor.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "spncci/lgi_unit_tensor.h"


#include <iostream>
#include <fstream>

#include "cppformat/format.h"

#include "sp3rlib/u3coef.h"
#include "sp3rlib/vcs.h"
#include "u3shell/moshinsky.h"
#include "u3shell/u3st_scheme.h"
#include "u3shell/two_body_operator.h"
#include "lsu3shell/lsu3shell_interface.h"

// #include "spncci/unit_tensor.h"

namespace spncci
{

	void ReadLGISU3Expansion(std::string filename, Eigen::MatrixXd& matrix, std::string type)
	{
    int lgi_index, num_nonzero_coefs, lsu3_index, dim;
    double coef;
		std::ifstream is(filename.c_str());
  	if(!is)
    	std::cout<<"Didn't open"<<std::endl;
    is>>dim;
    matrix(dim,dim);
		while(is)
			{
				is>>lgi_index>>num_nonzero_coefs;
				for(int num=0; num<num_nonzero_coefs; num++)
					{
						is>>lsu3_index>>coef;
						if(type=="ket")
							matrix(lsu3_index,lgi_index)=coef;
						else if (type=="bra")
							matrix(lgi_index,lsu3_index)=coef;
						else
							std::cout<<"invalid type, please indicate if bra or ket" <<std::endl;
					}
			}
	}

	void PopulateUnitTensorLGISectors(
					std::string filename,
					const lsu3shell::LSU3Vector& lsu3shell_vector_bra,
					const lsu3shell::LSU3Vector& lsu3shell_vector_ket,
					const spncci::LGIVectorType& lgi_vector,
					const Eigen::MatrixXd& bra_matrix,
					const Eigen::MatrixXd& ket_matrix,
	    		const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels,
	    		std::vector< std::pair<int,int> >& lgi_pair_vector,
				  std:: map< 
				    std::pair<int,int>,
				    std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache >
				    >& lgi_unit_tensor_rme_map
				)
	// lsu3shell_vector is a list of all the basis state labels, contained in same struture 
	// as is used for lgi's 
	{
		int i,j;
		double rme;
		u3::SU3 x0;
		std::pair<int,int>NpNSector(0,0);
		std::ifstream is(filename.c_str());
		if(!is)
	  	std::cout<<"Didn't open"<<std::endl;
	  Eigen::MatrixXd temp_matrix(1,1);
	  
		u3shell::RelativeUnitTensorLabelsU3ST tensor;
		int lsu3_dimp=lsu3shell_vector_bra.size(),lsu3_dim=lsu3shell_vector_ket.size();
		int lgi_dimp=bra_matrix.cols(), lgi_dim=ket_matrix.rows();
		for(int t=0; t<relative_unit_tensor_labels.size(); ++t)
			{
				//Generate matrix vector for each realtive_unit_tensor
				tensor=relative_unit_tensor_labels[t];
				x0=tensor.x0();
				std::vector<Eigen::MatrixXd> matrix_vector;
				
				lsu3shell::ReadLSU3ShellRMEs(
						is, lsu3_dimp, lsu3_dim, x0,lsu3shell_vector_bra,
						lsu3shell_vector_ket,matrix_vector
					);
				//TODO:: Make into function
				// Populating cache 
				Eigen::MatrixXd lgi_rme_matrix;
				for(int rho0=1; rho0<=matrix_vector.size(); ++rho0)
					{
						lgi_rme_matrix=bra_matrix*matrix_vector[rho0]*ket_matrix;
						for(int i=0; i<lgi_dimp; ++i)
							for(int j=0; j<lgi_dim; ++j)
								{
									u3::U3 sigmap(lgi_vector[i].sigma), sigma(lgi_vector[j].sigma);
									std::pair<int,int>lgi_pair(i,j);
									// spncci::UnitTensor relative_tensor=spncci::RelativeUnitTensor(tensor);
									spncci::UnitTensorU3Sector key(sigmap,sigma,spncci::RelativeUnitTensor(tensor),rho0);
									temp_matrix(0,0)=lgi_rme_matrix(i,j);
									lgi_unit_tensor_rme_map[lgi_pair][NpNSector][key]=temp_matrix;
									if (fabs(lgi_rme_matrix(i,j))>10e-7)
										lgi_pair_vector.push_back(lgi_pair);
								}
					}

			}
	}

} // end namespace