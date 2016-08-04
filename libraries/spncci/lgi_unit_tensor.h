/****************************************************************
  lgi_unit_tensor.h

  Generate unit tensors for LGI's 
  --Generate list of Unit tensor between LGI's
  --Generate operator files for LSU3shell
  --Reads in LSU3shell output and populates unit tensor cache
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  7/28/16 (aem,mac): Created.
****************************************************************/
#ifndef LGI_UNIT_TENSOR_H_
#define LGI_UNIT_TENSOR_H_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash_fwd.hpp>
#include <eigen3/Eigen/Eigen>

#include "spncci/sp_basis.h"
#include "sp3rlib/u3.h"
#include "sp3rlib/u3coef.h"
#include "spncci/unit_tensor.h"
#include "u3shell/relative_operator.h"


namespace spncci
{
void GenerateLSU3ShellOperators(
		int Nmax, 
		std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels
		);
	// Generate input files for LSUshell recoupler for all relative 
	// unit tensor's which may have non-zero matrix elements between
	// LGI's. 

void ReadLGISU3Expansion(std::string filename, Eigen::MatrixXd& matrix, std::string type);
		// Read in LGI Expansion in terms of LSU3shell basis states 
		// Stores in matix with dimenisons depending on type:
		// 		if type="ket", then matrix dimensions are (lsu3shell_dim, num_lgi)
		//		if type="bra", then matrix diemsnions are (num_lgi, lsu3shell_dim)

void PopulateUnitTensorLGISectors(
			std::string filename,
			const spncci::LGIVectorType& lsu3shell_vector_bra,
			const spncci::LGIVectorType& lsu3shell_vector_ket,
			const spncci::LGIVectorType& lgi_vector,
			const Eigen::MatrixXd& bra_matrix,
			const Eigen::MatrixXd& ket_matrix,
  		const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& relative_unit_tensor_labels,
		  std::vector< std::pair<int,int> >& lgi_pair_vector,
		  std:: map< 
		    std::pair<int,int>,
		    std::map<std::pair<int,int>,spncci::UnitTensorSectorsCache >
		    >& lgi_unit_tensor_rme_map
		);
		// Reads in lsu3shell rme's and calculates lgi rme's
		//Caches lgi rme's into the lgi_unit_tensor_rme_map

} // end namespace

#endif