/****************************************************************
  lgi_solver.h

  Interface for lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.

****************************************************************/

#ifndef LGI_SOLVER_H_
#define LGI_SOLVER_H_

#include <boost/functional/hash_fwd.hpp>
#include <eigen3/Eigen/Eigen>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/u3spn_scheme.h"
#include "lsu3shell/lsu3shell_rme.h"


namespace lgi
{
  void 
  GenerateNcmMatrixVector(      
    const std::string& nrel_filename,
    const lsu3shell::LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    basis::MatrixVector& matrix_vector 
    );
  // Generates vector of Ncm Matrix sectors in LSU3shell basis from 
  // Nrel matrix sectors. 

  void 
  GenerateBrelNcmMatrices(
      const std::string& brel_filename,
      const std::string& nrel_filename,
      const lsu3shell::LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      basis::MatrixVector& BrelNcm_vector 
    );
  // Constructs the Brel+Ncm Matrices for each ket subspace and
  // stores them in BrelNcm_vector

  void 
  GenerateLGIExpansion(
    const lsu3shell::LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    const std::string& brel_filename,
    const std::string& nrel_filename,
    basis::MatrixVector& lgi_expansion_matrix_vector  
  );
  // Generates the LGI Expansion in terms of lsu3shell SU(3)xSU(2) reduced basis states
  // by solving for the null space of the Brel+Ncm matrix

  void
  TransformOperatorToSpBasis(
    u3shell::OperatorLabelsU3S& operator_labels,
    std::string operator_file,
    basis::MatrixVector& spncci_operator_vector
  );
  // Similarity transformation from LSU3shell basis to Sp(3,R)xSU(2) basis

}
#endif
