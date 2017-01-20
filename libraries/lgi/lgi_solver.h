/****************************************************************
  lgi_solver.h

  Interface for lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.
  ...
  1/19/17 (mac): Replace null solver.
****************************************************************/

#ifndef LGI_PROB_H_
#define LGI_PROB_H_
#include <vector>
#include <eigen3/Eigen/Eigen>
#include "am/am.h"  
#include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_rme.h"

namespace lgi
{
  void 
    GenerateNcmMatrixVector(
        int A,     
        std::ifstream& is_nrel,
        const lsu3shell::LSU3BasisTable& lsu3_basis_table,
        const u3shell::SpaceU3SPN& space, 
        basis::MatrixVector& matrix_vector 
      );
  // Generates vector of Ncm Matrix sectors in LSU3shell basis from 
  // Nrel matrix sectors. 

  void 
    GenerateBrelNcmMatrices(
        int A,
        std::ifstream& is_brel,
        std::ifstream& is_nrel,
        const lsu3shell::LSU3BasisTable& lsu3_basis_table,
        const u3shell::SpaceU3SPN& space, 
        basis::MatrixVector& BrelNcm_vector 
      );
  // Constructs the Brel+Ncm Matrices for each ket subspace and
  // stores them in BrelNcm_vector

  void 
    GenerateLGIExpansion(
        int A,
        HalfInt Nsigma_0,
        const lsu3shell::LSU3BasisTable& lsu3_basis_table,
        const u3shell::SpaceU3SPN& space, 
        std::ifstream& is_brel,
        std::ifstream& is_nrel,
        lgi::LGIVector& lgi_vector,
        basis::MatrixVector& lgi_expansion_matrix_vector,
        bool keep_empty_subspaces=false  
      );
  // Generates the LGI Expansion in terms of lsu3shell SU(3)xSU(2) 
  // reduced basis states by solving for the null space of the Brel+Ncm matrix
  // 
  //  lsu3_basis_table and space are generated as output from lsu3shell::ReadLSU3Basis
  //  which reads in basis table generated by ncsmSU3xSU2BasisLSU3Tabular
  //
  //  Arguments: 
  //    A (input) : atomic mass number of nucleus
  //    lsu3_basis_table (input) : lookup table of LSU3Shell basis states between
  //                                LSU3Shell basis and U3SPN space. 
  //    space (input) : space on which brel and nrel rme's are calculated
  //    brel_filename (input) : file containing lsu3shell rme's of Brel
  //    nrel_filename (input) : file containing lsu3shell rme's of Nrel
  //    

  void
    TransformOperatorToSpBasis(
        const u3shell::SectorsU3SPN& sectors,
        const basis::MatrixVector& basis_transformation_matrices,
        const basis::MatrixVector& lsu3shell_operator_matrices,
        basis::MatrixVector& spncci_operator_matrices
      );
  // Similarity transformation from LSU3shell basis to Sp(3,R)xSU(2) basis

  void
    WriteLGILabels(const lgi::LGIVector& lgi_vector, std::ofstream& os);


}
#endif
