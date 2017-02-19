/****************************************************************
  lsu3shell_rme.h

  Input of RMEs in lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.
  9/8/16 (mac): Add operator RME comparison function.
  2/17/16 (mac): Change functions to take filename and do open
    failure checking.
****************************************************************/

#ifndef LSU3SHELL_RME_H_
#define LSU3SHELL_RME_H_

#include "boost/functional/hash_fwd.hpp"
#include "eigen3/Eigen/Dense"

#include "am/am.h"  
#include "basis/operator.h"
#include "u3shell/u3spn_scheme.h"

#include "lsu3shell/lsu3shell_basis.h"

namespace lsu3shell
{

  void 
  ReadLSU3ShellRMEs(
      const std::string& filename,
      const LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      basis::MatrixVector& matrix_vector
    );

  void 
  ReadLSU3ShellRMEs(
      std::ifstream& is,
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& sectors,
      basis::MatrixVector& matrix_vector // in operator.h and initial to zero
    );
  // DEPRECATED: requires external code to do stream open; when
  // upgrading calling code to use new version, beware that order of
  // arguments has also been adjusted (to keep all the variables
  // describing the operator together in the argument list)

  bool 
  CompareLSU3ShellRMEs(
      std::ostream& log_stream,
      const U3SPNBasisLSU3Labels& basis_provenance,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::SectorsU3SPN& sectors,
      const basis::MatrixVector& matrices1,
      const basis::MatrixVector& matrices2,
      double tolerance,
      bool verbose = false
    );
  // Make diagnostic comparison report
  //
  // Arguments:
  //   log_stream : output stream for diagnostics
  //   space : space on which operators are defined
  //   sectors : sectors on which operators are defined
  //   matrices1 : matrices for first operator
  //   matrices2 : matrices for second operator
  //   tolerance : max acceptable magnitude for entrywise residual
  //   verbose (optional): whether or not to print sector header lines
  //
  // Returns:
  //   (bool): comparison OK


  void GenerateLSU3ShellNcmRMEs(
    const u3shell::SpaceU3SPN& space,
    const u3shell::SectorsU3SPN& Nrel_sectors,
    const basis::MatrixVector& Nrel_matrices,
    int A,
    basis::MatrixVector& Ncm_matrices
  );
  // Generates the Ncm matrix elements from Nrel matrix elements
  // and stores them in a vector of lsu3shell basis sectors
  //
  // Arguments:
  //  space (input): space defined by lsu3shell basis
  //  Nrel_sectors (input): sectors for Nrel operator
  //  Nrel_matrices (input): matrices for Nrel operator
  //  A (input): atomic mass number
  //  matrix_vector (output): container for Ncm matrix sectors.

  void GenerateNcmMatrixVector(
    int A,      
    std::ifstream& is_Nrel,
    const lsu3shell::LSU3BasisTable& lsu3_basis_table,
    const u3shell::SpaceU3SPN& space, 
    basis::MatrixVector& matrix_vector 
  );
  // DEPRECATED: requires external code to do stream open
  //
  // Generates the Ncm matrix elements from Nrel matrix elements
  // and stores them in a vector of lsu3shell basis sectors
  //
  // Arguments:
  //  A (input): atomic mass number
  //  is_Nrel (input): stream from file containing Nrel rmes
  //  lsu3_basis_table (input): table of lsu3shell basis states
  //  space (input): space defined by lsu3shell basis
  //  matrix_vector (output): container for Ncm matrix sectors.

}
#endif
