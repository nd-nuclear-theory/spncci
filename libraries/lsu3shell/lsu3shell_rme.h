/****************************************************************
  lsu3shell_rme.h

  Input of RMEs in lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  - 8/1/16 (aem,mac): Created.
  - 9/7/16 (mac): Split from lsu3shell_interface.
  - 9/8/16 (mac): Add operator RME comparison function.
  - 2/17/16 (mac): Change functions to take filename and do open
    failure checking.
  - 6/4/17 (mac): Provide optional scale factor for rmes on input.
  - 6/11/17 (mac):
    + Update for binary rme file format.
    + Eliminate deprecated forms of ReadLSU3ShellRMEs and
      GenerateLSU3ShellNcmRMEs taking stream arguments.
  - 6/12/19 (mac): Make text/binary switchable by global mode flag
    g_rme_binary_format.
  - 6/17/19 (mac): Update binary mode output: add header, shorten indexing
    integer type, make storage single/double switchable.

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

  // global configuration
  extern bool g_rme_binary_format;

  void 
  ReadLSU3ShellRMEs(
      const std::string& filename,
      const LSU3BasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space, 
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      basis::MatrixVector& matrix_vector,
      double scale_factor = 1.
    );
  // Read LSU3Shell RMEs.
  //
  // Global mode:
  //   g_rme_binary_format: false = ascii format; true = binary format
  //
  // Arguments:
  //   filename (input): name of binary rme file
  //   lsu3_basis_table (input): information on LSU3shell basis states
  //   space (input): U3SPN basis
  //   operator_labels (input): U3ST tensorial labels for operator being input
  //   sectors (input): U3SPN sectors for operator being input
  //   matrix_vector (output): matrices for operator blocks
  //   scale_factor (input,optional): scale factor for rmes on input,
  //     e.g., for conversion from relative to intrinsic operator


  // For reference, here is the old syntax:
  //
  // void 
  // ReadLSU3ShellRMEs(
  //     std::ifstream& is,
  //     const u3shell::OperatorLabelsU3ST& operator_labels,
  //     const LSU3BasisTable& lsu3_basis_table,
  //     const u3shell::SpaceU3SPN& space, 
  //     const u3shell::SectorsU3SPN& sectors,
  //     basis::MatrixVector& matrix_vector,
  //     double scale_factor = 1.
  //   );
  //
  // This required that the calling code take care of the stream open,
  // which also had the drawback that there is no standardized error
  // message for file open failures.
  //
  // When upgrading calling code to use new version, beware that order
  // of the other arguments has also been adjusted.


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

  // For reference, here is the old syntax:
  //
  // void GenerateNcmMatrixVector(
  //   int A,      
  //   std::ifstream& is_Nrel,
  //   const lsu3shell::LSU3BasisTable& lsu3_basis_table,
  //   const u3shell::SpaceU3SPN& space, 
  //   basis::MatrixVector& matrix_vector 
  // );

}
#endif
