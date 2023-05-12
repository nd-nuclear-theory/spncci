/****************************************************************
  lsu3shell_rme.h

  Input of RMEs in lsu3shell basis.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

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
  - 10/11/17 (aem):
      + Add MatrixFloatType and switched MatrixVector to OperatorBlocks
      + Extracted i/o for symplectic generator su3rmes from spncci_io
****************************************************************/

#ifndef LSU3SHELL_RME_H_
#define LSU3SHELL_RME_H_

#include <Eigen/Dense>

#include "am/am.h"
#include "basis/operator.h"
#include "u3shell/u3spn_scheme.h"

#include "lsu3shell/lsu3shell_basis.h"

namespace lsu3shell
{

  // global configuration
  extern bool g_rme_binary_format;
  typedef unsigned int RMEIndexType;

  // matrix precision
  typedef double MatrixFloatType;
  typedef basis::OperatorBlock<MatrixFloatType> OperatorBlock;
  typedef basis::OperatorBlocks<MatrixFloatType> OperatorBlocks;


  void
  ReadLSU3ShellRMEs(
      bool sp3r_generators,
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table_bra,
      const u3shell::SpaceU3SPN& space_bra,
      const LSU3ShellBasisTable& lsu3_basis_table_ket,
      const u3shell::SpaceU3SPN& space_ket,
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks,
      double scale_factor
    );

  // Read LSU3Shell RMEs.
  //
  // Global mode:
  //   g_rme_binary_format: false = ascii format; true = binary format
  //
  // Arguments:
  //   sp3r_generators (input) : tempory fix for error arising when generator rmes are not
  //      truncated as a string.  If true, then string truncation applied, else, read in normally
  //   filename (input): name of binary rme file
  //   lsu3_basis_table_bra,lsu3_basis_table_ket (input): information on LSU3shell basis states
  //   space_bra,space_bra (input): U3SPN bases
  //   operator_labels (input): U3ST tensorial labels for operator being input
  //   sectors (input): U3SPN sectors for operator being input
  //   matrix_vector (output): matrices for operator blocks
  //   scale_factor (input,optional): scale factor for rmes on input,
  //     e.g., for conversion from relative to intrinsic operator


  void
  ReadLSU3ShellRMEs(
      bool sp3r_generators,
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space,
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks,
      double scale_factor=1.
    );

  //Wrapper function for when bra space and ket space are the same.


  void
  ReadLSU3ShellRMEs(
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space,
      const u3shell::OperatorLabelsU3ST& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& matrix_vector,
      double scale_factor = 1.
    );
  // Wraper function for when sp3r_generators not specified.

//******************************************** Added by J.H. **************************************
  void ReadLSU3ShellRMEs(
      const std::string& filename,
      const LSU3ShellBasisTable& lsu3_basis_table,
      const u3shell::SpaceU3SPN& space,
      const u3shell::OperatorLabelsU3S& operator_labels,
      const u3shell::SectorsU3SPN& sectors,
      lsu3shell::OperatorBlocks& blocks
    );
//*************************************************************************************************

  void
    ReadLSU3ShellSymplecticOperatorRMEs(
        const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space,
        const std::string& Brel_filename, u3shell::SectorsU3SPN& Bintr_sectors, lsu3shell::OperatorBlocks& Bintr_matrices,
        const std::string& Nrel_filename, u3shell::SectorsU3SPN& Nintr_sectors, lsu3shell::OperatorBlocks& Nintr_matrices,
        int A
      );
    // Reads in Brel and Nrel su3rme files and applies appropriate scaling to convert to intrinsic operators using
    // mechanics convention.

  void
    ReadLSU3ShellSymplecticRaisingOperatorRMEs(
        const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space,
        const std::string& Arel_filename, u3shell::SectorsU3SPN& Aintr_sectors, lsu3shell::OperatorBlocks& Aintr_matrices,
        int A
      );
    // Reads in Arel su3rme files and applies appropriate scaling to convert to Aintr under mechanics convention.

  bool
  CompareLSU3ShellRMEs(
      std::ostream& log_stream,
      const U3SPNBasisLSU3Labels& basis_provenance,
      const u3shell::SpaceU3SPN& space,
      const u3shell::SectorsU3SPN& sectors,
      const lsu3shell::OperatorBlocks& matrices1,
      const lsu3shell::OperatorBlocks& matrices2,
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
    const lsu3shell::OperatorBlocks& Nrel_matrices,
    int A,
    lsu3shell::OperatorBlocks& Ncm_matrices
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
  //   const lsu3shell::LSU3ShellBasisTable& lsu3_basis_table,
  //   const u3shell::SpaceU3SPN& space,
  //   lsu3shell::OperatorBlocks& matrix_vector
  // );

}
#endif
