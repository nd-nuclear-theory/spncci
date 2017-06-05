/****************************************************************
  io_control.h

  High-level control code for data I/O in SpNCCI calculation programs.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).
  6/4/17 (mac): Revise to rescale rmes from relative to intrinsic operators
    on input
****************************************************************/

#ifndef SPNCCI_SPNCCI_IO_CONTROL_H_
#define SPNCCI_SPNCCI_IO_CONTROL_H_

#include <string>

#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "u3shell/relative_operator.h"
#include "u3shell/u3spn_scheme.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // reading lsu3shell RMEs
  ////////////////////////////////////////////////////////////////

  void
    ReadLSU3ShellSymplecticOperatorRMEs(
        const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space, 
        const std::string& Brel_filename, u3shell::SectorsU3SPN& Bintr_sectors, basis::MatrixVector& Bintr_matrices,
        const std::string& Arel_filename, u3shell::SectorsU3SPN& Aintr_sectors, basis::MatrixVector& Aintr_matrices,
        const std::string& Nrel_filename, u3shell::SectorsU3SPN& Nintr_sectors, basis::MatrixVector& Nintr_matrices,
        int A
      );
  // Read lsu3shell RMEs for basic symplectic generator
  // operators.
  //
  // Reads RMEs for Brel, Arel, and Nrel, plus deduces Ncm.
  //
  // Arguments:
  //   lsu3shell_basis_table (input): lsu3shell basis data
  //   lsu3shell_space (input): lsu3shell basis
  //   For <operator> in {A,B,N}:
  //     <operator>rel_filename (input): rme filename
  //     <operator>intr_sectors (output): U3SPN sectors
  //     <operator>intr_matrices (output): matrices of RMEs
  //   A (input): nucleon number

  void
  ReadLSU3ShellSeedUnitTensorRMEs(
      const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels,
      const std::string& filename,
      u3shell::SectorsU3SPN& unit_tensor_sectors,
      basis::MatrixVector& unit_tensor_lsu3shell_matrices
    );

  // void
  //   ReadLSU3ShellSeedUnitTensorRMEs(
  //       const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
  //       const u3shell::SpaceU3SPN& lsu3shell_space, 
  //       const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
  //       const std::string& relative_unit_tensor_filename_template,
  //       std::vector<u3shell::SectorsU3SPN>& lgi_unit_tensor_sectors,
  //       std::vector<basis::MatrixVector>& lgi_unit_tensor_lsu3shell_matrices
  //     );
  // Read lsu3shell RMEs for seed unit tensors.
  //
  // Arguments:
  //   lsu3shell_basis_table (input): lsu3shell basis data
  //   lsu3shell_space (input): lsu3shell basis
  //   lgi_unit_tensor_labels (input): labels for the unit tensors to read
  //   relative_unit_tensor_filename_template (input): filename template for use with fmt::format
  //   lgi_unit_tensor_sectors (output): U3SPN sectors (for each unit tensor)
  //   lgi_unit_tensor_lsu3shell_matrices (output): matrices of RMEs (for each unit tensor)

}  // namespace

#endif
