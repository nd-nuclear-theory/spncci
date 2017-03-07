/****************************************************************
  io_control.h

  High-level control code for data I/O in SpNCCI calculation programs.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).

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
        const std::string& Brel_filename,
        u3shell::SectorsU3SPN& Brel_sectors,
        basis::MatrixVector& Brel_matrices,
        const std::string& Arel_filename,
        u3shell::SectorsU3SPN& Arel_sectors,
        basis::MatrixVector& Arel_matrices,
        const std::string& Nrel_filename,
        u3shell::SectorsU3SPN& Nrel_sectors,
        basis::MatrixVector& Nrel_matrices
      );
  // Read lsu3shell RMEs for basic symplectic generator
  // operators.
  //
  // Reads RMEs for Brel, Arel, and Nrel, plus deduces Ncm.
  //
  // Arguments:
  //   lsu3shell_basis_table (input): lsu3shell basis data
  //   lsu3shell_space (input): lsu3shell basis
  //   <operator>_filename (input): rme filename
  //     for <operator> in Arel, Brel, Nrel
  //   <operator>_sectors (output): U3SPN sectors
  //   <operator>_matrices (output): matrices of RMEs


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
