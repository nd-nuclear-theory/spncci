/****************************************************************
  io_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/io_control.h"

#include "cppformat/format.h"

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
      )
  {
    // read Brel
    u3shell::OperatorLabelsU3ST Brel_labels(-2,u3::SU3(0,2),0,0,0);
    Brel_sectors = u3shell::SectorsU3SPN(lsu3shell_space,Brel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(
        Brel_filename,
        lsu3shell_basis_table,lsu3shell_space,
        Brel_labels,Brel_sectors,Brel_matrices
      );

    // read Arel
    u3shell::OperatorLabelsU3ST Arel_labels(2,u3::SU3(2,0),0,0,0);
    Arel_sectors = u3shell::SectorsU3SPN(lsu3shell_space,Arel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(
        Arel_filename,
        lsu3shell_basis_table,lsu3shell_space,
        Arel_labels,Arel_sectors,Arel_matrices
      );

    // read Nrel
    u3shell::OperatorLabelsU3ST Nrel_labels(0,u3::SU3(0,0),0,0,0);
    Nrel_sectors = u3shell::SectorsU3SPN(lsu3shell_space,Nrel_labels,true);
    lsu3shell::ReadLSU3ShellRMEs(
        Nrel_filename,
        lsu3shell_basis_table,lsu3shell_space,
        Nrel_labels,Nrel_sectors,Nrel_matrices
      );

  }

  void
  ReadLSU3ShellSeedUnitTensorRMEs(
      const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
      const u3shell::SpaceU3SPN& lsu3shell_space, 
      const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels,
      const std::string& filename,
      u3shell::SectorsU3SPN& unit_tensor_sectors,
      basis::MatrixVector& unit_tensor_lsu3shell_matrices
    )
  {
    // lgi_unit_tensor_sectors.resize(lgi_unit_tensor_labels.size());
    // lgi_unit_tensor_lsu3shell_matrices.resize(lgi_unit_tensor_labels.size());
    // for (int unit_tensor_index=0; unit_tensor_index<lgi_unit_tensor_labels.size(); ++unit_tensor_index)
    //   {
    // set up aliases for current unit tensor
    
    // const u3shell::RelativeUnitTensorLabelsU3ST& unit_tensor_labels = lgi_unit_tensor_labels[unit_tensor_index];
    // u3shell::SectorsU3SPN& unit_tensor_sectors = lgi_unit_tensor_sectors[unit_tensor_index];
    // basis::MatrixVector& unit_tensor_lsu3shell_matrices = lgi_unit_tensor_lsu3shell_matrices[unit_tensor_index];
  
    // read rmes
    const bool spin_scalar = false;
    // std::string filename = fmt::format(relative_unit_tensor_filename_template,unit_tensor_index);
    unit_tensor_sectors = u3shell::SectorsU3SPN(lsu3shell_space,unit_tensor_labels,spin_scalar);
    lsu3shell::ReadLSU3ShellRMEs(
        filename,
        lsu3shell_basis_table,lsu3shell_space,
        unit_tensor_labels,unit_tensor_sectors,unit_tensor_lsu3shell_matrices
      );
      // }
  }


}  // namespace
