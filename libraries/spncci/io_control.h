/****************************************************************
  io_control.h

  High-level control code for data I/O in SpNCCI calculation programs.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).
  6/4/17 (mac): Revise to rescale rmes from relative to intrinsic operators
    on input
  6/16/17 (aem) : Extracted from spncci.cpp :ReadRelativeObservables 
  6/16/17 (aem) : Factored i/o for symplectic generators into function for
      Brel+Nrel and one for Arel
****************************************************************/

#ifndef SPNCCI_SPNCCI_IO_CONTROL_H_
#define SPNCCI_SPNCCI_IO_CONTROL_H_

#include <string>
#include <unordered_set>
#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "u3shell/relative_operator.h"
#include "u3shell/u3spn_scheme.h"
#include "u3shell/unit_tensor_space_u3s.h"
#include "u3shell/upcoupling.h"

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
  //   For <operator> in {B,N}:
  //     <operator>rel_filename (input): rme filename
  //     <operator>intr_sectors (output): U3SPN sectors
  //     <operator>intr_matrices (output): matrices of RMEs
  //   A (input): nucleon number


  void
    ReadLSU3ShellSymplecticRaisingOperatorRMEs(
        const lsu3shell::LSU3BasisTable& lsu3shell_basis_table,
        const u3shell::SpaceU3SPN& lsu3shell_space, 
        const std::string& Arel_filename, u3shell::SectorsU3SPN& Aintr_sectors, basis::MatrixVector& Aintr_matrices,
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
  //   For symplectic raising operator Arel:
  //     <operator>rel_filename (input): rme filename
  //     <operator>intr_sectors (output): U3SPN sectors
  //     <operator>intr_matrices (output): matrices of RMEs
  //   A (input): nucleon number



  void
  ReadRelativeObservables(
    int Nmax, int N1v, const std::vector<double>& hw_values,
    const std::string& observable_directory,const std::vector<std::string>& observable_filenames, 
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    std::vector<std::vector<u3shell::RelativeRMEsU3SSubspaces>>& observables_relative_rmes,
    std::vector<std::vector<u3shell::IndexedOperatorLabelsU3S>>& relative_observable_labels
    );
  // Reads in relative observable U(3)xSU(2)xSU(2) tensor components
  // in from files and generates a list of observable symmetries
  // (u3shell::IndexedOperatorLabelsU3S) for each observable.  Tensor
  // components are accumulated over all hw values (interactions may
  // have different tensor components depending on hw value).
  //
  //  Inputs:
  //    Nmax (int) : Oscillator truncation
  //    N1v (int) : Valence shell 
  //    hw_values : list of hw values in mesh
  //    observable_directory : location of observable files
  //    observable_filenames : vector of base names for observable files
  //    unit_tensor_space : set of unit tensor subspaces defined by the x0,S0,etap,eta unit tensor labels
  //
  //  Output:
  //    observables_relative_rmes : array of containers for observable rmes, indexed by hw, then by observable
  //    relative_observable_labels : vector of lists of (u3s,kappa0,L0) labels for each observable, 
  //    accumulated over all hw values.  Labels used to construct U3S sectors. 




}  // namespace

#endif
