/****************************************************************
  lgi_solver.h

  Interface for lsu3shell basis.

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.
  ...
  1/19/17 (mac): Replace null solver.
  2/17/17 (mac):
    - Remove internal function GenerateBrelNcmMatrices from header.
    - Remove leftover function declaration GenerateNcmOperatorBlocks<double>
      from header.
    - Change GenerateLGIExpansion to take matrices rather than
      streams.
    - Extract WriteLGILabels to lgi.h.
  6/11/17 (mac): Remove deprecated form of GenerateLGIExpansion with
    stream arguments.
  10/11/17 (aem) : Extract GetLGIExpansion from spncci/computation_control
  2/8/18 (aem) : removed lgi_families with gamma_max=0
  6/24/19 (aem) : Moved TransformOperatorTpSpBasis to lgi_unit_tensor
****************************************************************/

#ifndef LGI_SOLVER_H_
#define LGI_SOLVER_H_

#include <vector>
#include <Eigen/Eigen>

#include "am/am.h"
#include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_rme.h"

namespace lgi
{

  // output mode
  extern int binary_format_code;
  extern int binary_float_precision;
  typedef short unsigned int RMEIndexType;
  typedef uint32_t LGIIndexType;


  void
    GenerateLGIExpansion(
        const u3shell::SpaceU3SPN& space,
        const u3shell::SectorsU3SPN& Brel_sectors,
        const basis::OperatorBlocks<double>& Brel_matrices,
        const u3shell::SectorsU3SPN& Ncm_sectors,
        const basis::OperatorBlocks<double>& Ncm_matrices,
        HalfInt Nsigma_0,
        lgi::MultiplicityTaggedLGIVector& lgi_families,
        basis::OperatorBlocks<double>& lgi_expansions,
        std::vector<int>& lsu3shell_index_lookup_table
      );
  // Generate the LGI Expansion in terms of lsu3shell SU(3)xSU(2)
  // reduced basis states by solving for the null space of the
  // Brel+Ncm matrix.
  //
  //  lsu3_basis_table and space are generated as output from lsu3shell::ReadLSU3ShellBasis
  //  which reads in basis table generated by ncsmSU3xSU2BasisLSU3Tabular
  //
  // Construct Brel and Ncm matrix in lsu3shell basis and solve for null space.
  // Columns of kernel are expansion coefficients for each lgi.
  //
  //  Arguments:  (TODO update me!)
  //    A (input) : atomic mass number of nucleus
  //    lsu3_basis_table (input) : lookup table of LSU3Shell basis states between
  //                                LSU3Shell basis and U3SPN space.
  //    space (input) : space on which Brel and Nrel rme's are calculated
  //    Brel_filename (input) : file containing lsu3shell rme's of Brel
  //    Nrel_filename (input) : file containing lsu3shell rme's of Nrel
  //
  // lsu3hsell_index_lookup_table: look up take relating lsu3shell subspace index to lgi family index

void GetLGIExpansion(
    const u3shell::SpaceU3SPN& lsu3shell_space,
    const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
    const std::string& Brel_filename,
    const std::string& Nrel_filename,
    int A, HalfInt Nsigma_0,
    lgi::MultiplicityTaggedLGIVector& lgi_families,
    lsu3shell::OperatorBlocks& lgi_expansions,
    std::vector<int>& lsu3shell_index_lookup_table
  );
  // Get list of LGI labels and multiplicities and lgi expansions in lsu3shell basis
  //
  // Inputs
  //  lsu3shell_basis_table,lsu3shell_space,  Filenames and A
  //
  // Outputs
  //   lgi::MultiplicityTaggedLGIVector lgi_families;
  //   basis::OperatorBlocks<double> lgi_expansions;
  //   lsu3hsell_index_lookup_table: look up table relating
  //    lsu3shell subspace index to lgi family index needed
  //    when there are lsu3shell subspaces with no lgi

void WriteLSU3ShellToLGIConversionTable(const std::vector<int>& lsu3shell_index_lookup_table);
  // Write out table that converts lsu3shell subspace index to lgi index

void ReadLSU3ShellToLGIConversionTable(std::vector<int>& lsu3shell_index_lookup_table);
  // Read in table that converts lsu3shell subspace index to lgi index



void WriteLGIExpansions(
    const std::string& filename,
    const lsu3shell::OperatorBlocks& lgi_expansions
  );
  // Write expansion of lgi's in lsu3shell basis as su3 reduced level
  // Inputs
  //  filename: output filename
  //  lgi_expansion : operator blocks containing lgi expansion

void WriteLGIExpansions(
  const std::string& filename,
  const lsu3shell::OperatorBlock& lgi_expansion
  );
  // Writes expansion of single LGI subspace in terms of lsu3shell basis at su3 reduced level
  //  File is binary.
  // Inputs:
  //  filename : output filename
  //  lgi_expansion : matrix containing LGI expansion.
  //      Columns represent different LGI, rows are U(3)S basis states

void WriteLGIExpansionsText(
  const std::string& filename,
  const lsu3shell::OperatorBlock& lgi_expansion
  );

  // Writes expansion of single LGI subspace in terms of lsu3shell basis at su3 reduced level
  //  File is text.
  // Inputs:
  //  filename : output filename
  //  lgi_expansion : matrix containing LGI expansion.
  //      Columns represent different LGI, rows are U(3)S basis states


void ReadLGIExpansion(
    int num_lgi_subspaces,
    const std::string& filename,
    basis::OperatorBlocks<double>& lgi_expansions
  );
  // Read LGI expansion from file
  // Inputs:
  //  num_lgi_subspaces
  //  filename
  //
  // Outputs
  //  lgi_expansions

}
#endif
