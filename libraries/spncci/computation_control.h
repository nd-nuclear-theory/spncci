/****************************************************************
  computation_control.h

  High-level control code for computations in SpNCCI calculation
  programs.
                                  
  DEPRECATED CODE STRUCTURE -- instead split up routines
  thematically...  rename to "recurrence.h"???

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  2/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).
  6/16/17 (aem) : Extracted from spncci.cpp
      + GetUnitTensorSeedBlocks
      + GetLGIExpansion
      + GenerateRecurenceUnitTensors
      + GenerateRecurrenceHypersectors
      + PopulateHypersectorsWithSeeds
  7/1/17 (aem) : Updated GetLGIExpansion for intrinsic basis
  7/6/17 (mac): Move out branching routines and eigenproblem control.
  10/4/17 (aem) : Updated CheckUnitTensorRecurrence to work for A<6
****************************************************************/

#ifndef SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_
#define SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_

#include "lgi/lgi.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "lsu3shell/lsu3shell_rme.h"
#include "spncci/branching_u3s.h"
#include "spncci/branching_u3lsj.h"
#include "spncci/branching.h"
#include "spncci/explicit_construction.h"
#include "spncci/io_control.h"
#include "spncci/spncci_basis.h"
#include "spncci/spncci_common.h"
#include "spncci/unit_tensor.h"
#include "sp3rlib/u3.h"
#include "u3shell/relative_operator.h"
#include "u3shell/u3spn_scheme.h"

namespace spncci
{

  // convenience typedef for use in iteration over J sectors
  typedef std::pair<HalfInt,HalfInt> JPair;


void GetLGIExpansion(
    const u3shell::SpaceU3SPN& lsu3shell_space, 
    const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
    const std::string& Brel_filename,
    const std::string& Nrel_filename,
    int A, HalfInt Nsigma_0,
    lgi::MultiplicityTaggedLGIVector& lgi_families,
    basis::MatrixVector& lgi_expansions
  );
  // Get list of LGI labels, multiplicities and expansion in lsu3shell basis
  //
  // Inputs
  //  lsu3shell_basis_table,lsu3shell_space,  Filenames and A
  // 
  // Outputs 
  //   lgi::MultiplicityTaggedLGIVector lgi_families;
  //    basis::MatrixVector lgi_expansions;

  void 
  GetUnitTensorSeedBlocks(
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const std::string& relative_unit_tensor_filename_template,
    const u3shell::SpaceU3SPN& lsu3shell_space, 
    const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
    const basis::MatrixVector& lgi_expansions,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    std::map< std::pair<int,int>, std::map<std::pair<int,int>, basis::OperatorBlocks<double>>>& lgi_unit_tensor_blocks
  );
  //  For each unit tensor:
  //    + Determine corresponding unit tensor subspace
  //    + Read in lsu3shell rmes
  //    + Transform lsu3shell rmes to spncci basis.  
  //
  //  The spncci unit tensor blocks are stored in lgi_unit_tensor_blocks.  
  //    - Outer map is keyed by irrep family pair 
  //        <irrep_family_index_1, irrep_family_index_2>
  //    - Inner map is keyed by a multiplicity tagged unit tensor subspace 
  //        <unit_tensor_subspace_index,outer_multipicity>
  //
  //  Inputs:
  //    lgi_unit_tensor_labels : list of unit tensor labels.  Order of labels expected to correspond to order of 
  //      su3rmes of operators 
  //    unit_tensor_space : Unit tensors broken up into subspaces by (x0,S0,etap,eta) labels
  //    relative_unit_tensor_filename_template : template of file name containing su3rmes for each unit tensor
  //    lsu3shell_space : basis indexing for lsu3shell basis
  //    lsu3shell_basis_table : provides lookup for converting between our indexing of lsu3shell and internal lsu3shell basis
  //      indexing
  //    lgi_expansions : expansion of each lgi in terms of lsu3shell SU(3) irreps
  //    baby_spncci_space : basis indexing for Spncci basis branched to U(3). 
  //
  //  Outputs:
  //    lgi_unit_tensor_blocks : container for rmes of unit tensors between lgi (seeds for unit tensor recurrence).
  //


void
  GenerateRecurrenceUnitTensors(
      int Nmax,
      const std::set<int>& lgi_operator_subset,
      const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
      std::map<spncci::NnPair,std::set<int>>& operator_subsets
    );
  // Creates a list of indices of unit_tensor_subspaces in unit_tensor_space
  // which will appear in the spncci recurrence for a given pair of symplectic irreps.
  //
  //  Inputs:
  //    Nmax : usual basis truncation parameters
  //    lgi_operator_subset : a list of unit tensor subspace indices
  //        with non-zero rme's between the lgis
  //    unit_tensor_space : space organizing unit tensors by symmetry labels 
  //   Output: 
  //    operator_subsets : list of indices of unit tensor subspaces, organized by
  //      


void GenerateRecurrenceHypersectors(
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const std::map<std::pair<int,int>,std::set<int>>& lgi_unit_tensor_subset,
    int Nmax, int irrep_family_index_bra, int irrep_family_index_ket,
    std::vector<std::vector<int>>& unit_tensor_hypersector_subsets,
    spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors

  );
  // Generate set of baby_spnnci_hypersectors for recurrence for given lgi pair
  // Also generates sets of unit tensor hypersector indices grouped by Nsum=Nnp+Nn
  //
  // Inputs:
  //    unit_tensor_space : Unit tensors broken up into subspaces by (x0,S0,etap,eta) labels
  //    baby_spncci_space : Indexing for spncci U(3)xSU(2) irreps
  //    lgi_unit_tensor_subset : set of unit tensor subspaces corresponding to lgis
  //    ...
  //
  // Outputs:
  //    unit_tensor_hypersector_subsets :
  //    baby_spncci_hypersectors :  


void PopulateHypersectorsWithSeeds(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const std::vector<int>& unit_tensor_hypersector_subset,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    const std::map< std::pair<int,int>, std::map<std::pair<int,int>, basis::OperatorBlocks<double>>>& 
      lgi_unit_tensor_blocks,
    basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  );
  // Populate unit_tensor_hyperblocks with seeds obtained from lsu3shell 
  //
  //  Inputs:
  //    irrep_family_index_bra, irrep_family_index_ket : indices of irrep family in spncci
  //    unit_tensor_hypersector_subsets[0] : subset of hypersectors for given irrep family pair for Nex=0
  //    baby_spncci_space : index spncci basis
  //    baby_spncci_hypersectors : hypersectors between given SU(3)xSU(2) irreps in irrep families 
  //    lgi_unit_tensor_blocks : container holding seeds to put into unit_tensor_hyperblocks
  //
  //  outputs:
  //    unit_tensor_hyperblocks : container holding rmes of unit tensor hypersectors (to be populated with seeds)
  //

  void CheckUnitTensorRecurrence(
    int irrep_family_index_bra, int irrep_family_index_ket,
    const u3shell::RelativeUnitTensorSpaceU3S& unit_tensor_space,
    const std::vector<u3shell::RelativeUnitTensorLabelsU3ST>& lgi_unit_tensor_labels,
    const std::string& relative_unit_tensor_filename_template,
    const u3shell::SpaceU3SPN& lsu3shell_space, 
    const lsu3shell::LSU3ShellBasisTable& lsu3shell_basis_table,
    const spncci::SpNCCISpace& spncci_space,
    const spncci::BabySpNCCISpace& baby_spncci_space,
    const basis::MatrixVector& spncci_expansions,
    const spncci::BabySpNCCIHypersectors& baby_spncci_hypersectors,
    const basis::OperatorHyperblocks<double>& unit_tensor_hyperblocks
  );
  // Computes unit tensor hyperblocks in spncci basis which is constructed explicitly in terms of lsu3shell irreps
  // and compares unit tensor hyperblocks computed recursively to explicitly computed unit tensor hyperblocks 
  // If difference between hyperblocks is found exceeding tolerance of 1e-4, print error message
  // Otherwise, print "no errors".
  //
  // Inputs:
  //    
  //
  // Outputs:

  void 
  SolveHamiltonian(
      const spncci::OperatorBlock& hamiltonian_matrix,
      const HalfInt& J,
      int num_eigenvalues,
      int eigensolver_num_convergence,  // whatever exactly this is...
      int eigensolver_max_iterations,
      double eigensolver_tolerance,
      spncci::Vector& eigenvalues,  // eigenvalues for J-subspace
      spncci::Matrix& eigenvectors  // eigenvectors for J-subspace
    );
  // Solve the Hamiltonian matrix (in a single J-space) for energy eigenvalues and vectors 


}  // namespace

#endif
