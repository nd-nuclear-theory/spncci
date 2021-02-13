/****************************************************************
  computation_control.h

  High-level control code for computations in SpNCCI calculation
  programs.
                                  
  DEPRECATED CODE STRUCTURE -- instead split up routines
  thematically...  rename to "recurrence.h"???

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT

  02/19/17 (mac): Extracted from unit tensor test codes
    (compute_unit_tensor_rmes.cpp and explicit.cpp).
  06/16/17 (aem) : Extracted from spncci.cpp
      + GetUnitTensorSeedBlocks
      + GetLGIExpansion
      + GenerateRecurenceUnitTensors
      + GenerateRecurrenceHypersectors
      + PopulateHypersectorsWithSeeds
  07/01/17 (aem) : Updated GetLGIExpansion for intrinsic basis
  07/06/17 (mac): Move out branching routines and eigenproblem control.
  10/04/17 (aem) : Updated CheckUnitTensorRecurrence to work for A<6
  10/11/17 (aem) : Moved GetLGIExpansion to lgi_solver
  01/30/18 (aem) :
      + Moved GenerateRecurrenceUnitTensor to spncci/recurrence
      + Eliminated GenerateRecurrenceHypersectors and moved
        hypersector construction to spncci/spncci_basis
      + Moved PopulateHypersectorsWithSeeds to spncci/recurrence
      + Elimiated GetUnitTensorSeedBlocks. Functionality moved to
        stand-alone programs/lgi/get_spncci_seed_blocks
  08/26/20 (aem) : Update basis setup and removed dependencies on
        spncci branched basis classes (SpaceSpU3S, SpaceSpLS,SpaceSpJ)
  01/24/21 (aem) : Deleted obselete SetUpSpNCCISpaces function using
        SpaceSpU3S, SpaceSpLS,SpaceSpJ
 
****************************************************************/

#ifndef SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_
#define SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_

#include <fstream>

#include "lgi/lgi.h"
#include "spncci/parameters.h"
#include "spncci/spncci_basis.h"
#include "spncci/parameters.h"
#include "spncci/vcs_cache.h"


namespace spncci
{
  void InitializeSpNCCI();

  void SetUpSpNCCISpaces(
      spncci::RunParameters& run_parameters,
      lgi::MultiplicityTaggedLGIVector& lgi_families,
      spncci::SpNCCISpace& spncci_space,
      spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::BabySpNCCISpace& baby_spncci_space,
      std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
      bool verbose=false
    );
  // Sets up spaces used in spncci calculations.  
  // input:
  //    run_parameters : Code run parameters including Nmax, Nsigma0, A and J_values
  //    verbose (optional): If true, print basis str (default false)
  //
  // output:
  //    lgi_families : List of lgi labels and associated multiplicities read in from file 
  //    spncci_space : 
  //    sigma_irrep_map :
  //    baby_spncci_space :
  //    spaces_spbasis : 


  void SetUpSpNCCISpaces(
      spncci::RunParameters& run_parameters,
      lgi::MultiplicityTaggedLGIVector& lgi_families,
      spncci::SpNCCISpace& spncci_space,
      spncci::SigmaIrrepMap& sigma_irrep_map,
      spncci::BabySpNCCISpace& baby_spncci_space,
      std::vector<spncci::SpaceSpBasis>& spaces_spbasis,
      spncci::KMatrixCache& k_matrix_cache, 
      spncci::KMatrixCache& kinv_matrix_cache,
      bool verbose=false
    );
    // Sets up spaces used in spncci calculations and computes K matrices for basis orthogonalization
    // 

}  // namespace

#endif
