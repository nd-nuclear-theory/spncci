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
  10/11/17 (aem) : Moved GetLGIExpansion to lgi_solver
  1/30/18 (aem) :
      + Moved GenerateRecurrenceUnitTensor to spncci/recurrence
      + Eliminated GenerateRecurrenceHypersectors and moved
        hypersector construction to spncci/spncci_basis
      + Moved PopulateHypersectorsWithSeeds to spncci/recurrence
      + Elimiated GetUnitTensorSeedBlocks. Functionality moved to
        stand-alone programs/lgi/get_spncci_seed_blocks
****************************************************************/

#ifndef SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_
#define SPNCCI_SPNCCI_COMPUTATION_CONTROL_H_

// #include "lgi/lgi.h"
// #include "lsu3shell/lsu3shell_basis.h"
// #include "lsu3shell/lsu3shell_rme.h"
// #include "spncci/branching_u3s.h"
// #include "spncci/branching_u3lsj.h"
// #include "spncci/branching.h"
// #include "spncci/explicit_construction.h"
// #include "spncci/io_control.h"
// #include "spncci/spncci_basis.h"
// #include "spncci/spncci_common.h"
// #include "spncci/parameters.h"
// #include "spncci/unit_tensor.h"
// #include "sp3rlib/u3.h"
// #include "u3shell/relative_operator.h"
// #include "u3shell/u3spn_scheme.h"

namespace spncci
{}  // namespace

#endif
