/****************************************************************
  lgi_solver.h

  Interface for lsu3shell basis.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  8/1/16 (aem,mac): Created.
  9/7/16 (mac): Split from lsu3shell_interface.

****************************************************************/

#ifndef LGI_SOLVER_H_
#define LGI_SOLVER_H_

#include <boost/functional/hash_fwd.hpp>
#include <eigen3/Eigen/Eigen>

#include "am/am.h"  
#include "sp3rlib/sp3r.h"
#include "u3shell/relative_operator.h"
#include "u3shell/tensor_labels.h"
#include "u3shell/two_body_operator.h"
#include "u3shell/u3spn_scheme.h"

namespace lgi
{

  void GenerateLSU3ShellExpansionLGI(
      int Nsigma_0,
      int Nsigma_min,
      int Nsigma_max, 
      std::string basis_file, 
      std::string brel_filename, 
      std::string nrel_filename,
      std::string lgi_filename,
      std::string lgi_expansion_filename
    );

}
#endif
