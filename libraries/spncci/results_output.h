/****************************************************************
  results_output.h

  Code to generate results tabulations.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/17/17 (mac): Created.
****************************************************************/

#ifndef SPNCCI_SPNCCI_RESULTS_OUTPUT_H_
#define SPNCCI_SPNCCI_RESULTS_OUTPUT_H_

#include <iostream>
#include <string>

#include "spncci/branching.h"
#include "spncci/parameters.h"

namespace spncci
{

  void WriteResultsHeader(std::ostream& out_stream, const spncci::RunParameters& run_parameters);

  void WriteResultsBasis(
      std::ostream& out_stream,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceSpU3S& spu3s_space,
      const spncci::SpaceSpLS spls_space
    );

  void WriteResultsEigenvalues(
      std::ostream& out_stream
      // TODO
    );

}  // namespace

#endif
