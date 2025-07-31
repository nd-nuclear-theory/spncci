/****************************************************************
  results_output.h

  Generates results for Sp(3,R) eigenproblem

  Colin V. Coane
  
  SPDX-License-Identifier: MIT 

  07/08/24 (cvc): Created.

****************************************************************/

#ifndef SP3RLIB_RESULTS_OUTPUT_H_
#define SP3RLIB_RESULTS_OUTPUT_H_

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>

#include "fmt/format.h"
#include "mcutils/eigen.h"
#include "mcutils/io.h"
#include "sp3rlib/sp3r.h"
#include "sp3rlib/sp3rsj.h"

#include "sp3rlib/config_parameters.h"


namespace sp3r
{


  void StartNewSection(std::ostream& out_stream, const std::string& title);
  // Start new file section.
  //
  // Arguments:
  //   out_stream (input): output stream
  //   title (input): section title


  // Write table of eigenvalues to output sorted by S,J
  void WriteEigenvalues(
    std::ostream& out_stream,
    const std::vector<std::pair<HalfInt,Halfint>>& SJ_values,
    const std::vector<sp3rlib::Vector>& eigenvalues
  );

  // Write eigenvectors for a given S,J to out_file stream. 
  void WriteEigenvectors(
    sp3r::Matrix& eigenvectors,
    const HalfInt& S,
    const HalfInt& J,
    std::ofstream& out_file,
    const int& binary_float_precision
  );

  // Write eigenvectors for a given S,J to file
  void WriteEigenvectors(
    sp3r::Matrix& eigenvectors,
    const HalfInt& S,
    const HalfInt& J,
    const std::string& filename,
    const int& binary_float_precision
  );

  void ReadEigenvectors(  
    const std::string& filename,
    sp3r::Matrix& eigenvectors
  );
  
}



#endif