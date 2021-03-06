/****************************************************************
  results_output.h

  Code to generate results tabulations.
                                  
  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

  6/17/17 (mac): Created.
  6/24/17 (mac): Add SpU3S subspace listing.
  6/27/17 (mac): Suppress output for empty subspaces or sectors.
****************************************************************/

#ifndef SPNCCI_SPNCCI_RESULTS_OUTPUT_H_
#define SPNCCI_SPNCCI_RESULTS_OUTPUT_H_

#include <iostream>
#include <string>

#include "cppformat/format.h"
#include "spncci/branching.h"
#include "spncci/parameters.h"

namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // version stamp
  ////////////////////////////////////////////////////////////////
  const int g_results_version = 201706270; // "yyyymmddv" (v=version w/in day)

  ////////////////////////////////////////////////////////////////
  // output utilities
  ////////////////////////////////////////////////////////////////

  void StartNewSection(std::ostream& out_stream, const std::string& title);
  // Start new file section.
  //
  // Arguments:
  //   out_stream (input): output stream
  //   title (input): section title

  template <typename tValue>
    void WriteKeyValue(std::ostream& out_stream, const std::string& keyword, const std::string& format, tValue value)
  // Write key-value pair.
  //
  // Arguments:
  //   out_stream (input): output stream
  //   title (input): keyword
  //   format (input): format code for value (no braces!)
  //   value (input): value to write
  {
    std::string full_format = fmt::format("{{{}}}",format);  // encapsulate format code in braces
    out_stream << fmt::format("{} = {}",keyword,fmt::format(full_format,value)) << std::endl;
  }

  template <typename tValueIterable>
    void WriteKeyValueList(
        std::ostream& out_stream, const std::string& keyword, const std::string& format,
        tValueIterable data
      )
    // Write key-value-list pair.
    //
    // The data may be any iterable container of values, such as std::vector or std::array.
    //
    // Arguments:
    //   out_stream (input): output stream
    //   title (input): keyword
    //   format (input): format code for value (no braces!)
    //   date (input): values to write
    {
    std::string full_format = fmt::format("{{{}}}",format);  // encapsulate format code in braces
    out_stream << fmt::format("{} =",keyword);
    for (const auto& value : data)
      out_stream << " " << fmt::format(full_format,value);
    out_stream << std::endl;
      
  }


  ////////////////////////////////////////////////////////////////
  // output code
  ////////////////////////////////////////////////////////////////

  void WriteCodeInformation(std::ostream& out_stream, const spncci::RunParameters& run_parameters);

  void WriteRunParameters(std::ostream& out_stream, const spncci::RunParameters& run_parameters);

  void WriteBasisStatistics(
      std::ostream& out_stream,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceSpU3S& spu3s_space,
      const spncci::SpaceSpLS& spls_space,
      const spncci::SpaceSpJ& spj_space
    );

  void WriteSpU3SSubspaceListing(
      std::ostream& out_stream,
      const spncci::SpaceSpU3S& spu3s_space,
      HalfInt Nsigma0
    );
  // Write listing of SpU3S basis subspaces.
  //
  // This is technically superfluous, since all information on SpU3S
  // subspaces can be deduced from the listing of BabySpNCCI
  // subspaces.  Therefore, reliance on this listing is discouraged.
  //
  // Arguments:
  //     out_stream (input): output stream
  //     spu3s_space (input): SpU3S basis
  //     Nsigma0 (input): U(1) label of lowest configuration

  void WriteBabySpNCCISubspaceListing(
      std::ostream& out_stream,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      HalfInt Nsigma0
    );

  void WriteCalculationParameters(
      std::ostream& out_stream,
      double hw
    );

  void WriteEigenvalues(
      std::ostream& out_stream,
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Vector>& eigenvalues,
      int gex
    );
  // Write table of eigenvalues.
  //
  // CAVEAT: parity quantum number is now simply taken as a given from
  // parameter rather than determined in some way from space
  //
  // Arguments:
  //   out_stream (input): output stream
  //   ...

  void WriteDecompositions(
      std::ostream& out_stream,
      const std::string& decomposition_name,
      const std::string& format_string,
      const spncci::SpaceSpJ& spj_space,
      const std::vector<spncci::Matrix>& decompositions,
      int gex
    );
  // Write tables of eigenfunction probability distributions.
  //
  // CAVEAT: For now, parity quantum number is now simply taken as a
  // given from parameter rather than determined in some way from the
  // SpJ subspaces.
  //
  // Arguments:
  //   out_stream (input): output stream
  //   decomposition_name (input): decompsition name for use in res file section title
  //   format_string (input): floating point format string for matrix entries
  //   spj_space (input): SpJ basis (for use in retrieving J value for use in results
  //     file comment line)
  //   decompositions (input): the decompositions, by J subspace
  //   gex (input): excitation parity grade (for use in results file comment line)  

  void WriteObservables(
      std::ostream& out_stream,
      const std::vector<spncci::SectorsSpJ>& observable_sectors,
      const std::vector<spncci::OperatorBlocks>& observable_results_matrices,
      int gex
    );
  // Write matrices of RMEs for observables.
  //
  // CAVEAT: For now, parity quantum number is taken as a given from
  // parameter rather than determined in some way from the bra and ket
  // SpJ subspaces.
  //
  // Arguments:
  //   out_stream (input): output stream
  //   ...

}  // namespace

#endif
