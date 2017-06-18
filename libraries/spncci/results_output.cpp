/****************************************************************
  io_control.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/results_output.h"

#include "cppformat/format.h"
#include "spncci/parameters.h"
#include "spncci/spncci_common.h"


namespace spncci
{

  ////////////////////////////////////////////////////////////////
  // output utilities (for internal use)
  ////////////////////////////////////////////////////////////////

  void StartNewSection(std::ostream& out_stream, const std::string& title)
  // Start new file section.
  //
  // Arguments:
  //   out_stream (input): output stream
  //   title (input): section title
  {
    out_stream << std::endl;
    out_stream << fmt::format("[{}]",title) << std::endl;
  }

  template <typename tValue>
  void WriteKeyValue(std::ostream& out_stream, const std::string& keyword, const std::string& format, tValue value)
  // Start new file section.
  //
  // Arguments:
  //   out_stream (input): output stream
  //   title (input): keyword
  //   format (input): format code for value (no braces!)
  //   value (input): value to write
  {
    // in principle, can instead use nested format string, but syntax is unclear from cppformat docs
    out_stream << fmt::format("{} = {"+format+"}",keyword,value) << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  // output code
  ////////////////////////////////////////////////////////////////

  void WriteResultsHeader(std::ostream& out_stream, const spncci::RunParameters& run_parameters)
  {
    StartNewSection(out_stream,"CODE");
    // StartNewSection(out_stream,"Internals");
    // WriteKeyValue(out_stream,"g_zero_tolerance",":e",g_zero_tolerance);

    StartNewSection(out_stream,"PARAMETERS");

    StartNewSection(out_stream,"Basis");
    WriteKeyValue(out_stream,"A",":d",run_parameters.A);
    WriteKeyValue(out_stream,"Nsigma0",":.1f",float(run_parameters.Nsigma0));
    WriteKeyValue(out_stream,"Nsigmamax",":d",run_parameters.Nsigmamax);
    WriteKeyValue(out_stream,"N1v",":d",run_parameters.N1v);
    WriteKeyValue(out_stream,"Nmax",":d",run_parameters.Nmax);

  }

  void WriteResultsBasis(
      std::ostream& out_stream,
      const spncci::SpNCCISpace& spncci_space,
      const spncci::BabySpNCCISpace& baby_spncci_space,
      const spncci::SpaceSpU3S& spu3s_space,
      const spncci::SpaceSpLS spls_space
    )
  {
    StartNewSection(out_stream,"BASIS");

    StartNewSection(out_stream,"Irreps");

    StartNewSection(out_stream,"BabySpNCCI");

    StartNewSection(out_stream,"SpU3S");
    WriteKeyValue(out_stream,"subspaces",":d",spu3s_space.size());
    WriteKeyValue(out_stream,"dimension",":d",spu3s_space.Dimension());
    WriteKeyValue(out_stream,"full_dimension",":d",spu3s_space.FullDimension());

    StartNewSection(out_stream,"SpLS");

  }


}  // namespace
