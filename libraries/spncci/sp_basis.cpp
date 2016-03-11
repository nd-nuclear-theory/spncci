/****************************************************************
  sp_basis.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "spncci/sp_basis.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "utilities/parsing.h"


namespace spncci
{

  std::string LGI::Str() const
  {
    std::ostringstream ss;

    ss << "[" 
       << " " << Nex
       << " " << Sp
       << " " << Sn
       << " " << S
       << " " << sigma.Str()
       << " " << "]";
    return ss.str();
  }

  void GenerateLGIVector(std::vector<LGI>& basis, const std::string& lgi_filename, const HalfInt& Nsigma_0)
  {

    // open input file
    std::ifstream lgi_stream(lgi_filename.c_str());
    OpenCheck(lgi_stream,lgi_filename);

    // scan input file
    std::string line;
    int line_count = 0;
    while ( std::getline(lgi_stream, line) )
    {
      // count line
      ++line_count;

      // set up for parsing
      std::istringstream line_stream(line);

      // parse line
      //   Nex 2Sp 2Sn 2S lambda mu count
      int Nex, twice_Sp, twice_Sn, twice_S, lambda, mu, count;
      line_stream >> Nex >> twice_Sp >> twice_Sn >> twice_S >> lambda >> mu >> count;
      ParsingCheck (line_stream, line_count, line);

      // conversions
      HalfInt Sp = HalfInt(twice_Sp,2);
      HalfInt Sn = HalfInt(twice_Sn,2);
      HalfInt S = HalfInt(twice_S,2);
      HalfInt Nsigma = Nsigma_0 + Nex;
      u3::U3 sigma(Nsigma,u3::SU3(lambda,mu));
      
      // replicate LGI according to count
      for (int i=1; i<=count; ++i)
	basis.emplace_back(Nex,Sp,Sn,S,sigma);
    }

    // close input file
    lgi_stream.close();
  }


}  // namespace
