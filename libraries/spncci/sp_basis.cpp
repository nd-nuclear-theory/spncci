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

  void GenerateLGIVector(std::vector<spncci::LGI>& lgi_vec, const std::string& lgi_filename, const HalfInt& Nsigma_0)
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
        //   Nex lambda mu 2Sp 2Sn 2S count
        int Nex, twice_Sp, twice_Sn, twice_S, lambda, mu, count;
        line_stream >> Nex >> lambda >> mu >> twice_Sp >> twice_Sn >> twice_S >> count;
        ParsingCheck (line_stream, line_count, line);

        // conversions
        HalfInt Nsigma = Nsigma_0 + Nex;
        u3::U3 sigma(Nsigma,u3::SU3(lambda,mu));
        HalfInt Sp = HalfInt(twice_Sp,2);
        HalfInt Sn = HalfInt(twice_Sn,2);
        HalfInt S = HalfInt(twice_S,2);
      
        // replicate LGI according to count
        for (int i=1; i<=count; ++i)
          lgi_vec.emplace_back(Nex,sigma,Sp,Sn,S);
      }

    // close input file
    lgi_stream.close();
  }

  int NmaxTruncator::Nn_max(const u3::U3S& sigmaS) const
  {
    HalfInt Nsigma = sigmaS.U3().N();
    int value = Nmax_ - int(Nsigma - Nsigma_0_);
    return value;
  }

  void GenerateSp3RSpaces(
                          const std::vector<spncci::LGI>& lgi_vec,
                          std::map<u3::U3S,sp3r::Sp3RSpace>& sigma_map,
                          const TruncatorInterface& truncator
                          )
  {
    // traverse LGI vector
    for (auto it = lgi_vec.begin(); it != lgi_vec.end(); ++it)
      {
        // find sigmaS of LGI
        u3::U3S sigmaS(it->sigma,it->S);
        
        // add sigmaS to map if needed
        if (sigma_map.count(sigmaS) == 0)
          {
            int Nn_max = truncator.Nn_max(sigmaS);
            // sigma_map[sigmaS] = sp3r::Sp3RSpace(sigmaS,Nn_max);  // FAILS: since compiler requires default constructor in case key is not present 
            // FAILS in g++ 4.5 -- map.emplace not yet defined?
            // sigma_map.emplace(
            //                sigmaS,  // key
            //                sigmaS,Nn_max  // constructor arguments for value
            //                );
            sp3r::Sp3RSpace sp3r_space(sigmaS.U3(),Nn_max);
            sigma_map.insert(std::make_pair(sigmaS,sp3r_space));

          }
      }
  }


}  // namespace
