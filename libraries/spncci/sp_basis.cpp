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

  void GenerateLGIVector(LGIVectorType& lgi_vector, const std::string& lgi_filename, const HalfInt& Nsigma_0)
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
        ParsingCheck(line_stream, line_count, line);

        // conversions
        HalfInt Nsigma = Nsigma_0 + Nex;
        u3::U3 sigma(Nsigma,u3::SU3(lambda,mu));
        HalfInt Sp = HalfInt(twice_Sp,2);
        HalfInt Sn = HalfInt(twice_Sn,2);
        HalfInt S = HalfInt(twice_S,2);
      
        // replicate LGI according to count
        for (int i=1; i<=count; ++i)
          lgi_vector.emplace_back(Nex,sigma,Sp,Sn,S);
      }

    // close input file
    lgi_stream.close();
  }

  int NmaxTruncator::Nn_max(const u3::U3& sigma) const
  {
    HalfInt Nsigma = sigma.N();
    int value = Nmax_ - int(Nsigma - Nsigma_0_);
    return value;
  }

  void GenerateSp3RIrreps(
                          LGIVectorType& lgi_vector,
                          SigmaIrrepMapType& sigma_irrep_map,
                          const TruncatorInterface& truncator
                          )
  {
    // traverse LGI vector
    for (auto it = lgi_vector.begin(); it != lgi_vector.end(); ++it)
      {
        // retrieve sigma of LGI
        // u3::U3 sigma(it->sigma); // CRYPTIC
        spncci::LGI& lgi = *it;
        u3::U3 sigma = lgi.sigma;

        // find truncation
        //
        // Note: Even if a subspace is truncated into oblivion (Nn_max<0), we
        // still need to create an empty Sp3RSpace and store the
        // relevant information with the LGI, e.g., that its dimension
        // is zero.

        int Nn_max = truncator.Nn_max(sigma);

        // construct and add Sp3RSpace (if needed)

        // FAILS: since compiler anticipates the case where sigma
        // is not found as a key and map would (hypothetically)
        // need to insert an SP3RSpace using the (nonexistent)
        // default constructor sp3r::Sp3RSpace()
        //
        // sigma_irrep_map[sigma] = sp3r::Sp3RSpace(sigma,Nn_max);  

        // FAILS in g++ 4.5: map.emplace not yet defined?
        //
        // sigma_irrep_map.emplace(
        //                sigma,  // key
        //                sigma,Nn_max  // constructor arguments for value
        //                );

        if (sigma_irrep_map.count(sigma) == 0)
          sigma_irrep_map.insert(
                                 //                                 std::make_pair(sigma,sp3r::Sp3RSpace(sigma,Nn_max))
                                 );

        // save info back to LGI
        
        int dimension = 0;
        
        //lgi.SaveSubspaceInfo(Nn_max,dimension,sigma_irrep_map[sigma]);

      }
  }


}  // namespace
