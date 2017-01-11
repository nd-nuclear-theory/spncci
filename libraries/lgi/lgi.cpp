/****************************************************************
  lgi.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "lgi/lgi.h"

#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"
#include "cppformat/format.h"


namespace lgi
{

  std::string LGI::Str() const
  {
    std::ostringstream ss;

    ss << "[" 
       << " " << Nex()
       << " " << U3().Str()
       << " (" << Sp()
       << "," << Sn()
       << "," << S()
       << ") " << "]";
    return ss.str();
  }

  void ReadLGISet(LGIVector& lgi_vector, const std::string& lgi_filename)
  {
    // open input file
    std::ifstream lgi_stream(lgi_filename.c_str());
    OpenCheck(bool(lgi_stream),lgi_filename);  // explicit cast to bool required in gcc 5

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
        //   Nex 2N lambda mu 2Sp 2Sn 2S count
        int Nex, twice_N, twice_Sp, twice_Sn, twice_S, lambda, mu, count;
        line_stream >> Nex
                    >> twice_N  >> lambda >> mu >> twice_Sp >> twice_Sn >> twice_S
                    >> count;
        ParsingCheck(line_stream, line_count, line);
        // conversions
        HalfInt Nsigma = HalfInt(twice_N,2);
        // std::cout<<fmt::format("{} {} {}", Nsigma, lambda,mu)<<std::endl;
        // assert(Nsigma == Nsigma_0 + Nex);
        u3::U3 sigma(Nsigma,u3::SU3(lambda,mu));
        HalfInt Sp = HalfInt(twice_Sp,2);
        HalfInt Sn = HalfInt(twice_Sn,2);
        HalfInt S = HalfInt(twice_S,2);
        u3shell::U3SPN omegaSPN(u3::U3S(sigma,S),Sp,Sn);
      
        // replicate LGI according to count
        lgi_vector.emplace_back(lgi::LGI(omegaSPN,Nex),count);
      }

    // close input file
    lgi_stream.close();
  }
}
