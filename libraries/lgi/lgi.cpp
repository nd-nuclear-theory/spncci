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

  void 
  WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_families,std::ostream& os)
  {
    int Nex;
    u3::U3 sigma;
    HalfInt Sp,Sn,S;
    // std::unordered_map<lgi::LGI,int,boost::hash<lgi::LGI>> lgi_counter;
    for(auto lgi_count : lgi_families)
      {
        std::tie(Nex,sigma,Sp,Sn,S)=lgi_count.irrep.Key();
        int count=lgi_count.tag;
        os
          <<Nex
          <<"  "<<TwiceValue(sigma.N())<<"  "<<sigma.SU3().lambda()<<"  "<<sigma.SU3().mu()
          <<"  "<<TwiceValue(Sp)<<"  "<<TwiceValue(Sn)<<"  "<<TwiceValue(S)
          <<"  "<<count<<std::endl;     
      }
    }

  void ReadLGISet(MultiplicityTaggedLGIVector& lgi_vector, const std::string& lgi_filename)
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
