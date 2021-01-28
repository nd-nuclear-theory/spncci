/****************************************************************
  lgi.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame

****************************************************************/
#include "lgi/lgi.h"

#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"
#include "mcutils/eigen.h"
#include "fmt/format.h"


namespace lgi
{

  HalfInt Nsigma0ForNuclide(const NuclideType& nuclide, bool intrinsic)
  {
    // each major shell eta=2*n+l (for a spin-1/2 fermion) contains (eta+1)*(eta+2) substates
    HalfInt Nsigma0 = 0;
    for (int species_index : {0,1})
      {
        int num_particles = nuclide[species_index];
        for (int eta=0; num_particles>0; ++eta)
          {
            // add contribution from particles in shell
            int shell_degeneracy = (eta+1)*(eta+2);
            int num_particles_in_shell = std::min(num_particles,shell_degeneracy);
            // want num_particles_in_shell*(eta+HalfInt(3,2)), but HalfInt does not provide multiplication
            Nsigma0 += HalfInt(num_particles_in_shell*(2*eta+3),2);

            // discard particles in shell
            num_particles -= num_particles_in_shell;
          }
      }
    // If intrinsic remove cm zero point energy 3/2
    if(intrinsic)
      Nsigma0=Nsigma0-HalfInt(3,2);
    return Nsigma0;
  }


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

  std::string LGIOutputStr(const MultiplicityTagged<lgi::LGI>& lgi_family)
    {
      int Nex;
      u3::U3 sigma;
      HalfInt Sp,Sn,S;

      std::tie(Nex,sigma,Sp,Sn,S)=lgi_family.irrep.Key();
      int gamma_max=lgi_family.tag;

      std::ostringstream ss;
      ss<<fmt::format("{}  {}  {}  {}  {}  {}  {}",
          Nex,sigma.SU3().lambda(), sigma.SU3().mu(),
          TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S),gamma_max);
      return ss.str();
    }

  void 
  WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_families, std::ostream& os)
  {
    int Nex;
    u3::U3 sigma;
    HalfInt Sp,Sn,S;
    // std::unordered_map<lgi::LGI,int,boost::hash<lgi::LGI>> lgi_counter;
    for(auto lgi_count : lgi_families)
      {
        std::tie(Nex,sigma,Sp,Sn,S)=lgi_count.irrep.Key();
        int count=lgi_count.tag;
        os<<LGIOutputStr(lgi_count)<<std::endl;
      }
  }

  void 
  WriteLGILabels(const lgi::MultiplicityTaggedLGIVector& lgi_families, const std::string& filename)
  {
    std::ofstream os;
    os.open(filename);
    WriteLGILabels(lgi_families, os); 
    os.close();
  }

  void 
  WriteLGIExpansionHeader(int Z, int N, int Nmax, std::ostream& os)
    {
      os<<"# Z  N  Nmax"<<std::endl;
      os<<"#  Nex lambda mu 2Sp 2Sn 2S gamma_max"<<std::endl;
      os<<fmt::format("{}  {}  {}", Z, N, Nmax)<<std::endl;
    }

  void 
  WriteLGIExpansion(
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    const lsu3shell::OperatorBlocks&lgi_expansions,
    std::ostream& os
  )
  {
    for(int i=0; i<lgi_families.size(); ++i)
      {
        const auto& lgi_family=lgi_families[i];
        os<<LGIOutputStr(lgi_family)<<std::endl;
        os<<mcutils::FormatMatrix(lgi_expansions[i], ".8f")<<std::endl;
      }
  }


  void 
  WriteLGIExpansion(
    int Z, int N, int Nmax,
    const lgi::MultiplicityTaggedLGIVector& lgi_families,
    lsu3shell::OperatorBlocks&lgi_expansion,
    const std::string& filename
  )
    {
      // Open file
      std::ofstream os;
      os.open(filename);

      // Write header 
      WriteLGIExpansionHeader(Z,N,Nmax,os);
      
      //Write expansions 
      WriteLGIExpansion(lgi_families,lgi_expansion,os);

      os.close();
    }


  void ReadLGISet(
    const std::string& lgi_filename, 
    const HalfInt& Nsigma0,
    MultiplicityTaggedLGIVector& lgi_vector
    )
  {
    // open input file
    std::ifstream lgi_stream(lgi_filename);
    mcutils::StreamCheck(bool(lgi_stream),lgi_filename,"Failed to open LGI input file");

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
        line_stream >> Nex >> lambda >> mu
                    >> twice_Sp >> twice_Sn >> twice_S
                    >> count;
        
        mcutils::ParsingCheck(line_stream, line_count, line);
        
        // conversions
        u3::U3 sigma(Nsigma0+Nex,u3::SU3(lambda,mu));
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

 void ReadLGILookUpTable(std::vector<int>& lgi_full_space_lookup_table, int num_irrep_families)
  // Reading in and filling out table of lgi indices in basis and lgi indices in full space by 
  // which the seed files are labeled. 
    {
      lgi_full_space_lookup_table.resize(num_irrep_families);
      // open input file
      std::string filename="lgi_full_space_lookup_table.dat";
      std::ifstream stream(filename);
      mcutils::StreamCheck(bool(stream),filename,"Failed to open LGI input file");

      // scan input file
      std::string line;
      int line_count = 0;
      while ( std::getline(stream, line) )
        {
          // count line
          ++line_count;

          // set up for parsing
          // std::cout<<line<<std::endl;
          std::istringstream line_stream(line);

          // parse line
          int index, full_space_index;
          line_stream>>index>>full_space_index;
          mcutils::ParsingCheck(line_stream, line_count, line);

          lgi_full_space_lookup_table[index]=full_space_index;
        }
    }


}
