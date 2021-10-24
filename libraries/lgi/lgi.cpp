/****************************************************************
  lgi.cpp

  Anna E. McCoy and Mark A. Caprio
  University of Notre Dame and TRIUMF

  SPDX-License-Identifier: MIT
****************************************************************/
#include "lgi/lgi.h"

#include <fstream>
#include <iostream>

#include "mcutils/parsing.h"
#include "mcutils/eigen.h"
#include "fmt/format.h"


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

  std::string LGIFamilyStr(const MultiplicityTagged<lgi::LGI>& lgi_family)
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
        os<<LGIFamilyStr(lgi_count)<<std::endl;
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


MultiplicityTaggedLGIVector get_lgi_vector(
    const nuclide::NuclideType& nuclide, 
    const HalfInt& Nsigma0,
    const int& Nmax
  )
{
  // Initialize lgi_dimension with cmf U3SPN dimensions from lsu3shell basis
  std::map<u3shell::U3SPN, unsigned int>
  lgi_dimensions = lsu3shell::lsu3shell_cmf_basis_dimensions(nuclide,Nsigma0,Nmax);


  // Iterate through basis and identify LGI dimension by substracting
  // U(3) irreps obtained by laddering from lower grade LGI.
  for(const auto& [lgi,dimension] : lgi_dimensions)
    {
      HalfInt Sp(lgi.Sp()),Sn(lgi.Sn()), S(lgi.S());
      int Nn_max = Nmax - int(lgi.N() - Nsigma0);
      std::vector<u3::U3> raising_polynomial_labels = sp3r::RaisingPolynomialLabels(Nn_max);

      for(const u3::U3& n : raising_polynomial_labels)
        {
          if (n.N()==0)
            continue;

          MultiplicityTagged<u3::U3>::vector omegas_tagged = u3::KroneckerProduct(lgi.U3(), n);
          for(const auto& [omega,rho_max] : omegas_tagged)
            {
              u3shell::U3SPN omegaSpSnS(omega,Sp,Sn,S);
              lgi_dimensions[omegaSpSnS] -= rho_max*dimension;
            }
        }
    }

  //Create LGI vector used in SpNCCI basis construction
  MultiplicityTaggedLGIVector lgi_vector;
  for(const auto& [lgi_u3spn,dim] : lgi_dimensions)
    {
      int Nsex=int(lgi_u3spn.N()-Nsigma0);
      if(Nsex%2==Nmax%2)
      {
        lgi::LGI lgi(lgi_u3spn,Nsex);
        lgi_vector.emplace_back(lgi,dim);

      }
    }

  return lgi_vector;
}

}
