/****************************************************************
  dimensions_test.cpp

  Anna E. McCoy
  Institute for Nuclear Theory

  SPDX-License-Identifier: MIT

  11/2/21 (aem): Created.
****************************************************************/
#include <fstream>
#include <ostream>  

#include "fmt/format.h"
#include "lgi/lgi.h"
#include "u3ncsm/dimensions.h"
#include "lsu3shell/lsu3shell_basis.h"
#include "utilities/utilities.h"
#include "utilities/nuclide.h"
#include "sp3rlib/u3coef.h"

void test_lsu3shell_basis_generation()
  {
    nuclide::NuclideType nuclide={3,3};
    bool intrinsic=false;
    HalfInt Nsigma0=nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
    int Nmax=6;
    std::cout<<"Nsigma0 : "<<Nsigma0<<std::endl;
    std::map<u3shell::U3SPN, unsigned int> u3spn_dimensions
      =lsu3shell::generate_lsu3shell_basis_dimensions(nuclide,Nsigma0,Nmax);

    std::map<u3shell::U3SPN, unsigned int> u3spn_cmf_dimensions 
      =lsu3shell::lsu3shell_cmf_basis_dimensions(Nsigma0,Nmax,u3spn_dimensions);


    std::map<u3shell::U3SPN, unsigned int> u3spn_cmf_dimensions2 
      =lsu3shell::lsu3shell_cmf_basis_dimensions(nuclide,Nsigma0,Nmax);


    std::string spncci_root_dir=utils::get_spncci_project_root_dir();

    const std::string input_filename=fmt::format("{}/spncci/data/Z03-N03-Nmax14_u3s-dim.dat",spncci_root_dir);
    std::cout<<input_filename<<std::endl;
    std::map<u3shell::U3SPN, lsu3shell::Dimensions> u3spn_dimensions_test;
    lsu3shell::read_lsu3shell_basis_dimensions(
        input_filename,Nsigma0,Nmax,u3spn_dimensions_test
      );

    for(const auto& [u3spn,dimension]: u3spn_dimensions)
      {
        if(dimension !=u3spn_dimensions_test[u3spn].total)
          {
            fmt::print("Crap {}:  {}  {}\n",u3spn.Str(),dimension,u3spn_dimensions_test[u3spn].total);
            exit(EXIT_FAILURE);
          }
      }
    for(const auto& [labels,dim] : u3spn_dimensions)
      {
        fmt::print("{}:  {:6d}  {:6d}. {:6d}\n",labels.Str(),dim,u3spn_cmf_dimensions[labels],u3spn_cmf_dimensions2[labels]);
      }
  }

void test_lgi_generation()
{
  // Read in lgi vector from file
  std::string spncci_root_dir=utils::get_spncci_project_root_dir();
  std::string filename=fmt::format("{}/spncci/data/lgi_set/lgi_test_full.dat",spncci_root_dir);
  lgi::MultiplicityTaggedLGIVector lgi_vector;
  HalfInt Nsigma0=nuclide::Nsigma0ForNuclide({3,3});
  lgi::ReadLGISet(filename, Nsigma0,lgi_vector);

  //// Construct lgi vector directly from lsu3shell basis using counting arguments
  fmt::print("List of lgi's generated using lsu3shell basis constructors 6Li\n");
  nuclide::NuclideType nuclide({3,3});
  int Nmax=2;
  lgi::MultiplicityTaggedLGIVector lgi_vector2 = lgi::get_lgi_vector(nuclide, Nsigma0,Nmax);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //To compare lgi input file generated with old ordering, need to change ordering, 
  // so turn lgi_vector into map which does sorting and then write back to vector in new order.
  std::map<lgi::LGI,int> temp_map;
  for(const auto& [lgi,gamma_max] : lgi_vector)
    temp_map[lgi]=gamma_max;

  lgi_vector.resize(0);
  for(auto [lgi,gamma_max] : temp_map)
    lgi_vector.emplace_back(lgi,gamma_max);

  for(int i=0; i<lgi_vector.size(); ++i)
    {
      const auto&[lgi,gamma_max]=lgi_vector[i];
      const auto&[lgi2,gamma2_max]=lgi_vector2[i];
      if(!(lgi==lgi2) || !(gamma_max==gamma2_max))
        fmt::print("ERROR");

      fmt::print("{}  {:4d}   {}  {:4d}\n", lgi.Str(),gamma_max,lgi2.Str(),gamma2_max);
    }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

void compare_cmf_dimensions(const int Z, const int N, const int Nmax)
  {
    nuclide::NuclideType nuclide({Z,N});
    HalfInt Nsigma0=nuclide::Nsigma0ForNuclide(nuclide);
    std::map<u3shell::U3SPN, unsigned int> u3spn_dimensions
      =lsu3shell::generate_lsu3shell_basis_dimensions(nuclide,Nsigma0,Nmax);

    std::map<u3shell::U3SPN, unsigned int> u3spn_cmf_dimensions 
      =lsu3shell::lsu3shell_cmf_basis_dimensions(Nsigma0,Nmax,u3spn_dimensions);

    lgi::MultiplicityTaggedLGIVector lgi_vector = lgi::get_lgi_vector(nuclide,Nsigma0,Nmax);

    int num_gamma_max_zero = 0;
    for(const auto& [lgi,gamma_max] : lgi_vector)
      {
        if(gamma_max==0)
          {
            num_gamma_max_zero++;
            continue;
          }
        const auto u3spn = lgi.u3spn();
        fmt::print("{}:   {:6d}  {:6d}  {:6d}    {:6.3f}\n",
          u3spn.Str(),u3spn_dimensions[u3spn],
          u3spn_cmf_dimensions[u3spn],
          gamma_max,
          double(u3spn_cmf_dimensions[u3spn])/u3spn_dimensions[u3spn]
          );
      }
    fmt::print("Number of U3SPN subspaces which have no LGI: {}\n",num_gamma_max_zero); 
  }


void test_modified_branching_rule()
  {
    HalfInt Nsigma0(9,2);
    int Nmax=6;

    std::map<u3shell::U3SPN, unsigned int>
    lgi_dimensions = lsu3shell::lsu3shell_basis_dimensions({2,1},Nsigma0,Nmax);
    for(const auto& [lgi,dim] : lgi_dimensions)
      {
        fmt::print("{}  {}\n",lgi.Str(),dim);
      }


    for(const auto& [lgi,multiplicity] : lgi_dimensions)
      if(multiplicity>0)
        {
          const auto& s = lgi.U3();
          auto Nn_max = int(Nmax-(s.N()-Nsigma0));

          if(!sp3r::IsUnitary(s))
            {
              fmt::print(" {} is not unitary. multiplicity {}\n",s,multiplicity);
              continue;
            }

          std::cout<<s.Str()<<"  "<<multiplicity<<std::endl;
          // auto irrep = sp3r::Sp3RSpace(s,Nn_max);
          sp3r::Sp3RSpace irrep(s,Nn_max);
          // std::cout<<irrep.DebugStr()<<std::endl;

          for(const auto& subspace : irrep)
            {
              const u3::U3& omega = subspace.U3();
              if(s==omega)
                continue;
              u3shell::U3SPN labels(omega,lgi.Sp(),lgi.Sn(),lgi.S());
              // std::cout<<"   "<<labels.Str()<<std::endl;
              if(subspace.size()*multiplicity > lgi_dimensions[labels])
                {
                  std::cout<<labels.Str()<<"  "<<subspace.size()*multiplicity
                  <<"  "<<lgi_dimensions[labels]<<std::endl;
                }
              assert(subspace.size()*multiplicity<=lgi_dimensions[labels]);
              lgi_dimensions[labels]-=subspace.size()*multiplicity;
            }
        }
    std::cout<<"----------------------------------------"<<std::endl;
    for(const auto& [lgi,multiplicity] : lgi_dimensions)
      fmt::print("{}  {}\n",lgi.Str(),multiplicity);

    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<"----------------------------------------"<<std::endl;
  }

int main(int argc, char **argv)
{
  // test_lsu3shell_basis_generation();
  // test_lgi_generation();
  int Z=6;
  int N=6;
  // int Nmax=14;
  int Nmax=4;
  // compare_cmf_dimensions(N,Z,Nmax);
  u3::U3CoefInit(39);
  test_modified_branching_rule();

}//end main
