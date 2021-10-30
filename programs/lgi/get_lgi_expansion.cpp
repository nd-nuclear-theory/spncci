/****************************************************************
  get_lgi_expansion.cpp

  Code generating expansion of lgi subspaces in terms of 
  lsu3shell su3 reduced basis states 

  Inputs:
    Command line parameters:  "Syntax: Z N Nmax file_partition(optional)"
      Z: number of protons
      N: number of neutrons
      Nmax: many-body truncation parameter
      file_partition: Indicates whether or not expansion of all LGI written to single file or
        each LGI gets its own file.  Options are
          --"single" (default)
          --"individual"

    requires directory lsu3shell_rme containing files:
      * lsu3shell_basis.dat
      * Brel.rme
      * Nrel.rme

  Output if "single" file partition chosen: 
    "lgi_families.dat" : files containing lgi labels 
        Nex  lambda mu 2Sp 2Sn 2S gamma

    "lgi_expansions.dat" : file containing expansion of lgi's in lsu3shell basis


  Output if "individual" file partition chosen, creates input files for LSU3Shell code SU3RME_Sp3R_MPI
      Lgi_families.dat
        Contains Nex lambda mu 2Sp 2Sn 2S <lgi_expansions_filename>
      
      For each lgi, there is a separate binary file containing the lgi expansion called    
          lgi_expansions_Nex_lambda_mu_2Sp_2Sn_2S.dat 
        
      The file contains 
        binary_float_precision<int> : either 4 or 8.  Gives binary float precision
        nrow<uint32_t> : Corresponds to dim[U(3) subspace]
        ncol<uint32_t> : Corresponds to num LGI dim[Sp(3,R) subspace]
        Followed by vectors containing U(3) expansion of each LGI.
  

  Anna E. McCoy
  TRIUMF

  SPDX-License-Identifier: MIT

  3/30/19 (aem): Recreated. 
  10/10/19 (aem): Added option to write lgi expansion to separate files 
****************************************************************/
#include <fstream>  
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "sp3rlib/u3coef.h"
#include "lgi/lgi_solver.h"
#include "mcutils/eigen.h"
#include "lgi/lgi_unit_tensors.h"
#include "spncci/spncci_basis.h"
#include "u3shell/relative_operator.h"

#include "mcutils/io.h"

// Testing function 
namespace lgi{}//namespace

int main(int argc, char **argv)
{
  if(argc<4)
  {
    std::cout<<"Syntax: Z N Nmax <file_partition>(optional) "<<std::endl;
    std::cout<<"  options for <file_partition> are: "<<std::endl;
    std::cout<<"    single: expansions for all LGI written to a single file"<<std::endl;
    std::cout<<"    individual: expansions for each LGI subspace written separate file"<<std::endl;
    std::cout<<"    individual-text: expansions for each LGI subspace written separate file in text mode"<<std::endl;
    std::cout<<std::endl
      <<"Requires directory named lsu3shell_rme in cwd which containing files:"<<std::endl
      <<"   Brel.rme"<<std::endl
      <<"   Nrel.rme"<<std::endl
      <<"   lsu3shell_basis.dat"<<std::endl;

    std::exit(EXIT_FAILURE);
  }

  // nuclide 
  int Z=std::stoi(argv[1]);
  int N=std::stoi(argv[2]);

  // Basis parameters
  int Nmax=std::stoi(argv[3]);  
  
  // Optional file partition parameter
  std::string file_partition="single";
  if(argc==5)
    file_partition=argv[4];

  std::cout<<"argc "<<argc<<"  file partition "<<file_partition<<std::endl;
  if(file_partition!="single" and file_partition!="individual" and  file_partition!="individual-text")
    {
      std::cout<<"invalid value for optional parameter <file_partition>.  Options are:"<<std::endl;
      std::cout<<"   single: expansions for all LGI written to a single file"<<std::endl;
      std::cout<<"   individual: expansions for each LGI subspace written separate file"<<std::endl;
      std::cout<<"   individual-text: expansions for each LGI subspace written separate file in text mode"<<std::endl;
      std::exit(EXIT_FAILURE);
    }

  // zero tolerance 
  lgi::zero_tolerance=1e-6;
  
  // output mode
  lgi::binary_format_code = 1;
  lgi::binary_float_precision=8;

  std::array<int,2> nuclide; // proton and neutron numbers
  nuclide[0]=Z;
  nuclide[1]=N; 

  bool intrinsic=true;

  // su3rme output files
  std::string su3rme_filename_base="lsu3shell_rme";
  std::string lsu3shell_basis_filename=su3rme_filename_base+"/lsu3shell_basis.dat"; // Will need to include path to file

  // Generate Nsigma0 and N1v from nuclei and type 
  HalfInt Nsigma0 = nuclide::Nsigma0ForNuclide(nuclide,intrinsic);
  int N1v=nuclide::ValenceShellForNuclide(nuclide);
  
  // Operator parameters
  std::string Brel_filename=su3rme_filename_base+"/Brel.rme";
  std::string Nrel_filename=su3rme_filename_base+"/Nrel.rme";

  // Unit tensor parameters
  int J0=-1;
  int T0=-1;

  ////////////////////////////////////////////////////////////////
  // read lsu3shell basis
  ////////////////////////////////////////////////////////////////
  std::cout << "Read lsu3shell basis..." << std::endl;
  // read lsu3shell basis (regroup into U3SPN subspaces)
  lsu3shell::LSU3ShellBasisTable lsu3shell_basis_table;
  lsu3shell::U3SPNBasisLSU3Labels lsu3shell_basis_provenance;
  u3shell::SpaceU3SPN lsu3shell_space;
  lsu3shell::ReadLSU3ShellBasis(
      Nsigma0, lsu3shell_basis_filename,lsu3shell_basis_table,
      lsu3shell_basis_provenance,lsu3shell_space
    );

  ////////////////////////////////////////////////////////////////
  // solve for LGIs
  ////////////////////////////////////////////////////////////////
  std::cout << "Solve for LGIs..." << std::endl;
  lgi::MultiplicityTaggedLGIVector lgi_families;
  lsu3shell::OperatorBlocks lgi_expansions;
  std::vector<int> lsu3shell_index_lookup_table;

  lgi::GetLGIExpansion(
      lsu3shell_space,lsu3shell_basis_table,
      Brel_filename,Nrel_filename,Z+N, Nsigma0,
      lgi_families, lgi_expansions,
      lsu3shell_index_lookup_table
    );

  if(file_partition=="single")
    {
      std::string lgi_filename="lgi_families.dat";
      lgi::WriteLGILabels(lgi_families, lgi_filename);

      // std::cout<<"write expansion to file "<<std::endl;
      std::string lgi_expansion_filename="lgi_expansions.dat";
      lgi::WriteLGIExpansions(lgi_expansion_filename,lgi_expansions);
    }
  else if(file_partition=="individual")
    {

      std::string lgi_filename="lgi_families_individual.dat";
      std::ofstream lgi_stream;
      lgi_stream.open(lgi_filename);
    
      int Nex;
      u3::U3 sigma;
      HalfInt Sp,Sn,S;
      
      //For each lgi subspace, write expansion to file and record in logfile lgi_families_individial
      for(int lgi_index=0; lgi_index<lgi_families.size(); ++lgi_index)
        {
          
          //Get LGI labels and expansion
          auto& lgi_count=lgi_families[lgi_index];
          const lsu3shell::OperatorBlock& lgi_expansion=lgi_expansions[lgi_index];

          //Extract labels 
          std::tie(Nex,sigma,Sp,Sn,S)=lgi_count.irrep.Key();
          int count=lgi_count.tag;
          
          //Define output filename for expansion
          std::string expansion_filename=fmt::format("lgi_expansions_{:02d}_{}_{}_{}_{}_{}.dat",
            Nex,sigma.SU3().lambda(), sigma.SU3().mu(),
            TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S)
          );
          
          lgi_stream<<fmt::format("{}  {}  {}  {}  {}  {}  {}",
            Nex,sigma.SU3().lambda(), sigma.SU3().mu(),
            TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S),expansion_filename)
          <<std::endl;

          lgi::WriteLGIExpansions(expansion_filename,lgi_expansion);
        }
    }
  else
    {
      std::string lgi_filename="lgi_families_individual-text.dat";
      std::ofstream lgi_stream;
      lgi_stream.open(lgi_filename);
    
      int Nex;
      u3::U3 sigma;
      HalfInt Sp,Sn,S;
      
      //For each lgi subspace, write expansion to file and record in logfile lgi_families_individial
      for(int lgi_index=0; lgi_index<lgi_families.size(); ++lgi_index)
        {
          
          //Get LGI labels and expansion
          auto& lgi_count=lgi_families[lgi_index];
          const lsu3shell::OperatorBlock& lgi_expansion=lgi_expansions[lgi_index];

          //Extract labels 
          std::tie(Nex,sigma,Sp,Sn,S)=lgi_count.irrep.Key();
          int count=lgi_count.tag;
          
          //Define output filename for expansion
          std::string expansion_filename=fmt::format("lgi_expansions_{:02d}_{}_{}_{}_{}_{}_text.dat",
            Nex,sigma.SU3().lambda(), sigma.SU3().mu(),
            TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S)
          );
          
          lgi_stream<<fmt::format("{}  {}  {}  {}  {}  {}  {}",
            Nex,sigma.SU3().lambda(), sigma.SU3().mu(),
            TwiceValue(Sp),TwiceValue(Sn),TwiceValue(S),expansion_filename)
          <<std::endl;

          lgi::WriteLGIExpansionsText(expansion_filename, lgi_expansion);
        }
    }
    
}
